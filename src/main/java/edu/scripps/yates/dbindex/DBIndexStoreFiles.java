package edu.scripps.yates.dbindex;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.ws.rs.NotSupportedException;

import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.utilities.fasta.Fasta;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.fasta.dbindex.MassRange;
import edu.scripps.yates.utilities.fasta.dbindex.ResidueInfo;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;

/**
 *
 * Fasta sequence index storing peptides by mass using files, a single file per
 * sequence mass (there could be many files per ppm)
 *
 * Could be more suitable for huge databases, but not as fasta as full SQLite
 * based index, it may scale better for huge indexes.
 *
 * @author Adam
 */
public class DBIndexStoreFiles implements DBIndexStore {

	private ProteinCache protCache;
	private String indexDir;
	private static final Logger logger = Logger.getLogger(DBIndexStoreFiles.class.getName());
	private long indexed = -1;
	private boolean inited = false;
	private File[] indexFiles;
	private boolean indexExists;
	private final String INDEX_DIR_SUFFIX = ".idx2";
	private TIntObjectHashMap<Void> indexedMassKeys; // mass * MAX_FILES_PER_PPM
	private static double MAX_FILES_PER_PPM = 0.1f; // 10 ppm per file
	// caches open handles, maps mass key to a handle
	private final TIntObjectHashMap<OutputStream> writerCache = new TIntObjectHashMap<OutputStream>();
	private final static Charset ENCODING = Charset.forName("UTF-8");
	private final byte[] SEP;
	private final byte[] NLSEP;
	private static final int MAX_MASS = (int) Constants.MAX_PRECURSOR_MASS;
	// private static final int MAX_MASS = (int) + 400;
	private static final int READ_BUFF_SIZE = 10 * 1024 * 1024;

	DBIndexStoreFiles() {
		SEP = " ".getBytes(ENCODING);
		NLSEP = "\n".getBytes(ENCODING);
	}

	@Override
	public void init(String databaseID) throws DBIndexStoreException {
		indexDir = databaseID + INDEX_DIR_SUFFIX;

		indexExists = indexExists();

		if (indexExists) {
			indexedMassKeys = new TIntObjectHashMap<Void>();
			// populate mass keys from files
			for (final File f : indexFiles) {
				final String fName = f.getName();
				final String mass = stripExtension(fName);
				try {
					indexedMassKeys.put(Integer.parseInt(mass), null);
				} catch (final NumberFormatException e) {
					logger.log(Level.SEVERE, "Unexpected file in index: " + fName, e);
					throw new DBIndexStoreException("Unexpected file in index: " + fName);
				}
			}

		} else {
			final File indexDirF = new File(indexDir);
			indexDirF.mkdir();
			if (!indexDirF.exists()) {
				logger.log(Level.SEVERE, "Error creating index dir: " + indexDir);
				throw new DBIndexStoreException("Error creating index dir: " + indexDir);
			}
		}

		inited = true;
	}

	@Override
	public boolean indexExists() throws DBIndexStoreException {

		final File indexDirF = new File(indexDir);
		if (!indexDirF.exists()) {
			return false;
		}

		indexFiles = indexDirF.listFiles();
		return indexFiles.length > 0;

	}

	@Override
	public void startAddSeq() throws DBIndexStoreException {
		if (indexExists) {
			throw new DBIndexStoreException("Index already exists, overwriting not supported");
		}

		indexed = 0;
	}

	@Override
	public void stopAddSeq() throws DBIndexStoreException {
		logger.log(Level.INFO, "Clearing writer cache");

		for (final OutputStream fw : writerCache.valueCollection()) {
			try {
				fw.flush();
				fw.close();
			} catch (final IOException e) {
				logger.log(Level.WARNING, "Error closing writer", e);
			}
		}
		writerCache.clear();

		// TODO merge sequences (optional)

		// TODO sort every file by mass (use sort.exe)

		indexExists = true;

		System.out.println("Total sequences added to index: " + Long.toString(indexed));

	}

	@Override
	public FilterResult filterSequence(double precMass, String sequence) {
		// no filtering, we add every sequence
		return FilterResult.INCLUDE;
	}

	private int getMassKey(double mass) {
		return (int) (mass * MAX_FILES_PER_PPM);
	}

	@Override
	public void addSequence(double precMass, int sequenceOffset, int sequenceLen, String sequence, String resLeft,
			String resRight, long proteinId) throws DBIndexStoreException {
		if (indexed % 10000000 == 0) {
			System.out.println(Long.toString(indexed));
		}
		++indexed;

		final int massKey = getMassKey(precMass);

		// get writer from cache
		OutputStream indexWriter = writerCache.get(massKey);
		if (indexWriter == null) {
			try {
				indexWriter = new BufferedOutputStream(
						new FileOutputStream(indexDir + File.separator + massKey + INDEX_DIR_SUFFIX));
			} catch (final IOException ex) {
				logger.log(Level.SEVERE, "Cannot open index file writer for mass key: " + massKey, ex);
				throw new DBIndexStoreException("Cannot open index file writer for mass key: " + massKey, ex);
			}
			// if (massKey < 500) {
			writerCache.put(massKey, indexWriter);
			// }
		}

		try {
			indexWriter.write(Double.toString(precMass).getBytes(ENCODING));
			indexWriter.write(SEP);
			indexWriter.write(Integer.toString(sequenceOffset).getBytes(ENCODING));
			indexWriter.write(SEP);
			indexWriter.write(Integer.toString(sequenceLen).getBytes(ENCODING));
			indexWriter.write(SEP);
			indexWriter.write(Integer.toString((int) proteinId).getBytes(ENCODING));
			indexWriter.write(NLSEP);

		} catch (final IOException ex) {
			logger.log(Level.SEVERE, "Cannot write to index file writer for mass key: " + massKey, ex);
			throw new DBIndexStoreException("Cannot write to index file writer for mass key: " + massKey, ex);
		} finally {
			// if (massKey >= 500) {
			// try {
			// indexWriter.close();
			// } catch (IOException ex) {
			// Logger.getLogger(DBIndexStoreFiles.class.getName()).log(Level.SEVERE,
			// null, ex);
			// }
			// }
		}

	}

	private int[] getFileRangeForMassRange(double minMass, double maxMass) {
		if (maxMass > MAX_MASS) {
			logger.log(Level.SEVERE, "Trying to get a sequence for more than supported mass: " + maxMass);
		}

		final int[] ret = new int[2];
		ret[0] = getMassKey(minMass);
		ret[1] = getMassKey(maxMass);

		return ret;
	}

	@Override
	public List<IndexedSequence> getSequences(double precMass, double tolerance) throws DBIndexStoreException {
		if (!indexExists) {
			throw new DBIndexStoreException("Can't get sequences, index does not exists");
		}

		double massLow = precMass - tolerance;
		if (massLow < 0f) {
			massLow = 0f;
		}
		final double massHigh = precMass + tolerance;

		return getSequencesMass(massLow, massHigh);
	}

	private List<IndexedSequence> getSequencesMass(double massLow, double massHigh) throws DBIndexStoreException {

		final int[] fileRange = getFileRangeForMassRange(massLow, massHigh);

		final List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

		// logger.log(Level.INFO, "Mass range query: " + massLow + " - " +
		// massHigh + ", files: " + fileRange[0] + " - " + fileRange[1]);

		// go for every file
		for (int massKey = fileRange[0]; massKey <= fileRange[1]; ++massKey) {

			BufferedReader reader = null;

			// check map that stores mass keys if entry exists in index
			if (!indexedMassKeys.containsKey(massKey)) {
				// no sequences matching for this mass, skip the masskey
				continue;
			}

			final File subIndexFile = new File(indexDir + File.separator + massKey + INDEX_DIR_SUFFIX);

			try {
				reader = new BufferedReader(new FileReader(subIndexFile), READ_BUFF_SIZE);
				// reader = new RandomAccessFile(subIndexFile, "r");
				// readerCache.put(massKey, reader);
			} catch (final FileNotFoundException ex) {
				logger.log(Level.SEVERE, "Error initializing reader for mass key: " + massKey, ex);
				throw new DBIndexStoreException("Error initializing reader for mass key: " + massKey, ex);
			}

			// temp sequences, before merge by unique sequence, for a single
			// file
			final Map<String, List<IndexedSeqInternal>> temp = new THashMap<String, List<IndexedSeqInternal>>();

			// read all sequences, get the ones with masses in the range

			// TODO when we have sorted sequences, we can bail out quicker when
			// passed max mass
			// and we can binary search to get to the start sequence (future
			// optimization)

			// try {
			// reader.seek(0);
			// } catch (IOException ex) {
			// logger.log(Level.SEVERE, "Error resetting the stream, massKey",
			// ex);
			// throw new
			// DBIndexStoreException("Error resetting the stream, massKey " +
			// massKey, ex);
			// }

			try {
				String sequenceLine = reader.readLine();
				for (; sequenceLine != null; sequenceLine = reader.readLine()) {
					// parse mass first to qualify
					final int firstSep = sequenceLine.indexOf(' ');
					if (firstSep == -1) {
						logger.log(Level.SEVERE, "Invalid sequence line, expecting 4 tokens: " + sequenceLine);
						throw new DBIndexStoreException("Invalid sequence line, expecting 4 tokens: " + sequenceLine);
					}
					final double seqMass = Double.parseDouble(sequenceLine.substring(0, firstSep));
					if (!(seqMass >= massLow && seqMass <= massHigh)) {
						// skip the sequence
						continue;
						// TODO when sorted, we can skip entire file if passed
						// massHigh

					}

					final String restLine = sequenceLine.substring(firstSep + 1);

					final String[] toks = restLine.split(" ");

					if (toks.length != 3) {
						logger.log(Level.SEVERE, "Invalid sequence line, expecting 4 tokens: " + sequenceLine);
						throw new DBIndexStoreException("Invalid sequence line, expecting 4 tokens: " + sequenceLine);
					}

					// within mass range
					final int offset = Integer.parseInt(toks[0]);
					final int length = Integer.parseInt(toks[1]);

					// System.out.println("----" + toks[2]);
					final String[] indexArr = toks[2].split(",");
					final int proteinId = Integer.parseInt(indexArr[0]);

					// construct temp sequence, will merge after the loop
					// TODO when sequences are unique in the index, we don't
					// need to merge at search time, add merge-at-index later
					final String peptideSequence = protCache.getPeptideSequence(proteinId, offset, length);
					// System.out.println("pep: " + peptideSequence);

					final IndexedSeqInternal tempSequence = new IndexedSeqInternal(seqMass, offset, length, proteinId,
							null);
					tempSequence.setProteinIdArr(indexArr);

					List<IndexedSeqInternal> sequences = temp.get(peptideSequence);
					if (sequences == null) {
						sequences = new ArrayList<IndexedSeqInternal>();
						temp.put(peptideSequence, sequences);
					}
					sequences.add(tempSequence);
				} // for every line

			} catch (final IOException e) {
				logger.log(Level.SEVERE, "Error reading sequence line for massKey: " + massKey, e);
				throw new DBIndexStoreException("Error reading sequence line for massKey " + massKey, e);
			} catch (final NumberFormatException e) {
				logger.log(Level.SEVERE, "Invalid format in sequence line for massKey: " + massKey, e);
				throw new DBIndexStoreException("Invalid format in sequence line for massKey: " + massKey, e);
			} finally {
				try {
					reader.close();
				} catch (final IOException ex) {
					logger.log(Level.SEVERE, "Error closing reader for mass key: " + massKey, ex);
				}
			}

			// group the same peptides from many proteins into single peptide
			// with protein id list
			for (final String pepSeqKey : temp.keySet()) {
				// for each sequence str

				final List<IndexedSeqInternal> sequences = temp.get(pepSeqKey);

				final List<Integer> proteinIds = new ArrayList<Integer>();
				for (final IndexedSeqInternal tempSeq : sequences) {
					proteinIds.add(tempSeq.proteinId);
				}

				// make sure the 1st protein id is that of the first sequence
				final IndexedSeqInternal firstSeq = sequences.get(0);

				final IndexedSequence firstSequence = new IndexedSequence(0, firstSeq.mass, pepSeqKey, "", "");
				firstSequence.setProteinIds(proteinIds);
				// set residues, NOTE only from first sequence

				final String proteinSequence = protCache.getProteinSequence(firstSeq.proteinId);
				final ResidueInfo residues = Util.getResidues(null, firstSeq.offset, firstSeq.length, proteinSequence);
				firstSequence.setResidues(residues);

				ret.add(firstSequence);

			} // end of this unique merged sequence

		} // end for every file

		return ret;
	}

	@Override
	public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {
		if (!indexExists) {
			throw new DBIndexStoreException("Can't get sequences, index does not exists");
		}

		final int size = ranges.size();
		if (size == 0) {
			return Collections.<IndexedSequence>emptyList();
		} else if (size == 1) {
			final MassRange range1 = ranges.get(0);
			return getSequences(range1.getPrecMass(), range1.getTolerance());
		}

		final List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

		// merge mass ranges to reduce dups before query time
		final ArrayList<Interval> intervals = new ArrayList<Interval>();

		for (final MassRange range : ranges) {
			final Interval ith = Interval.massRangeToInterval(range);
			intervals.add(ith);
		}

		final List<Interval> merged = MergeIntervals.mergeIntervals(intervals);

		// submit merged intervals
		for (final Interval interval : merged) {
			ret.addAll(getSequencesMass(interval.getStart(), interval.getEnd()));
		}

		// for (MassRange range : ranges) {
		// ret.addAll(getSequences(range.getPrecMass(), range.getTolerance()));
		// }

		return ret;
	}

	@Override
	public long addProteinDef(long num, String accession, String protSequence) throws DBIndexStoreException {
		// using protein cache entirely
		// nothing here
		return num;
	}

	@Override
	public void setProteinCache(ProteinCache proteinCache) {
		protCache = proteinCache;
	}

	@Override
	public boolean supportsProteinCache() {
		return true;
	}

	@Override
	public List<IndexedProtein> getProteins(IndexedSequence sequence) throws DBIndexStoreException {
		// get prot id from index and prot ids from cache
		final List<IndexedProtein> ret = new ArrayList<IndexedProtein>();

		for (final Integer protId : sequence.getProteinIds()) {
			final IndexedProtein ip = new IndexedProtein(protCache.getProteinDef(protId), protId);
			ret.add(ip);
		}

		return ret;
	}

	@Override
	public long getNumberSequences() throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer not inited, can't check num sequences");
		}
		if (!indexExists) {
			return 0;
		}

		if (indexed == -1) {
			logger.log(Level.INFO, "Counting sequences in index");
			indexed = 0;
			for (final File f : indexFiles) {
				final String absPath = f.getAbsolutePath();
				try {
					indexed += Util.countLines(absPath);
				} catch (final IOException ex) {
					logger.log(Level.SEVERE, "Error counting sequences", ex);
					throw new DBIndexStoreException("Error counting sequences", ex);

				}
			}
			logger.log(Level.INFO, "Done counting sequences in index");

		}

		return indexed;
	}

	@Override
	public ResidueInfo getResidues(IndexedSequence peptideSequence, IndexedProtein protein)
			throws DBIndexStoreException {
		// in this implementation residues are already stored in the index
		return new ResidueInfo(peptideSequence.getResLeft(), peptideSequence.getResLeft());
	}

	private static String stripExtension(String fileName) {
		final int dot = fileName.indexOf('.');
		if (dot == -1) {
			return fileName;
		}

		return fileName.substring(0, dot);

	}

	/**
	 * TODO this method should be private and invoked only internally in this class
	 * after initial indexing is done
	 *
	 * @param sparam saerch params
	 * @throws DBIndexerException
	 * @throws IOException
	 */
	protected void merge(DBIndexSearchParams sparam) throws DBIndexerException, IOException {
		// reset prot number

		final String indexName = sparam.getDatabaseName();

		final File dbFile = new File(indexName);
		if (!dbFile.exists() || !dbFile.canRead()) {
			throw new DBIndexerException("Cannot read (and index) the Fasta database file: " + indexName);
		}

		// get param specific name for the inde0x
		// indexName = this.getFullIndexFileName();

		Iterator<Fasta> itr = null;
		try {
			itr = IndexUtil.getFastaReader(sparam).getFastas();
		} catch (final IOException ex) {
			logger.log(Level.SEVERE, "Could not set protein cache", ex);
			return;
		}

		while (itr.hasNext()) {
			final Fasta fasta = itr.next();
			protCache.addProtein(fasta.getSequestLikeAccession(), fasta.getSequence());
		}

		final String[] flist = dbFile.list();
		final String path = indexName + File.separator;

		for (final String f : flist) {
			// System.out.println(path + f);
			final Hashtable<String, IndexedPeptideModel> pht = new Hashtable<String, IndexedPeptideModel>();

			final BufferedReader br = new BufferedReader(new FileReader(path + f));
			String eachLine = null;
			while ((eachLine = br.readLine()) != null) {

				// System.out.println("=== " + eachLine);
				final String[] arr = eachLine.split(" ");
				final int pid = Integer.parseInt(arr[3]);
				if (pid < 0) {
					br.close();
					throw new IllegalArgumentException("Probably integer overflow!");
				}
				final double mass = Double.parseDouble(arr[0]);
				final int startIndex = Integer.parseInt(arr[1]);
				final int offsetIndex = Integer.parseInt(arr[2]);

				final String seq = protCache.getPeptideSequence(pid, startIndex, offsetIndex);
				final IndexedPeptideModel pm = new IndexedPeptideModel(mass, arr[3], startIndex, offsetIndex);

				final IndexedPeptideModel tpm = pht.get(seq);
				if (null == tpm) {
					pht.put(seq, pm);
				} else {
					// System.out.println("merge " + seq);
					tpm.merge(pm);
					// pht.put(seq, tpm);
				}

				// System.out.println("======" + pid + " " + startIndex + " " +
				// offsetIndex);
				// System.out.println(seq);
			}

			// System.out.println("======" + pht.size());

			final List<IndexedPeptideModel> list = new ArrayList<IndexedPeptideModel>(pht.values());

			// for(int i=0;i<10;i++)
			// System.out.println(list.get(i).getMass() + " " +
			// list.get(i).getProteinIds());

			Collections.sort(list);
			final File oldFile = new File(path + f);
			oldFile.delete();

			final BufferedWriter fileWriter = new BufferedWriter(new FileWriter(path + f), 1024 * 128);
			for (final Iterator<IndexedPeptideModel> pitr = list.iterator(); pitr.hasNext();) {
				final IndexedPeptideModel pm = pitr.next();
				fileWriter.write(String.valueOf(pm.getMass()));
				fileWriter.write(" ");
				fileWriter.write(String.valueOf(pm.getStartPos()));
				fileWriter.write(" ");
				fileWriter.write(String.valueOf(pm.getOffset()));
				fileWriter.write(" ");
				fileWriter.write(pm.getProteinIds());
				fileWriter.write("\n");
			}
			br.close();
			fileWriter.close();

			// System.out.println("======" + pht.size() + " " + list.size());

			/*
			 * for(int i=0;i<10;i++) { String[] pid =
			 * list.get(i).getProteinIds().split(","); if(pid.length>1)
			 * System.out.println(list.get(i).getMass() + " " +
			 * list.get(i).getProteinIds()); }
			 */

			// System.exit(0);

		}

		// System.out.println("aaaaaaaaaaaaaaaa" + protCache);
		// System.out.println("aaaaaaaaaaaaaaaa" +
		// protCache.getProteinSequence(5));

	}

	/**
	 * test driver
	 *
	 * @param args
	 */
	public static void main(String[] args) {
		DBIndexStore store = new DBIndexStoreFiles();

		final String protDef1 = "4R79.2 CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";
		final String protDef2 = "Reverse_4R79.2  CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";

		try {
			logger.log(Level.INFO, "Testing init");
			final File idx = new File("./test.fasta.idx2");
			if (idx.exists() && idx.isDirectory()) {
				for (final File f : idx.listFiles()) {
					f.delete();
				}
				idx.delete();
			}

			final ProteinCache protCache = new ProteinCache();
			store.init("./test.fasta");
			store.setProteinCache(protCache);

			// test a transaction
			logger.log(Level.INFO, "Testing adds 1");
			store.startAddSeq();

			final long protId = store.addProteinDef(0, protDef1, "ABCDEFGHIJKL");
			protCache.addProtein("4R79.2", "ABCDEFGHIJKL");

			store.addSequence(1f, 0, 1, "A", null, null, protId);
			store.addSequence(2f, 0, 2, "AB", null, null, protId);
			store.addSequence(3f, 0, 3, "ABC", null, null, protId);
			store.addSequence(4f, 0, 4, "ABCD", null, null, protId);
			store.addSequence(6000.42323f, 0, 5, "ABCDE", null, null, protId);
			store.addSequence(6999.42323f, 0, 6, "ABCDEZ", null, null, protId);
			store.addSequence(3, 6, 3, "GHI", null, null, protId);
			// store.addSequence(3, 6, 3, "GHI", null, null, protId);
			// store.stopAddSeq();
		} catch (final DBIndexStoreException ex) {
			logger.log(Level.SEVERE, null, ex);
		}

		try {
			// test more adds with the same statement in a new transaction
			logger.log(Level.INFO, "Testing adds 2");
			// store.startAddSeq();

			final long protId = store.addProteinDef(1, protDef2, "GHIJKLMNOPR");
			final ProteinCache protCache = new ProteinCache();
			protCache.addProtein("Reverse_4R79.2", "GHIJKLMN");

			store.addSequence(3, 1, 3, "HIJ", null, null, protId);
			store.addSequence(5, 2, 5, "IJKLM", null, null, protId);
			// store.addSequence(7, 2, 7, "IJKLMNO", null, null, protId);
			// store.addSequence(7.1f, 2, 8, "IJKLMNOO", null, null, protId);
			store.addSequence(3, 0, 3, "GHI", null, null, protId); // test dup
			store.addSequence(3, 0, 3, "GHI", null, null, protId); // test dup
			store.stopAddSeq();
		} catch (final DBIndexStoreException ex) {
			logger.log(Level.SEVERE, null, ex);
		}

		try {
			// test gets
			logger.log(Level.INFO, "Testing gets 1");
			logger.log(Level.INFO, "NUM SEQUENCES: " + store.getNumberSequences());
			final List<IndexedSequence> res1 = store.getSequences(10, 8.9f);
			for (final IndexedSequence seq : res1) {
				System.out.println(seq);
				final List<IndexedProtein> proteins = store.getProteins(seq);
				for (final IndexedProtein protein : proteins) {
					final ResidueInfo res = store.getResidues(seq, protein);
					System.out.println(protein);
					System.out.println(res);
				}
			}
		} catch (final DBIndexStoreException ex) {
			logger.log(Level.SEVERE, null, ex);
		}

		try {
			// test gets
			System.out.println("Testing range gets 1");
			final MassRange r1 = new MassRange(2, 1f);
			final MassRange r2 = new MassRange(6, 1f);
			final MassRange r3 = new MassRange(6, 1f); // test redundant
			final MassRange r4 = new MassRange(6, 1.2f); // test overlapping
			final List<MassRange> ranges = new ArrayList<MassRange>();
			ranges.add(r2);
			ranges.add(r1);
			ranges.add(r3);
			ranges.add(r4);
			final List<IndexedSequence> res1 = store.getSequences(ranges);
			for (final IndexedSequence seq : res1) {
				System.out.println(seq);
				final List<IndexedProtein> proteins = store.getProteins(seq);
				for (final IndexedProtein protein : proteins) {
					final ResidueInfo res = store.getResidues(seq, protein);
					System.out.println(protein);
					System.out.println(res);
				}
			}
		} catch (final DBIndexStoreException ex) {
			logger.log(Level.SEVERE, null, ex);
		}

		if (true) {
			return;
		}

		try {
			store = new DBIndexStoreFiles();
			store.init("EBI-IPI_Human_IPI_Human_3_85_06-29-2011_reversed.fasta");
			for (int i = 0; i < 1000; ++i) {
				logger.log(Level.INFO, "Testing getting sequence " + i);
				final List<IndexedSequence> res = store.getSequences(i, 3f);
				if (res.size() > 0) {
					final IndexedSequence seq = res.get(0);
					logger.log(Level.INFO, "Got " + res.size() + " sequences, first sequence: " + seq);
					final List<IndexedProtein> proteins = store.getProteins(seq);
					logger.log(Level.INFO, "Got " + res.size() + " sequences, first sequence's proteins: ");
					for (final IndexedProtein protein : proteins) {
						System.out.println(protein);
					}
				}

			}

		} catch (final DBIndexStoreException ex) {
			logger.log(Level.SEVERE, null, ex);
		}

		try {
			logger.log(Level.INFO, "Testing getting multiple sequence proteins ");
			store = new DBIndexStoreFiles();
			store.init("test_02.fasta");

			final List<IndexedSequence> res = store.getSequences(1000, 3000);
			for (final IndexedSequence seq : res) {
				System.out.println(seq);
				final List<IndexedProtein> prots = store.getProteins(seq);
				for (final IndexedProtein prot : prots) {
					System.out.println("\t" + prot);
				}
			}

		} catch (final DBIndexStoreException ex) {
			logger.log(Level.SEVERE, null, ex);
		}

	}

	@Override
	public void lastBuffertoDatabase() {
		throw new UnsupportedOperationException("Not supported yet."); // To
																		// change
																		// body
																		// of
																		// generated
																		// methods,
																		// choose
																		// Tools
																		// |
																		// Templates.
	}

	@Override
	public Iterator<IndexedSequence> getSequencesIterator(List<MassRange> ranges) throws DBIndexStoreException {
		return getSequences(ranges).iterator();
	}

	@Override
	public List<Integer> getEntryKeys() throws DBIndexStoreException {
		throw new NotSupportedException("Method not implemented");
	}
}
