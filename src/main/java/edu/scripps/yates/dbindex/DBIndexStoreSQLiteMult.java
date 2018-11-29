package edu.scripps.yates.dbindex;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.fasta.dbindex.MassRange;
import edu.scripps.yates.utilities.fasta.dbindex.ResidueInfo;
import gnu.trove.set.hash.THashSet;

/**
 * Using multiple databases to divide index into multiple manageable chunks
 *
 * @author Adam
 *
 */
public class DBIndexStoreSQLiteMult implements DBIndexStore {

	protected long totalSeqCount = 0;
	protected boolean inited = false;
	protected ProteinCache proteinCache;
	protected String dbPathBase = null;
	private static final org.apache.log4j.Logger logger = org.apache.log4j.Logger
			.getLogger(DBIndexStoreSQLiteMult.class);
	protected final DBIndexStoreSQLiteByteIndexMerge[] buckets;
	protected static final String IDX_SUFFIX = ".idx";
	// we divide into buckets based on this
	protected static final int MAX_MASS = (int) Constants.MAX_PRECURSOR_MASS;
	// adjust number of buckets to make sure no single bucket ends up with more
	// than 100million sequences

	// 1 ppm per row
	// private static final int PPM_PER_ENTRY = 1;

	// offset (4), length (4),
	// protein id (4)

	protected final DBIndexSearchParams sparam;
	protected boolean inMemoryIndex = false;

	public DBIndexStoreSQLiteMult(edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams sparam, boolean inMemoryIndex,
			DBIndexStoreSQLiteByteIndexMerge[] buckets) {
		inited = false;
		this.sparam = sparam;
		totalSeqCount = 0;

		this.inMemoryIndex = inMemoryIndex;

		Constants.NUM_BUCKETS = sparam.getIndexFactor();
		Constants.BUCKET_MASS_RANGE = MAX_MASS / Constants.NUM_BUCKETS;

		this.buckets = buckets;

		if (inMemoryIndex) {
			// throw new
			// IllegalArgumentException("In-memory index not supported");
		}
	}

	public DBIndexStoreSQLiteMult(edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams sparam, boolean inMemoryIndex) {
		this(sparam, inMemoryIndex, new DBIndexStoreSQLiteByteIndexMerge[sparam.getIndexFactor()]);
	}

	@Override
	public boolean indexExists() throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Not intialized");
		}

		final File indexDir = new File(dbPathBase);
		if (!indexDir.exists()) {
			return false;
		}

		final File[] indexFiles = indexDir.listFiles();
		if (indexFiles.length == 0) {
			return false;
		}

		return hasSequences();

	}

	@Override
	public void init(String databaseID) throws DBIndexStoreException {
		if (databaseID == null || databaseID.equals("")) {
			throw new DBIndexStoreException("Index path is missing, cannot initialize the indexer.");
		}

		if (inited) {
			throw new DBIndexStoreException("Already intialized");
		}

		final File dbBase = new File(databaseID).getAbsoluteFile();
		String baseName = dbBase.getName();
		if (!baseName.endsWith(IDX_SUFFIX)) {
			baseName = baseName + IDX_SUFFIX;
		}
		final File parentDir = dbBase.getParentFile();

		final File indexDir = new File(parentDir.getAbsolutePath() + File.separator + baseName);
		if (indexDir.exists() && !indexDir.isDirectory()) {
			logger.info("Trying to delete old index file: " + indexDir.getAbsolutePath());
			indexDir.delete();

		}

		try {
			indexDir.mkdir();
		} catch (final SecurityException e) {
		}

		dbPathBase = indexDir.getAbsolutePath();

		logger.info("Using database index dir: " + dbPathBase);

		if (!indexDir.exists()) {
			throw new DBIndexStoreException(
					"Index dir does not exist: " + indexDir + ", cannot initialize the indexer.");
		}

		// initialize buckets
		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			// buckets[i] = new DBIndexStoreSQLiteByte(i); //version with merge
			// during search
			buckets[i] = new DBIndexStoreSQLiteByteIndexMerge(sparam, inMemoryIndex, i, proteinCache); // version
			// with
			// merge
			// during
			// index
			final String bucketId = dbPathBase + File.separator + i + IDX_SUFFIX;
			try {
				buckets[i].init(bucketId);
			} catch (final DBIndexStoreException e) {
				logger.error("Error initializing bucket " + i, e);
				throw e;
			}
		}

		inited = true;

	}

	@Override
	public void startAddSeq() throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		logger.info("Starting adding sequences");

		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			buckets[i].startAddSeq();
		}
	}

	@Override
	public void stopAddSeq() throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		logger.info("Finalizing adding sequences");

		logger.info("Will flush sequences from cache, merge sequences, index and compact the index database.");

		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			buckets[i].stopAddSeq();
		}

		logger.info("Total sequences added: " + totalSeqCount);
	}

	@Override
	public long getNumberSequences() throws DBIndexStoreException {
		// NOTE, this only returns total number of mass entries, not really
		// useful anymore
		// we don't know num of sequenes in index at search time, because we
		// would have to readit all
		long total = 0;
		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			total += buckets[i].getNumberSequences();
		}
		return total;
	}

	@Override
	public List<Integer> getEntryKeys() throws DBIndexStoreException {
		final List<Integer> ret = new ArrayList<Integer>();

		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			ret.addAll(buckets[i].getEntryKeys());
		}
		return ret;
	}

	public boolean hasSequences() throws DBIndexStoreException {

		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			final long total = buckets[i].getNumberSequences();
			if (total > 0) {
				return true;
			}
		}
		return false;
	}

	private int getBucketForMass(double precMass) {
		return (int) precMass / Constants.BUCKET_MASS_RANGE;
	}

	/**
	 * get min, max bucket inclusive for the orig. mass range
	 *
	 * @param minMass
	 * @param maxMass
	 * @return
	 */
	private int[] getBucketsForMassRange(double minMass, double maxMass) {
		if (maxMass > MAX_MASS) {
			logger.error("Trying to get a sequence for more than supported mass: " + maxMass);
		}

		final int[] ret = new int[2];
		ret[0] = (int) minMass / Constants.BUCKET_MASS_RANGE;
		ret[1] = (int) maxMass / Constants.BUCKET_MASS_RANGE;

		// if (ret[0] != ret[1]) {
		// checking how often we are between buckets
		// logger.info("Between buckets: masses: " + minMass + ", " +
		// maxMass + ", buckets: "
		// + ret[0] + ", " + ret[1]);
		// }
		return ret;
	}

	@Override
	public FilterResult filterSequence(double precMass, String sequence) {
		// check if the sequence has the mandatory AAs
		final char[] mandatoryInternalAAs = sparam.getMandatoryInternalAAs();
		if (mandatoryInternalAAs != null && mandatoryInternalAAs.length > 0) {
			// to not repeat
			final Set<String> aas = new THashSet<String>();
			for (final char mandatoryInternalAAChar : mandatoryInternalAAs) {
				final String mandatoryInternalAA = String.valueOf(mandatoryInternalAAChar);
				if (!aas.contains(mandatoryInternalAA)) {
					aas.add(mandatoryInternalAA);
					// exclude the last AA, which is the cleavage site
					final String sequenceTMP = sequence.substring(0, sequence.length() - 1);
					if (sequenceTMP.contains(mandatoryInternalAA)) {
						return FilterResult.INCLUDE;
					}
				}
			}
			return FilterResult.SKIP;
		}
		if (sparam.getMaxPrecursorMass() < precMass || sparam.getMinPrecursorMass() > precMass) {
			return FilterResult.SKIP;
		}
		return FilterResult.INCLUDE;
	}

	@Override
	public void addSequence(double precMass, int seqOffset, int seqLength, String sequence, String resLeft,
			String resRight, long proteinId) throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		totalSeqCount++;
		if (totalSeqCount % 100000 == 0) {
			logger.info("Total peptides indexed so far: " + totalSeqCount);
		}

		final int bucket = getBucketForMass(precMass);
		if (bucket > Constants.NUM_BUCKETS - 1) {
			logger.error("Cannot add to index, unsupported precursor mass: " + precMass + " for sequence: " + sequence
					+ ". Max supported mass is: " + MAX_MASS);
			return;
		}
		// logger.info("Inserting into bucket: " + bucket);
		buckets[bucket].addSequence(precMass, seqOffset, seqLength, sequence, resLeft, resRight, proteinId);

	}

	@Override
	public ResidueInfo getResidues(IndexedSequence peptideSequence, IndexedProtein protein)
			throws DBIndexStoreException {
		// get the true residues for this peptide
		final String proteinSequence = proteinCache.getProteinSequence(protein.getId());

		int seqOffset = peptideSequence.getSequenceOffset();
		if (seqOffset == IndexedSequence.OFFSET_UNKNOWN) {
			final String pepSeq = peptideSequence.getSequence();
			seqOffset = proteinSequence.indexOf(pepSeq);
		}
		if (seqOffset == -1) {
			throw new RuntimeException("Could not get subsequence, unexpected error: peptide: " + peptideSequence
					+ ", protein: " + protein);
		}
		final int seqLen = peptideSequence.getSequenceLen();

		return Util.getResidues(peptideSequence, seqOffset, seqLen, proteinSequence);

	}

	@Override
	public List<IndexedSequence> getSequences(double precMass, double tolerance) throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		final List<IndexedSequence> ret = new ArrayList<IndexedSequence>();
		// logger.log(Level.FINE, "Starting peptide sequences query");
		// long start = System.currentTimeMillis();

		double minMass = precMass - tolerance;

		if (minMass < 0) {
			minMass = 0;
		}
		final double maxMass = precMass + tolerance;

		// System.out.println("==============" + precMass + " ==" + tolerance +
		// " " + minMass + " " + maxMass);
		final int[] bucketRange = getBucketsForMassRange(minMass, maxMass);
		if (bucketRange[0] > Constants.NUM_BUCKETS - 1 || bucketRange[1] > Constants.NUM_BUCKETS - 1) {
			logger.warn(
					"Cannot query, unsupported precursor mass: " + precMass + ". Max supported mass is: " + MAX_MASS);
			return ret;
		}

		for (int bucket = bucketRange[0]; bucket <= bucketRange[1]; ++bucket) {
			final List<IndexedSequence> sequencesPerBucket = buckets[bucket].getSequences(precMass, tolerance);
			ret.addAll(sequencesPerBucket);
		}

		// long end = System.currentTimeMillis();
		// logger.info("Peptide sequences query took: " + (end -
		// start) + " milisecs, results: " + ret.size());

		return ret;
	}

	@Override
	public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {
		if (ranges.size() == 1) {
			// simple case - single range, use the single range version
			final MassRange range = ranges.get(0);
			return getSequences(range.getPrecMass(), range.getTolerance());
		}

		// long start = System.currentTimeMillis();

		// handle multiple ranges
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		// convert mass ranges to intverals and merge overlapping mass ranges
		final ArrayList<Interval> intervals = new ArrayList<Interval>();
		for (final MassRange range : ranges) {
			final Interval ith = Interval.massRangeToInterval(range);
			intervals.add(ith);
		}
		final List<Interval> mergedIntervals = MergeIntervals.mergeIntervals(intervals);

		final List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

		// contruct composite queries per bucket
		// store mass ranges falling in every bucket
		// using non-generic list, because arrays don't support generic objects
		// and using [][] requires extra conversion later to adapt to API
		final List bucketIntervals[] = new ArrayList[buckets.length];

		for (final Interval massInterval : mergedIntervals) {

			final double minMass = massInterval.getStart();
			final double maxMass = massInterval.getEnd();

			// figure out which buckets the range falls in
			final int[] bucketRange = getBucketsForMassRange(minMass, maxMass);
			if (bucketRange[0] > Constants.NUM_BUCKETS - 1 || bucketRange[1] > Constants.NUM_BUCKETS - 1) {
				logger.error("Cannot query, unsupported precursor mass: " + minMass + "-" + maxMass
						+ ". Max supported mass is: " + MAX_MASS);
				return ret;
			}

			// add massrange to every bucket that qualifies
			for (int bucket = bucketRange[0]; bucket <= bucketRange[1]; ++bucket) {
				List<Interval> bucketR = bucketIntervals[bucket];
				if (bucketR == null) {
					// lazy-initialize the bucket first time
					bucketR = new ArrayList<Interval>();
					bucketIntervals[bucket] = bucketR;
				}

				// append mass range
				if (!bucketR.contains(massInterval)) {
					bucketR.add(massInterval);
				}
			}
		}

		// query the buckets that have mass ranges assigned
		for (int bucket = 0; bucket < buckets.length; ++bucket) {
			final List<Interval> queryBucketIntervals = bucketIntervals[bucket];
			if (queryBucketIntervals == null) {
				continue; // skip, no intervals for this bucket
			}

			final List<IndexedSequence> sequencesPerBucket = buckets[bucket]
					.getSequencesIntervals(queryBucketIntervals);
			ret.addAll(sequencesPerBucket);
		}

		// long end = System.currentTimeMillis();
		// logger.info("Peptide sequences query took: " + (end -
		// start) + " milisecs, results: " + ret.size());

		return ret;

	}

	@Override
	public boolean supportsProteinCache() {
		return true;
	}

	@Override
	public void setProteinCache(ProteinCache proteinCache) {
		this.proteinCache = proteinCache;
		for (int i = 0; i < Constants.NUM_BUCKETS; ++i) {
			buckets[i].setProteinCache(proteinCache);
		}
	}

	@Override
	public long addProteinDef(long num, String definition, String proteinSequence) throws DBIndexStoreException {
		// using protein cache entirely
		// nothing here
		return num;
	}

	@Override
	public List<IndexedProtein> getProteins(IndexedSequence sequence) throws DBIndexStoreException {
		final List<IndexedProtein> ret = new ArrayList<IndexedProtein>();

		for (final Integer protId : sequence.getProteinIds()) {
			final IndexedProtein ip = new IndexedProtein(proteinCache.getProteinDef(protId), protId);
			ret.add(ip);
		}

		return ret;

	}

	/**
	 * test driver
	 *
	 * @param args
	 */
	public static void main(String[] args) {
		DBIndexStore store = new DBIndexStoreSQLiteMult(SearchParams.getInstance(), false);

		final String protDef1 = "4R79.2 CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";
		final String protDef2 = "Reverse_4R79.2  CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";

		try {
			logger.info("Testing init");
			final File idx = new File("./test.fasta.idx");
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
			logger.info("Testing adds 1");
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
			logger.error(null, ex);
		}

		try {
			// test more adds with the same statement in a new transaction
			logger.info("Testing adds 2");
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
			logger.error(null, ex);
		}

		try {
			// test gets
			logger.info("Testing gets 1");
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
			logger.error(null, ex);
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
			logger.error(null, ex);
		}

		if (true) {
			return;
		}

		try {
			store = new DBIndexStoreSQLiteMult(SearchParams.getInstance(), false);
			store.init("EBI-IPI_Human_IPI_Human_3_85_06-29-2011_reversed.fasta");
			for (int i = 0; i < 1000; ++i) {
				logger.info("Testing getting sequence " + i);
				final List<IndexedSequence> res = store.getSequences(i, 3f);
				if (res.size() > 0) {
					final IndexedSequence seq = res.get(0);
					logger.info("Got " + res.size() + " sequences, first sequence: " + seq);
					final List<IndexedProtein> proteins = store.getProteins(seq);
					logger.info("Got " + res.size() + " sequences, first sequence's proteins: ");
					for (final IndexedProtein protein : proteins) {
						System.out.println(protein);
					}
				}

			}

		} catch (final DBIndexStoreException ex) {
			logger.error(null, ex);
		}

		try {
			logger.error("Testing getting multiple sequence proteins ");
			store = new DBIndexStoreSQLiteMult(SearchParams.getInstance(), false);
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
			logger.error(null, ex);
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
}
