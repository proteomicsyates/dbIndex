package edu.scripps.yates.dbindex;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.ws.rs.NotSupportedException;

import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.core.SimpleAnalyzer;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.DoubleField;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.IntField;
import org.apache.lucene.index.AtomicReader;
import org.apache.lucene.index.AtomicReaderContext;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.index.IndexWriterConfig.OpenMode;
import org.apache.lucene.index.IndexableField;
import org.apache.lucene.search.Collector;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.NumericRangeQuery;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.Scorer;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.util.Version;

import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.fasta.dbindex.MassRange;
import edu.scripps.yates.utilities.fasta.dbindex.ResidueInfo;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.TIntHashSet;

/**
 *
 * Index store implementation based on Lucene. To support large indexes and
 * search of large databases.
 *
 * @author Adam
 */
public class DBIndexStoreLucene implements DBIndexStore {

	private ProteinCache protCache;
	private final Document document = new Document();

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

	private enum FIELDS {

		MASS, PROT_ID, SEQUENCE, LR, RR, LEN, OFFSET
	};

	private final Field massField = new DoubleField(FIELDS.MASS.toString(), 0.0, Field.Store.YES);
	private final Field protIdField = new IntField(FIELDS.PROT_ID.toString(), 0, Field.Store.YES);
	private final Field lenField = new IntField(FIELDS.LEN.toString(), 0, Field.Store.YES);
	private final Field offsetField = new IntField(FIELDS.OFFSET.toString(), 0, Field.Store.YES);
	// private final Field seqField = new
	// StringField(FIELDS.SEQUENCE.toString(), "", Field.Store.NO);
	// private final Field lrField = new StringField(FIELDS.LR.toString(), "",
	// Field.Store.NO);
	// private final Field rrField = new StringField(FIELDS.RR.toString(), "",
	// Field.Store.NO);
	private final String INDEX_DIR_SUFFIX = ".idx2";
	private String indexDir;
	// indexing
	private IndexWriter indexWriter;
	// searching
	private IndexReader indexReader;
	private IndexSearcher indexSearcher;
	private Analyzer writeAnalyzer;
	// private QueryParser queryParser;
	private static final Logger logger = Logger.getLogger(DBIndexStoreLucene.class.getName());
	private long indexed = 0;
	private static final Version LUCENE_VERSION = Version.LUCENE_41;
	private static final String searchField = FIELDS.MASS.toString();

	@Override
	public void init(String databaseID) throws DBIndexStoreException {
		indexDir = databaseID + INDEX_DIR_SUFFIX;

		initWriter();

		initReader();

	}

	private void initWriter() throws DBIndexStoreException {
		// init fields
		document.add(massField);
		document.add(protIdField);
		document.add(lenField);
		document.add(offsetField);

		// analyzer = new StandardAnalyzer(LUCENE_VERSION);
		writeAnalyzer = new SimpleAnalyzer(LUCENE_VERSION);
		final IndexWriterConfig config = new IndexWriterConfig(LUCENE_VERSION, writeAnalyzer);

		try {
			final Directory dir = FSDirectory.open(new File(indexDir));
			config.setOpenMode(OpenMode.CREATE);
			config.setRAMBufferSizeMB(2048.0);
			// TODO disable auto commit
			indexWriter = new IndexWriter(dir, config);

			// document.add(new Field(FIELD_PATH, path, Field.Store.YES,
			// Field.Index.UN_TOKENIZED));

		} catch (final IOException e) {
			logger.log(Level.SEVERE, "Error setting up index", e);
			throw new DBIndexStoreException("Error setting up index", e);
		}
	}

	private void initReader() throws DBIndexStoreException {

		try {
			indexReader = DirectoryReader.open(FSDirectory.open(new File(indexDir)));
			indexSearcher = new IndexSearcher(indexReader);
		} catch (final org.apache.lucene.index.IndexNotFoundException e) {
			logger.log(Level.WARNING, "Index not found, OK if running in index mode " + indexDir, e);
			// TODO maybe in index mode, TOOD fix, may need to know the mode how
			// to initialize
		} catch (final IOException ex) {
			logger.log(Level.SEVERE, "Cannot initialize reader for indexDir: " + indexDir, ex);
			throw new DBIndexStoreException("Cannot initialize reader for indexDir: " + indexDir, ex);
		}

	}

	@Override
	public void startAddSeq() throws DBIndexStoreException {
	}

	@Override
	public void stopAddSeq() throws DBIndexStoreException {
		try {
			// indexWriter.addDocument(document);

			if (indexWriter != null) {
				logger.log(Level.INFO, "Commiting index");
				indexWriter.commit();
				logger.log(Level.INFO, "Closing writing index");
				indexWriter.close();
			}
		} catch (final IOException ex) {
			logger.log(Level.SEVERE, "Error closing the indexer", ex);
		}
	}

	@Override
	public boolean indexExists() throws DBIndexStoreException {
		final File dir = new File(indexDir);
		if (dir.exists() && dir.isDirectory() && dir.listFiles().length == 0) {
			return false;
		}

		return getNumberSequences() > 0;
	}

	@Override
	public FilterResult filterSequence(double precMass, String sequence) {
		// no filtering, we add every sequence
		return FilterResult.INCLUDE;
	}

	@Override
	public void addSequence(double precMass, int sequenceOffset, int sequenceLen, String sequence, String resLeft,
			String resRight, long proteinId) throws DBIndexStoreException {

		if (indexed % 1000000 == 0) {
			System.out.println(indexed);
		}
		indexed++;

		massField.setDoubleValue(precMass);
		protIdField.setIntValue((int) proteinId);
		lenField.setIntValue(sequenceLen);
		offsetField.setIntValue(sequenceOffset);

		try {
			indexWriter.addDocument(document);
			// indexWriter.commit();
		} catch (final IOException ex) {
			logger.log(Level.SEVERE, "Error indexing a sequence", ex);
			throw new DBIndexStoreException("Error indexing a sequence", ex);
		}

	}

	private List<IndexedSequence> runQuery(Query query) throws DBIndexStoreException {
		// System.out.println("Searching for: " + query.toString());
		// Date start = new Date();
		final ResultsCollector resultsCollector = new ResultsCollector();
		try {
			indexSearcher.search(query, resultsCollector);

		} catch (final IOException ex) {
			logger.log(Level.SEVERE, "Error searching with query: " + query.toString(), ex);
			throw new DBIndexStoreException("Error searching with query: " + query.toString(), ex);
		}

		// Date end = new Date();
		// System.out.println("Time: " + (end.getTime() - start.getTime()) +
		// "ms");

		return resultsCollector.getSequences();
	}

	@Override
	public List<IndexedSequence> getSequences(double precMass, double tolerance) throws DBIndexStoreException {

		if (indexSearcher == null) {
			throw new DBIndexStoreException(
					"Cannot search, searcher not initialized on the index, are you running in indexing mode?");
		}
		// QueryParser queryParser = new QueryParser(LUCENE_VERSION,
		// searchField, readAnalyzer);

		double massLow = precMass - tolerance;
		if (massLow < 0.0) {
			massLow = 0.0;
		}
		final double massHigh = precMass + tolerance;

		final Query query = NumericRangeQuery.newDoubleRange(searchField, massLow, massHigh, true, true);

		return runQuery(query);
	}

	@Override
	public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {

		if (indexSearcher == null) {
			throw new DBIndexStoreException(
					"Cannot search, searcher not initialized on the index, are you running in indexing mode?");
		}

		final org.apache.lucene.search.BooleanQuery composite = new org.apache.lucene.search.BooleanQuery();

		for (final MassRange range : ranges) {

			final double precMass = range.getPrecMass();
			final double tolerance = range.getTolerance();
			double massLow = precMass - tolerance;
			if (massLow < 0.0) {
				massLow = 0.0;
			}
			final double massHigh = precMass + tolerance;

			final Query query = NumericRangeQuery.newDoubleRange(searchField, massLow, massHigh, true, true);
			composite.add(query, org.apache.lucene.search.BooleanClause.Occur.SHOULD);
		}

		return runQuery(composite);
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
		if (indexReader == null) {
			return 0;
		}
		return indexReader.numDocs();
	}

	@Override
	public ResidueInfo getResidues(IndexedSequence peptideSequence, IndexedProtein protein)
			throws DBIndexStoreException {
		// in this implementation residues are already stored in the index
		return new ResidueInfo(peptideSequence.getResLeft(), peptideSequence.getResLeft());
	}

	/**
	 * test driver
	 *
	 * @param args
	 */
	public static void main(String[] args) {
		DBIndexStore store = new DBIndexStoreLucene();

		final String protDef1 = "4R79.2 CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";
		final String protDef2 = "Reverse_4R79.2  CE19650 WBGene00007067 Ras family status:Partially_confirmed TR:Q9XXA4 protein_id:CAA20282.1";

		try {
			logger.log(Level.INFO, "Testing init");
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
			store = new DBIndexStoreLucene();
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
			store = new DBIndexStoreLucene();
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

	private class ResultsCollector extends Collector {
		// to collapse multiple sequences into single one, with multproteins

		private final Map<String, List<IndexedSeqInternal>> temp = new THashMap<String, List<IndexedSeqInternal>>();
		private AtomicReaderContext curReader;

		@Override
		public void setScorer(Scorer scorer) throws IOException {
		}

		@Override
		public void collect(int i) throws IOException {
			final AtomicReader reader = curReader.reader();

			final Document doc = reader.document(i);

			// System.out.println("DOC: " + doc.toString());

			final List<IndexableField> fields = doc.getFields();
			if (fields.isEmpty()) {
				return;
			}

			final IndexableField mass = doc.getField(FIELDS.MASS.toString());
			final double massF = mass.numericValue().doubleValue();
			final IndexableField offset = doc.getField(FIELDS.OFFSET.toString());
			final int offsetI = offset.numericValue().intValue();
			final IndexableField protId = doc.getField(FIELDS.PROT_ID.toString());
			final int protIdI = protId.numericValue().intValue();
			final IndexableField len = doc.getField(FIELDS.LEN.toString());
			final int lenI = len.numericValue().intValue();

			String peptideSequence = null;
			// System.out.print("prot id: " + protId + " offset: " + offsetI +
			// " len: " + lenI);

			peptideSequence = protCache.getPeptideSequence(protIdI, offsetI, lenI);
			// System.out.println("pep: " + peptideSequence);

			// we are cheating and supplying protein id instead of peptide id
			// to set it temporarily, before we merge them into a list
			final IndexedSeqInternal tempSequence = new IndexedSeqInternal(massF, offsetI, lenI, protIdI,
					peptideSequence);
			List<IndexedSeqInternal> sequences = temp.get(peptideSequence);
			if (sequences == null) {
				sequences = new ArrayList<IndexedSeqInternal>();
				temp.put(peptideSequence, sequences);
			}
			sequences.add(tempSequence);

		}

		List<IndexedSequence> getSequences() {
			final List<IndexedSequence> ret = new ArrayList<IndexedSequence>();
			// group the same peptides from many proteins into single peptide
			// with protein id list
			for (final String pepSeqKey : temp.keySet()) {

				final List<IndexedSeqInternal> sequences = temp.get(pepSeqKey);

				final TIntHashSet proteinIds = new TIntHashSet();
				for (final IndexedSeqInternal tempSeq : sequences) {
					proteinIds.add(tempSeq.proteinId);
				}

				final IndexedSeqInternal firstSeq = sequences.get(0);
				final IndexedSequence mergedSequence = new IndexedSequence(0, firstSeq.mass, pepSeqKey, "", "");
				final ArrayList<Integer> proteinIds2 = new ArrayList<Integer>();
				for (final int integer : proteinIds._set) {
					proteinIds2.add(integer);
				}
				mergedSequence.setProteinIds(proteinIds2);
				// set residues
				final String protSequence = protCache.getProteinSequence(firstSeq.proteinId);
				final ResidueInfo residues = Util.getResidues(null, firstSeq.offset, firstSeq.length, protSequence);
				mergedSequence.setResidues(residues);
				ret.add(mergedSequence);
			}

			return ret;
		}

		@Override
		public void setNextReader(AtomicReaderContext arc) throws IOException {
			curReader = arc;
		}

		@Override
		public boolean acceptsDocsOutOfOrder() {
			return true;
		}
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
