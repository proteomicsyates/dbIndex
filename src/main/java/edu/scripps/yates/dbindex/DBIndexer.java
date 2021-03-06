package edu.scripps.yates.dbindex;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.FilenameUtils;

import com.compomics.util.general.UnknownElementMassException;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.proteoform.fasta.ProteoFormFastaReader;
import edu.scripps.yates.dbindex.DBIndexStore.FilterResult;
import edu.scripps.yates.dbindex.io.SearchParamReader;
import edu.scripps.yates.dbindex.model.Util;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.fasta.Fasta;
import edu.scripps.yates.utilities.fasta.FastaReader;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexType;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.fasta.dbindex.MassRange;
import edu.scripps.yates.utilities.fasta.dbindex.PeptideFilter;
import edu.scripps.yates.utilities.masses.AssignMass;
import edu.scripps.yates.utilities.masses.FormulaCalculator;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import gnu.trove.set.hash.THashSet;

/**
 *
 * Supports storing and retrieving proteins and sequences to/from index
 *
 * There can be different underlying store index implementations, implementing
 * DBIndexStore interface
 *
 * Supports modes: indexing only (store), search from indexed store, and search
 * without a store (cut on the fly) with temporary in-memory store
 *
 */
public class DBIndexer {

	/**
	 * Different operational modes supported
	 */
	public enum IndexerMode {

		INDEX, SEARCH_INDEXED, SEARCH_UNINDEXED,
	};

	private final IndexerMode mode;
	private IndexType indexType;
	protected DBIndexStore indexStore;
	protected final DBIndexSearchParams sparam;
	// private double massTolerance;
	private boolean inited;
	private String indexName; // index name that is params-specific
	protected long protNum; // protein number in sequence, starting at 1
	protected ProteinCache proteinCache;
	protected final PeptideFilter peptideFilter;
	private Pattern decoy;
	private static final org.apache.log4j.Logger logger = org.apache.log4j.Logger.getLogger(DBIndexer.class);
	protected static final int MAX_SEQ_LENGTH = 10000;
	private static final FormulaCalculator formulaCalculator = new FormulaCalculator();

	public static void main(String[] args) {

		if (args.length == 0) {
			System.out.println("Need a directory name");
			System.exit(1);
		}

		runDBIndexer(args[0]);
	}

	public static void runDBIndexer(String path) {
		runDBIndexer(path, SearchParamReader.DEFAULT_PARAM_FILE_NAME);
	}

	public static void mergeLargeDBIndex(String path) {
		SearchParamReader preader;
		try {
			// System.out.println("path" + path);
			preader = new SearchParamReader(path, SearchParamReader.DEFAULT_PARAM_FILE_NAME);
			// preader = new SearchParamReader(args[0],
			// SearchParamReader.DEFAULT_PARAM_FILE_NAME);
		} catch (final IOException ex) {
			logger.error("Error getting params", ex);
			return;
		}
		final DBIndexSearchParams sparam = preader.getParam();

		if (!sparam.getIndexType().equals(IndexType.INDEX_LARGE)) {
			throw new RuntimeException("Merge only works for large db index type, other indexed have merge built-in");
		}

		// indexer in indexing mode
		final DBIndexStoreFiles largeIndex = new DBIndexStoreFiles();
		try {
			largeIndex.merge(sparam);

		} catch (final Exception ex) {
			logger.error("Error running merge on large index.", ex);
		}

	}

	public static void runDBIndexer(String path, String paramFile) {
		SearchParamReader preader;
		try {
			preader = new SearchParamReader(path, paramFile);
			// preader = new SearchParamReader(args[0],
			// SearchParamReader.DEFAULT_PARAM_FILE_NAME);
		} catch (final IOException ex) {
			logger.error("Error getting params", ex);
			return;
		}
		final DBIndexSearchParams sparam = preader.getParam();

		// indexer in indexing mode
		final DBIndexer indexer = new DBIndexer(sparam, IndexerMode.INDEX);
		try {
			indexer.init();
			indexer.run();
		} catch (final DBIndexerException ex) {
			logger.error("Error running indexer in the index mode.", ex);
		}

	}

	public DBIndexer(edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams sparam, IndexerMode mode,
			DBIndexStore dbIndexStore) {
		IndexUtil.setNumRowsToLookup(0);
		this.sparam = sparam;
		this.mode = mode;
		peptideFilter = sparam.getPeptideFilter();
		indexStore = dbIndexStore;
		inited = false;
		protNum = -1;
		if (sparam.getDiscardDecoyRegexp() != null) {
			decoy = Pattern.compile(sparam.getDiscardDecoyRegexp());
		}
	}

	/**
	 * Create new db indexer, passing params and indexing mode, and the
	 *
	 * @param sparam
	 * @param indexMode       indexing mode (search indexed, search unindexed, index
	 *                        only)
	 * @param massGroupFactor this is the factor by which each double mass will be
	 *                        multiplied to get the key in the index. Using any
	 *                        other constructor, the massGroupFactor is 10000
	 */
	public DBIndexer(edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams sparam, IndexerMode mode) {

		this.sparam = sparam;
		this.mode = mode;
		peptideFilter = sparam.getPeptideFilter();

		// this.massTolerance = sparam.getRelativePeptideMassTolerance();

		if (!mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
			indexType = sparam.getIndexType();
			if (indexType.equals(IndexType.INDEX_NORMAL) && !sparam.isUsingMongoDB()) {
				boolean inMemoryIndex = sparam.isInMemoryIndex();
				// if in memory index, use in-memory for search only
				// in future we can try indexing in-memory, but it will fail for
				// larger dbs
				if (mode.equals(IndexerMode.INDEX)) {
					inMemoryIndex = false;
				}
				indexStore = new DBIndexStoreSQLiteMult(sparam, inMemoryIndex);
				// /this.indexStore = new DBIndexStoreSQLiteString();
			} else {
				// index_type=2
				// this.indexStore = new DBIndexStoreFiles();
				indexStore = new DBIndexStoreMongoDb(sparam);
				// //this.indexStore = new DBIndexStoreLucene();
			}
			logger.info("Using database index type: " + indexType);

		} else {
			// set up in memory temporary "index" that does the filtering
			// this is the nonindexed saerch
			indexStore = new MassRangeFilteringIndex();
			logger.info("Using no-index for the search. ");
		}

		inited = false;
		protNum = -1;
		if (sparam.getDiscardDecoyRegexp() != null) {
			decoy = Pattern.compile(sparam.getDiscardDecoyRegexp());
		}
	}

	/**
	 * Cut a Fasta protein sequence according to params spec and index the protein
	 * and generated sequences
	 *
	 * Should be called once per unique Fasta sequence
	 *
	 * @param fasta fasta to index
	 * @throws IOException
	 */
	private void cutSeq(final Fasta fasta) throws IOException {
		// change by SALVA in order to keep all fasta information in the index
		// final String protAccession = fasta.getSequestLikeAccession();
		final String protFastaHeader = fasta.getOriginalDefline();
		final String protSeq = fasta.getSequence();

		cutSeq(protFastaHeader, protSeq);

	}

	/**
	 * Cut a Fasta protein sequence according to params spec and index the protein
	 * and generated sequences
	 *
	 * Should be called once per unique Fasta sequence
	 *
	 * @param fasta fasta to index
	 * @throws IOException
	 */
	protected void cutSeq(final String protAccession, String protSeq) throws IOException {
		// Enzyme enz = sparam.getEnzyme();
		int length = protSeq.length();
		//
		// AssignMass aMass = AssignMass.getInstance(true);

		final char[] pepSeq = new char[MAX_SEQ_LENGTH]; // max seq length
		int curSeqI = 0;
		double aaMass = 0;
		final int maxIntCleavage = sparam.getMaxMissedCleavages();
		// at least one of this AAs has to be in the sequence:
		final char[] mandatoryInternalAAs = sparam.getMandatoryInternalAAs();

		try {
			final long proteinId = indexStore.addProteinDef(++protNum, protAccession, protSeq);
			// System.out.println(fasta.getSequestLikeAccession());
			// System.out.println(fasta.getDefline());
			// System.out.println(protSeq);

			for (int start = 0; start < length; ++start) {
				int end = start;

				// clear the preallocated seq byte array
				// Arrays.fill(seq, 0, curSeqI > 0?curSeqI-1:0, (byte) 0); //no
				// need, we copy up to curSeqI nowu
				curSeqI = 0;

				// double precMass = Constants.H2O_PROTON_SCALED_DOWN;
				double precMass = 0;

				// Salva added 24Nov2014
				if (sparam.isH2OPlusProtonAdded())
					precMass += AssignMass.H2O_PROTON;
				precMass += AssignMass.getcTerm();
				precMass += AssignMass.getnTerm();
				// System.out.println("===>>" + precMass + "\t" +
				// Constants.MAX_PRECURSOR);
				// System.out.println("==" + j + " " + length + " " + (j <
				// length));

				// int testC=0;
				int pepSize = 0;

				int intMisCleavageCount = -1;

				// while (precMass <= Constants.MAX_PRECURSOR_MASS && end <
				// length) {
				while (precMass <= sparam.getMaxPrecursorMass() && end < length) {
					pepSize++;

					final char curIon = protSeq.charAt(end);
					if (curIon == '[') {
						logger.debug("starting formula at " + end);
						final StringBuilder chemicalFormula = new StringBuilder();
						while (protSeq.charAt(++end) != ']') {
							chemicalFormula.append(protSeq.charAt(end));
						}

						aaMass = formulaCalculator.calculateMass(chemicalFormula.toString());

						// remove the formula from the protein sequence
						protSeq = protSeq.replace("[" + chemicalFormula + "]", "");
						// position index properly
						end -= chemicalFormula.length() + 2;
						length = protSeq.length();
						logger.info("Found PTM formula " + chemicalFormula.toString() + " equals to " + aaMass
								+ " at position " + (end + 1) + " of protein " + protAccession);
					} else {
						pepSeq[curSeqI++] = curIon;
						aaMass = AssignMass.getMass(curIon);
					}
					precMass = precMass + aaMass;
					final String peptideSeqString = String.valueOf(Arrays.copyOf(pepSeq, curSeqI));
					if (peptideFilter != null && !peptideFilter.isValid(peptideSeqString)) {
						logger.debug("Skipping peptide '" + peptideSeqString + "' by filter");
						break;
					}
					if (sparam.getEnzyme().isEnzyme(protSeq.charAt(end))) {
						intMisCleavageCount++;
					}

					final boolean cleavageStatus = sparam.getEnzyme().checkCleavage(protSeq, start, end,
							sparam.getEnzymeNocutResidues());
					if (cleavageStatus) {

						if (intMisCleavageCount > maxIntCleavage) {
							break;
						}

						// if (precMass > Constants.MAX_PRECURSOR_MASS) {
						if (precMass > sparam.getMaxPrecursorMass()) {
							break;
						}

						if (pepSize >= Constants.MIN_PEP_LENGTH && precMass >= sparam.getMinPrecursorMass()) { // Constants.MIN_PRECURSOR
																												// )

							if (mandatoryInternalAAs != null) {
								boolean found = false;
								for (final char internalAA : mandatoryInternalAAs) {
									if (peptideSeqString.indexOf(internalAA) >= 0) {
										found = true;
									}
								}
								if (!found) {
									break;
								}
							}

							// qualifies based on params

							// check if index will accept it
							final FilterResult filterResult = indexStore.filterSequence(precMass, peptideSeqString);

							if (filterResult.equals(FilterResult.SKIP_PROTEIN_START)) {
								// bail out earlier as we are no longer
								// interested in this protein starting at start
								break; // move to new start position
							} else if (filterResult.equals(FilterResult.INCLUDE)) {

								final int resLeftI = start >= Constants.MAX_INDEX_RESIDUE_LEN
										? start - Constants.MAX_INDEX_RESIDUE_LEN
										: 0;
								final int resLeftLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, start);
								final StringBuilder sbLeft = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
								for (int ii = 0; ii < resLeftLen; ++ii) {
									sbLeft.append(protSeq.charAt(ii + resLeftI));
								}
								final int resRightI = end + 1;
								final int resRightLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, length - end - 1);
								final StringBuilder sbRight = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
								if (resRightI < length) {
									for (int jj = 0; jj < resRightLen; ++jj) {
										sbRight.append(protSeq.charAt(jj + resRightI));
									}
								}

								// add -- markers to fill
								// Constants.MAX_INDEX_RESIDUE_LEN length
								final int lLen = sbLeft.length();
								for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - lLen; ++c) {
									sbLeft.insert(0, '-');
								}
								final int rLen = sbRight.length();
								for (int c = 0; c < Constants.MAX_INDEX_RESIDUE_LEN - rLen; ++c) {
									sbRight.append('-');
								}

								final String resLeft = sbLeft.toString();
								final String resRight = sbRight.toString();

								indexStore.addSequence(precMass, start, curSeqI, peptideSeqString, resLeft, resRight,
										proteinId);
								// System.out.println("\t" + peptideSeqString);
							} // end if add sequence
						}
					}
					++end;
				}
				//
			}
		} catch (final DBIndexStoreException e) {
			logger.error("Error writing sequence to db index store, ", e);
		} catch (final UnknownElementMassException e) {
			e.printStackTrace();
			logger.error("Error writing sequence to db index store, ", e);
		}

	}

	/**
	 * Initialize the indexer
	 *
	 * @throws IOException
	 */
	public void init() throws DBIndexerException {
		if (inited) {
			throw new IllegalStateException("Already inited");
		}

		// reset prot number
		protNum = -1;

		final String dbName = sparam.getDatabaseName();

		if (!sparam.isUsingMongoDB()) {
			final File dbFile = new File(dbName);
			if (!dbFile.exists() || !dbFile.canRead()) {
				throw new DBIndexerException("Cannot read (and index) the Fasta database file: " + dbName);
			}
		}
		// get param specific name for the inde0x
		indexName = createFullIndexFileName();

		try {
			// initialize index storage
			indexStore.init(indexName);
			if (mode.equals(IndexerMode.SEARCH_INDEXED)) {
				// if (mode.equals(IndexerMode.SEARCH_INDEXED) &&
				// sparam.getIndexType() == DBIndexer.IndexType.INDEX_NORMAL) {

				if (!indexStore.indexExists()) {
					logger.error("Error, cannot search, index does not exist: " + indexName);
					throw new DBIndexerException("Error, cannot search, index does not exist: " + indexName);
				}
				setProteinCache();
			} else if (mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
				setProteinCache();
			}
			inited = true;
		} catch (final DBIndexStoreException ex) {
			logger.error("Could not initialize index storage.", ex);
			throw new DBIndexerException("Could not initialize index storage.");
		}
	}

	protected String createFullIndexFileName() {
		return IndexUtil.createFullIndexFileName(sparam);
	}

	/**
	 * Set protein cache if storage supports it Useful if database has already been
	 * indexed Otherwise, the cache can be populated while indexing is done
	 */
	private void setProteinCache() {
		if (!indexStore.supportsProteinCache()) {
			return;
		}

		final ProteinCache protCache = getProteinCache();

		// ensure only 1 thread at time enters this
		synchronized (protCache) {

			logger.info("Populating protein cache");
			if (protCache.isPopulated() == false) {
				logger.info("Initializing protein cache");

				Iterator<Fasta> itr = null;
				try {
					itr = IndexUtil.getFastaReader(sparam).getFastas();
				} catch (final IOException ex) {
					logger.error("Could not set protein cache", ex);
					return;
				}

				while (itr.hasNext()) {
					final Fasta fasta = itr.next();
					// SALVA CHANGE
					// protCache.addProtein(fasta.getSequestLikeAccession(),
					// fasta.getSequence());
					protCache.addProtein(fasta.getDefline(), fasta.getSequence());
				}

				logger.info("Done initializing protein cache");
			} else {
				// logger.info(
				// "Protein cache already populated, reusing.");
			}
		}

		indexStore.setProteinCache(protCache);
		// logger.info("Done populating protein cache");

	}

	/**
	 * Indexes the database if database index does not exists
	 *
	 * @throws IOException
	 */
	public void run() throws DBIndexerException {

		if (!inited) {
			throw new IllegalStateException("Not initialized.");
		}

		// if search unindexed, do nothing, will use protein cache to cutSeq on
		// the fly
		if (mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
			return;
		}

		try {

			if (mode.equals(IndexerMode.INDEX)) {
				if (indexStore.indexExists()) {
					logger.info("Found existing index, skipping indexing.");

					// skip indexing
					return;
				} else {
					logger.warn("Found no sequences in the index, will now index.");
				}
			}
		} catch (final DBIndexStoreException ex) {
			logger.error("Could not query number of sequences", ex);
			return;
		}

		FileInputStream fis = null; // fasta reader

		// setup status writer
		final String statusFilePath = Util.getFileBaseName(sparam.getDatabaseName()) + "log";
		FileWriter statusWriter = null;
		int totalProteins = 0;
		String totalProteinsString = null;
		int indexedProteins = 0;
		final long t0 = System.currentTimeMillis();
		final FastaReader fastaReader = IndexUtil.getFastaReader(sparam);
		try {
			statusWriter = new FileWriter(statusFilePath);
		} catch (final IOException ex) {
			logger.error("Error initializing index progress writer for file path: " + statusFilePath, ex);
		}

		// start indexing
		try {

			// getUniprotAnnorations
			logger.info("Reading FASTA file to retrieve Uniprot annotations from Internet...");
			final long t1 = System.currentTimeMillis();

			final Set<String> accs = fastaReader.getACCsFromFasta();

			if (accs == null || accs.isEmpty()) {
				throw new DBIndexerException("Reading FASTA file '" + sparam.getDatabaseName()
						+ "' was not able to extract any single Uniprot protein accession.\nUse a uniprot formated fasta database.");
			}
			logger.info(accs.size()
					+ " proteins extracted from fasta file with Uniprot accession. Now looking in the local index.");

			final Iterator<Fasta> itr = fastaReader.getFastas();
			try {
				totalProteins = fastaReader.getNumberFastas();
				totalProteinsString = String.valueOf(totalProteins);
			} catch (final IOException ex) {
				totalProteins = accs.size();
			}
			final UniprotProteinLocalRetriever uplr = getUniprotProteinLocalRetriever();
			if (uplr != null) {
				final boolean retrieveFastaIsoforms = true;

				uplr.getAnnotatedProteins(sparam.getUniprotVersion(), accs, retrieveFastaIsoforms,
						isRetrieveFastaIsoformsFromMainForms());
				final long t2 = System.currentTimeMillis() - t1;
				logger.info("Uniprot annotations were saved to local file system at '"
						+ sparam.getUniprotReleasesFolder().getAbsolutePath() + "' in "
						+ DatesUtil.getDescriptiveTimeFromMillisecs(t2));
			}
			// index, set protein cache on the fly
			indexStore.startAddSeq();

			fis = new FileInputStream(sparam.getDatabaseName());

			// create prot cache, in case searcher is kicked off in the same
			// process after indexing is done
			final ProteinCache protCache = getProteinCache();
			indexStore.setProteinCache(protCache);
			final Set<String> uniprotAccs = new THashSet<String>();
			final ProgressCounter counter = new ProgressCounter(totalProteins, ProgressPrintingType.PERCENTAGE_STEPS, 0,
					true);
			int decoyDiscarded = 0;
			while (itr.hasNext()) {

				final Fasta fasta = itr.next();
				// change by SALVA in order to keep all the information of the
				// header in the index
				protCache.addProtein(fasta.getOriginalDefline(), fasta.getSequence());
				// final String uniprotACC =
				// FastaParser.getACC(fasta.getOriginalDefline()).getFirstelement();
				final String uniprotACC = fasta.getAccession();
				if (decoy != null) {
					final Matcher matcher = decoy.matcher(uniprotACC);
					if (matcher.find()) {
						decoyDiscarded++;
						continue;
					}
				}
				cutSeq(fasta);

				uniprotAccs.add(uniprotACC);
				if (fastaReader instanceof ProteoFormFastaReader) {
					counter.setProgress(uniprotAccs.size());
				} else {
					counter.increment();
				}

				++indexedProteins;
				if (statusWriter != null) {
					statusWriter.append(totalProteinsString).append("\t").append(Integer.toString(indexedProteins))
							.append("\n");
					if (indexedProteins % 100 == 0) {
						statusWriter.flush();
					}
				}

				final String printIfNecessary = counter.printIfNecessary();
				if (!"".equals(printIfNecessary)) {
					logger.info("Proteins indexed: " + printIfNecessary);

				}

			}
			// System.out.print("Printing the last buffer....");
			// indexStore.lastBuffertoDatabase();

		} catch (final IOException ex) {
			ex.printStackTrace();
			logger.error("Error initializing and adding sequences", ex);
			throw new DBIndexerException("Error initializing and adding sequences", ex);
		} catch (final DBIndexStoreException ex) {
			ex.printStackTrace();
			logger.error("Error initializing and adding sequences", ex);
			throw new DBIndexerException("Error initializing and adding sequences", ex);
		} catch (final Exception ex) {
			ex.printStackTrace();
			logger.error("Error initializing and adding sequences", ex);
			throw new DBIndexerException("Error initializing and adding sequences", ex);
		} finally {
			try {
				if (fis != null) {
					try {
						fis.close();
					} catch (final IOException ex) {
						ex.printStackTrace();
						logger.error("Cannot close stream", ex);
					}
				}
				indexStore.stopAddSeq();

				if (statusWriter != null) {
					try {
						statusWriter.flush();
						statusWriter.close();
					} catch (final IOException ex) {
						ex.printStackTrace();
						logger.warn("Errir closing index status writer", ex);
					}
				}

			} catch (final DBIndexStoreException ex) {
				ex.printStackTrace();
				logger.error("Error finalizing adding sequences", ex);
			}
		}

	}

	protected boolean isRetrieveFastaIsoformsFromMainForms() {
		return false;
	}

	protected UniprotProteinLocalRetriever getUniprotProteinLocalRetriever() {
		if (sparam.getUniprotReleasesFolder() != null) {
			final boolean useIndex = true;
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
					sparam.getUniprotReleasesFolder(), useIndex);
			return uplr;
		}
		return null;
	}

	protected ProteinCache getProteinCache() {
		if (proteinCache == null) {
			proteinCache = new ProteinCache();
		}
		return proteinCache;
	}

	private List<IndexedSequence> cutAndSearch(List<MassRange> massRanges) {
		if (!mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
			throw new RuntimeException("Cut and search only supported for SEARCH_UNINDEXED mode !");
		}

		// set up in memory temporary "index" that does the filtering
		((MassRangeFilteringIndex) indexStore).init(massRanges);

		// reset prot number
		protNum = -1;

		// cut sequences in prot cache and check if qualify

		final ProteinCache protCache = getProteinCache();

		final int numProteins = protCache.getNumberProteins();

		for (int prot = 0; prot < numProteins; ++prot) {
			final String protDef = protCache.getProteinDef(prot);
			final String protSequence = protCache.getProteinSequence(prot);

			try {
				// stores sequences that match only in our temporary "index"
				cutSeq(protDef, protSequence);
			} catch (final IOException ex) {
				logger.error("Error cutting sequence in no-index mode: " + protSequence, ex);
			}

		}

		List<IndexedSequence> sequences = null;
		try {
			sequences = indexStore.getSequences(massRanges);

		} catch (final DBIndexStoreException ex) {
			logger.error("Error getting sequences from in-memory filtering index", ex);
		}

		return sequences;

	}

	/**
	 * Get list of sequences in index for the mass, and using tolerance in
	 * SearchParams. Mass dependant tolerance is already scaled/precalculated
	 *
	 * Requires call to init() and run() first, which might perform indexing, if
	 * index for current set of sequest params does not exist
	 *
	 * @param precursorMass     mass to select by, in Da
	 * @param massToleranceInDa mass tolerance, already calculated for that mass as
	 *                          it is mass-dependant. In Da.
	 * @return list of matching sequences
	 * @throws DBIndexStoreException
	 */
	public List<IndexedSequence> getSequencesUsingDaltonTolerance(double precursorMass, double massToleranceInDa)
			throws DBIndexStoreException {
		if (mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
			final MassRange range = new MassRange(precursorMass, massToleranceInDa);
			final List<MassRange> ranges = new ArrayList<MassRange>();
			ranges.add(range);
			return cutAndSearch(ranges);
		} else {
			return indexStore.getSequences(precursorMass, massToleranceInDa);
		}
	}

	/**
	 * Get list of sequences in index for the mass, and using tolerance in
	 * SearchParams. Mass dependant tolerance is NOT already scaled/precalculated,
	 * since it is in PPM
	 *
	 * Requires call to init() and run() first, which might perform indexing, if
	 * index for current set of sequest params does not exist
	 *
	 * @param precursorMass      mass to select by, in Da
	 * @param massToleranceInPPM mass tolerance in PPM.
	 * @return list of matching sequences
	 * @throws DBIndexStoreException
	 */
	public List<IndexedSequence> getSequencesUsingPPMTolerance(double precursorMass, double massToleranceInPPM)
			throws DBIndexStoreException {

		final double massTolerance = IndexUtil.getToleranceInDalton(precursorMass, massToleranceInPPM);

		if (mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
			final MassRange range = new MassRange(precursorMass, massTolerance);
			final List<MassRange> ranges = new ArrayList<MassRange>();
			ranges.add(range);
			return cutAndSearch(ranges);
		} else {

			final List<IndexedSequence> sequences = indexStore.getSequences(precursorMass, massTolerance);

			// using a predefined PRECISION, go up in the mass range to
			// see
			// if lower bound of higher masses PPM range in lower than
			// the
			// actual precursor mass, and in that case, search for that
			// mass.

			double upperBound = precursorMass + massTolerance;

			while (true) {
				// calculate the tolerance of the lowerBound
				final double massTolerance2 = IndexUtil.getToleranceInDalton(upperBound, massToleranceInPPM);
				final double lowerBoundOfUpperBoundMass = upperBound - massTolerance2;
				if (lowerBoundOfUpperBoundMass < precursorMass) {
					final List<IndexedSequence> sequences2 = indexStore.getSequences(upperBound, 0.0f);
					if (sequences2.isEmpty()) {
						break;
					} else {

						int numNewSeqs = 0;
						for (final IndexedSequence indexedSequence2 : sequences2) {
							if (!sequences.contains(indexedSequence2)) {
								sequences.add(indexedSequence2);
								numNewSeqs++;
							}
						}
						if (numNewSeqs > 0)
							logger.info(numNewSeqs + " new sequences looking in upper mass " + upperBound
									+ " from precursor mass " + precursorMass);
					}
				} else {
					break;
				}
				final double newupperBound = upperBound + Constants.PRECISION;
				if (newupperBound == upperBound)
					break;
				upperBound = newupperBound;

			}

			return sequences;

		}
	}

	/**
	 * Get list of sequences in index for the ranges specified wit mass . Requires
	 * call to init() and run() first, which might perform indexing, if index for
	 * current set of sequest params does not exist
	 *
	 * @param massRanges mass ranges to query
	 * @return list of matching sequences
	 * @throws DBIndexStoreException
	 */
	public List<IndexedSequence> getSequences(List<MassRange> massRanges) throws DBIndexStoreException {
		List<IndexedSequence> sequences = null;
		final long t1 = System.currentTimeMillis();
		if (mode.equals(IndexerMode.SEARCH_UNINDEXED)) {
			sequences = cutAndSearch(massRanges);
		} else {

			sequences = indexStore.getSequences(massRanges);

		}
		final long t2 = System.currentTimeMillis();
		// System.out.println("DEBUG DBINdexer.getSequences(), seqs: " +
		// sequences.size() + ", time: " + ((t2-t1)) + "ms");

		// logger.info("getSequences(): " + sequences.size());
		return sequences;
	}

	/**
	 * Get proteins associated with the sequence from the index. Requires call to
	 * init() and run() first, which might perform indexing, if index for current
	 * set of sequest params does not exist
	 *
	 * @param seq peptide sequence
	 * @return list of indexed protein objects associated with the sequence
	 * @throws DBIndexStoreException
	 */
	public List<IndexedProtein> getProteins(IndexedSequence seq) throws DBIndexStoreException {

		return indexStore.getProteins(seq);

	}

	/*
	 * private String getFullIndexFileName(String baseName) { String uniqueIndexName
	 * = baseName + "_"; //generate a unique string based on current params that
	 * affect the index final StringBuilder uniqueParams = new StringBuilder();
	 * //uniqueParams.append(sparam.getEnzyme().toString());
	 * //uniqueParams.append(sparam.getEnzymeNumber());
	 * uniqueParams.append(sparam.getEnzymeOffset());
	 * uniqueParams.append(sparam.getEnzymeResidues());
	 * uniqueParams.append(sparam.getEnzymeNocutResidues());
	 * uniqueParams.append("\nCleav: ");
	 * uniqueParams.append(sparam.getMaxInternalCleavageSites());
	 * uniqueParams.append(sparam.getMaxMissedCleavages());
	 * uniqueParams.append(sparam.getMaxNumDiffMod());
	 * uniqueParams.append("\nMods:"); for (final ModResidue mod :
	 * sparam.getModList()) { uniqueParams.append(mod.toString()).append(" "); }
	 * uniqueParams.append("\nMods groups:"); for (final List<double> modGroupList :
	 * sparam.getModGroupList()) { for (final double f : modGroupList) {
	 * uniqueParams.append(f).append(" "); } } final String uniqueParamsStr =
	 * uniqueParams.toString(); //logger.info( "Unique params: " + uniqueParamsStr);
	 * }
	 */
	public void close() {
	}

	public boolean isInited() {
		return inited;
	}

	/**
	 * Get proteins associated with the sequence. Requires call to init() and run()
	 * first, which might perform indexing, if index for current set of sequest
	 * params does not exist
	 *
	 * @param seq peptide sequence
	 * @return list of indexed protein objects associated with the sequence
	 * @throws DBIndexStoreException
	 */
	public Set<IndexedProtein> getProteins(String seq) throws DBIndexStoreException {
		final Set<IndexedProtein> ret = new THashSet<IndexedProtein>();

		if (!sparam.isUsingSeqDB()) {
			// get the mass of the sequence
			final boolean h2oPlusProtonAdded = sparam.isH2OPlusProtonAdded();
			final double mass = IndexUtil.calculateMass(seq, h2oPlusProtonAdded);
			// get the indexed sequences in the database
			final List<IndexedSequence> indexedSequences = getSequencesUsingDaltonTolerance(mass, 0.0f);
			for (final IndexedSequence indexedSequence : indexedSequences) {
				final String sequence = indexedSequence.getSequence();
				if (sequence.equals(seq)) {
					ret.addAll(indexStore.getProteins(indexedSequence));
				}
			}
		} else {
			final IndexedSequence indexPeptide = new IndexedSequence();
			indexPeptide.setSequence(seq);
			ret.addAll(indexStore.getProteins(indexPeptide));
		}

		return ret;
	}

	/**
	 * The DBIndexer builds an index from a Fasta file in the same directory were
	 * the fasta file was located. This function will check if the entered parameter
	 * file is located in the default dbindex location (defined by a property in the
	 * dbindex.properties). Copy the file if it was not located there, and return
	 * the new File.
	 *
	 * @param fastaFile
	 * @return
	 */
	public static File moveFastaToIndexLocation(File fastaFile) {
		try {
			final String dbindexPath = DBIndexImpl.getDBIndexPath();
			final File dbIndexFolder = new File(dbindexPath);
			if (!fastaFile.getParentFile().equals(dbIndexFolder)) {
				final File newFile = new File(dbIndexFolder.getAbsolutePath() + File.separator
						+ FilenameUtils.getName(fastaFile.getAbsolutePath()));
				if (newFile.exists())
					return newFile;
				FileUtils.copyFileToDirectory(fastaFile, dbIndexFolder);

				if (newFile.exists())
					return newFile;
			}
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return fastaFile;
	}

	public long getNumParentMasses() throws DBIndexStoreException {
		return indexStore.getNumberSequences();

	}

	public List<Double> getParentMasses() throws DBIndexStoreException {
		final List<Integer> entryKeys = indexStore.getEntryKeys();
		final List<Double> ret = new ArrayList<Double>();
		for (final Integer entryKey : entryKeys) {
			ret.add(entryKey * 1.0 / sparam.getMassGroupFactor());
		}
		return ret;

	}
}
