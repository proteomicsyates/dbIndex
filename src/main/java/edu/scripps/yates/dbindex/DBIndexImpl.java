package edu.scripps.yates.dbindex;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.dbindex.DBIndexer.IndexerMode;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;
import edu.scripps.yates.dbindex.io.SearchParamReader;
import edu.scripps.yates.dbindex.util.PropertiesReader;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexInterface;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexStoreException;
import edu.scripps.yates.utilities.fasta.dbindex.IndexType;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.fasta.dbindex.MassRange;
import edu.scripps.yates.utilities.fasta.dbindex.PeptideFilter;
import edu.scripps.yates.utilities.masses.AssignMass;
import gnu.trove.map.hash.THashMap;

public class DBIndexImpl implements DBIndexInterface {
	private final static Logger log = Logger.getLogger(DBIndexImpl.class);
	private static String dbIndexPath;
	protected DBIndexer indexer;
	private static final String PINT_DEVELOPER_ENV_VAR = "PINT_DEVELOPER";

	private final Map<String, Set<IndexedProtein>> proteinsByPeptideSeqs = new THashMap<>();
	protected final static Map<File, DBIndexImpl> dbIndexByFile = new THashMap<>();
	protected final static Map<String, DBIndexImpl> dbIndexByParamKey = new THashMap<>();

	public static DBIndexImpl getByParamFile(File paramFile) {
		if (dbIndexByFile.containsKey(paramFile)) {
			return dbIndexByFile.get(paramFile);
		}
		return new DBIndexImpl(paramFile);
	}

	public static DBIndexImpl getByParam(DBIndexSearchParams sParam) {
		if (dbIndexByParamKey.containsKey(sParam.getFullIndexFileName())) {
			return dbIndexByParamKey.get(sParam.getFullIndexFileName());
		}
		return new DBIndexImpl(sParam);
	}

	public DBIndexImpl() {

	}

	/**
	 * Load the index of a fasta database stated on the paramFile. this will
	 * construct the index if it is not previously build.
	 *
	 * @param paramFile
	 */
	public DBIndexImpl(File paramFile) {
		final String paramFileName = SearchParamReader.DEFAULT_PARAM_FILE_NAME;

		try {

			final String dbIndexPath = getDBIndexPath();
			SearchParamReader pr = null;
			if (paramFile != null) {
				pr = new SearchParamReader(paramFile);
			} else {
				pr = new SearchParamReader(dbIndexPath, paramFileName);
			}

			final SearchParams sParam = pr.getSearchParams();
			final edu.scripps.yates.dbindex.DBIndexer.IndexerMode indexerMode = sParam.isUseIndex()
					? IndexerMode.SEARCH_INDEXED
					: IndexerMode.SEARCH_UNINDEXED;

			indexer = new DBIndexer(sParam, indexerMode);
			try {
				indexer.init();

			} catch (final DBIndexerException ex) {
				indexer = new DBIndexer(sParam, IndexerMode.INDEX);
				try {
					indexer.init();
					indexer.run();
				} catch (final DBIndexerException e) {
					e.printStackTrace();
					log.error("Could not initialize the indexer in search mode and init the worker thread");
				}
			}
			dbIndexByFile.put(paramFile, this);
		} catch (final IOException e) {
			e.printStackTrace();
			log.error(e.getMessage());

		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
		}

	}

	public DBIndexImpl(DBIndexSearchParams sParam) {

		this(sParam, sParam.isUseIndex() ? IndexerMode.SEARCH_INDEXED : IndexerMode.SEARCH_UNINDEXED);

	}

	/**
	 * Load the index of a fasta database stated on the paramFile. this will
	 * construct the index if it is not previously build.
	 *
	 * @param paramFile
	 */
	public DBIndexImpl(DBIndexSearchParams sParam, edu.scripps.yates.dbindex.DBIndexer.IndexerMode indexerMode) {
		try {

			// init the masses
			AssignMass.getInstance(sParam.isUseMonoParent());

			indexer = new DBIndexer(sParam, indexerMode);
			try {
				indexer.init();

			} catch (final DBIndexerException ex) {
				indexer = new DBIndexer(sParam, IndexerMode.INDEX);
				try {
					indexer.init();
					indexer.run();
				} catch (final DBIndexerException e) {
					e.printStackTrace();
					log.error("Could not initialize the indexer in search mode and init the worker thread");
				}
			}
			dbIndexByParamKey.put(sParam.getFullIndexFileName(), this);
		} catch (final Exception e) {
			e.printStackTrace();
			log.error(e.getMessage());
			throw e;
		}

	}

	public static String getDBIndexPath() {
		if (dbIndexPath == null) {
			final Map<String, String> env = System.getenv();
			String dbIndexPathProperty = PropertiesReader.DBINDEX_PATH_SERVER;
			if (env.containsKey(PINT_DEVELOPER_ENV_VAR) && env.get(PINT_DEVELOPER_ENV_VAR).equals("true")) {
				log.info("USING DEVELOPMENT MODE");
				dbIndexPathProperty = PropertiesReader.DBINDEX_PATH;
			}
			log.info("Using property: " + dbIndexPathProperty);
			try {
				dbIndexPath = (String) PropertiesReader.getProperties().get(dbIndexPathProperty);
				log.info("Using local folder: " + dbIndexPath);
			} catch (final Exception e) {
				e.printStackTrace();
			}
		}
		return dbIndexPath;
	}

	/**
	 * Get list of sequences in index for the mass, and using tolerance in
	 * SearchParams. Mass dependant tolerance is already scaled/precalculated
	 *
	 * Requires call to init() and run() first, which might perform indexing, if
	 * index for current set of sequest params does not exist
	 *
	 * @param precursorMass mass to select by, in ppm
	 * @param massTolerance mass tolerance, already calculated for that mass as it
	 *                      is mass-dependant
	 * @return list of matching sequences
	 * @throws DBIndexStoreException
	 */
	@Override
	public List<IndexedSequence> getSequences(double precursorMass, double massTolerance) throws DBIndexStoreException {
		return indexer.getSequencesUsingDaltonTolerance(precursorMass, massTolerance);
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
	@Override
	public List<IndexedSequence> getSequences(List<MassRange> massRanges) throws DBIndexStoreException {
		return indexer.getSequences(massRanges);
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
	@Override
	public List<IndexedProtein> getProteins(IndexedSequence seq) throws DBIndexStoreException {
		return indexer.getProteins(seq);
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
	@Override
	public Set<IndexedProtein> getProteins(String seq) throws DBIndexStoreException {
		// look into the cached proteins
		if (proteinsByPeptideSeqs.containsKey(seq) && !proteinsByPeptideSeqs.get(seq).isEmpty())
			return proteinsByPeptideSeqs.get(seq);
		final Set<IndexedProtein> proteins = indexer.getProteins(seq);
		if (proteins != null && !proteins.isEmpty()) {
			proteinsByPeptideSeqs.put(seq, proteins);
		}
		// log.debug(proteins.size() + " proteins contains peptide '" + seq +
		// "'");
		return proteins;
	}

	/**
	 * Gets the default search params.<br>
	 * <b>Note that this parameters will not be valid for crosslinker analysis,
	 * since H2O+PROTON mass is added to any peptide mass in the index.
	 *
	 * @param fastaFile
	 * @return
	 */
	public static DBIndexSearchParams getDefaultDBIndexParams(File fastaFile) {
		return getDefaultDBIndexParams(fastaFile.getAbsolutePath());
	}

	/**
	 * Gets the default search params.<br>
	 * <b>Note that this parameters will not be valid for crosslinker analysis,
	 * since H2O+PROTON mass is added to any peptide mass in the index.
	 *
	 * @param fastaFile
	 * @param inMemoryIndex overrides the default property of "in memory index"
	 * @return
	 */
	public static DBIndexSearchParams getDefaultDBIndexParams(String fastaFilePath, boolean inMemoryIndex) {
		try {
			final IndexType indexType = IndexType.valueOf(
					String.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_INDEX_TYPE)));

			final int indexFactor = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_INDEX_FACTOR));
			final String dataBaseName = fastaFilePath;
			final int maxMissedCleavages = Integer.valueOf(PropertiesReader.getProperties()
					.getProperty(PropertiesReader.DEFAULT_MAX_INTERNAL_CLEAVAGES_SITES));
			final double maxPrecursorMass = Double
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MAX_PRECURSOR_MASS));
			final double minPrecursorMass = Double
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MIN_PRECURSOR_MASS));
			final boolean useIndex = Boolean
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.USE_INDEX));
			final String enzymeNocutResidues = String
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.ENZYME_NOCUT_RESIDUES));

			final String enzymeResidues = String
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.ENZYME_RESIDUES));

			final int enzymeOffset = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.ENZYME_OFFSET));

			final boolean useMono = Boolean
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MASS_TYPE_PARENT));

			final boolean isH2OPlusProtonAdded = Boolean
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.ADD_H2O_PLUS_PROTON));
			final int massGroupFactor = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MASS_GROUP_FACTOR));

			final String uniprotVersion = String
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.UNIPROT_VERSION));

			final DBIndexSearchParamsImpl ret = new DBIndexSearchParamsImpl(indexType, inMemoryIndex, indexFactor,
					dataBaseName, maxMissedCleavages, maxPrecursorMass, minPrecursorMass, useIndex, enzymeNocutResidues,
					enzymeResidues, enzymeOffset, useMono, isH2OPlusProtonAdded, massGroupFactor, null, false, null,
					uniprotVersion, null);
			return ret;
		} catch (final NumberFormatException e) {
			e.printStackTrace();
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Gets the default search params.<br>
	 * <b>Note that this parameters will not be valid for crosslinker analysis,
	 * since H2O+PROTON mass is added to any peptide mass in the index.
	 *
	 * @param fastaFile
	 * @return
	 */
	public static DBIndexSearchParams getDefaultDBIndexParams(String fastaFilePath) {
		boolean inMemoryIndex;
		try {
			final Map<String, String> env = System.getenv();
			if (env.containsKey(PINT_DEVELOPER_ENV_VAR) && env.get(PINT_DEVELOPER_ENV_VAR).equals("true")) {
				log.info("USING DEVELOPMENT MODE");

				inMemoryIndex = Boolean.valueOf(
						PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_IN_MEMORY_INDEX));
			} else {
				inMemoryIndex = Boolean.valueOf(
						PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_IN_MEMORY_INDEX));
			}
			log.info("InMemoryIndex=" + inMemoryIndex);
			return getDefaultDBIndexParams(fastaFilePath, inMemoryIndex);
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Gets the default search params.<br>
	 * <b>Note that this parameters will not be valid for crosslinker analysis,
	 * since H2O+PROTON mass is added to any peptide mass in the index.
	 *
	 * @param fastaFile
	 * @return
	 */
	public static DBIndexSearchParams getDefaultDBIndexParamsForCrosslinkerAnalysis(File fastaFile) {
		return getDefaultDBIndexParamsForCrosslinkerAnalysis(fastaFile.getAbsolutePath());
	}

	/**
	 * Gets the default search params.<br>
	 * <b>Note that this parameters will not be valid for crosslinker analysis,
	 * since H2O+PROTON mass is added to any peptide mass in the index.
	 *
	 * @param fastaFile
	 * @return
	 */
	public static DBIndexSearchParams getDefaultDBIndexParamsForCrosslinkerAnalysis(String fastaFilePath) {
		boolean inMemoryIndex;
		try {
			final Map<String, String> env = System.getenv();
			if (env.containsKey(PINT_DEVELOPER_ENV_VAR) && env.get(PINT_DEVELOPER_ENV_VAR).equals("true")) {
				log.info("USING DEVELOPMENT MODE");
				inMemoryIndex = false;

			} else {
				inMemoryIndex = Boolean.valueOf(
						PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_IN_MEMORY_INDEX));
			}
			log.info("InMemoryIndex=" + inMemoryIndex);
			return getDefaultDBIndexParamsForCrosslinkerAnalysis(fastaFilePath, inMemoryIndex);
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Gets the {@link DBIndexSearchParams} needed for a proteoform analysis,
	 * getting annotations of sequence variances, mutations, splice variants and
	 * PTMs from uniprot and including them in the index.
	 * 
	 * @param fastaFilePath
	 * @param maxMissedCleavages
	 * @param minPrecursorMass
	 * @param maxPrecursorMass
	 * @param enzymeResidues
	 * @param enzymeNocutResidues
	 * @param enzymeOffset
	 * @param semiCleavage
	 * @param uniprotVersion
	 * @param discardDecoyRegexp
	 * @param uniprotReleasesFolderPath
	 * @return
	 */
	public static DBIndexSearchParams getDefaultDBIndexParamsForProteoformAnalysis(String fastaFilePath,
			int maxMissedCleavages, double minPrecursorMass, double maxPrecursorMass, String enzymeResidues,
			String enzymeNocutResidues, int enzymeOffset, boolean semiCleavage, String uniprotVersion,
			String discardDecoyRegexp, String uniprotReleasesFolderPath) {
		try {
			final IndexType indexType = IndexType.valueOf(
					String.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_INDEX_TYPE)));

			final int indexFactor = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_INDEX_FACTOR));
			final String dataBaseName = fastaFilePath;

			final boolean useIndex = Boolean
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.USE_INDEX));

			final boolean useMono = Boolean
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MASS_TYPE_PARENT));
			final boolean isH2OPlusProtonAdded = true;
			final int massGroupFactor = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MASS_GROUP_FACTOR));

			final char[] mandatoryInternalAAs = null;

			final boolean inMemoryIndex = false;
			final PeptideFilter peptideFilter = null;

			final DBIndexSearchParamsImpl ret = new DBIndexSearchParamsImpl(indexType, inMemoryIndex, indexFactor,
					dataBaseName, maxMissedCleavages, maxPrecursorMass, minPrecursorMass, useIndex, enzymeNocutResidues,
					enzymeResidues, enzymeOffset, useMono, isH2OPlusProtonAdded, massGroupFactor, mandatoryInternalAAs,
					semiCleavage, peptideFilter, uniprotVersion, discardDecoyRegexp);
			ret.setUniprotReleasesFolder(new File(uniprotReleasesFolderPath));
			ret.setLookProteoforms(true);
			return ret;
		} catch (final NumberFormatException e) {
			e.printStackTrace();
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * Gets the default search params.<br>
	 * <b>Note that this parameters will not be valid for crosslinker analysis,
	 * since H2O+PROTON mass is added to any peptide mass in the index.
	 *
	 * @param fastaFile
	 * @param inMemoryIndex overrides the property inMemoryIndex
	 * @return
	 */
	public static DBIndexSearchParams getDefaultDBIndexParamsForCrosslinkerAnalysis(String fastaFilePath,
			boolean inMemoryIndex) {
		try {
			final IndexType indexType = IndexType.valueOf(
					String.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_INDEX_TYPE)));

			final int indexFactor = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.DEFAULT_INDEX_FACTOR));
			final String dataBaseName = fastaFilePath;
			final int maxMissedCleavages = Integer.valueOf(PropertiesReader.getProperties()
					.getProperty(PropertiesReader.DEFAULT_MAX_INTERNAL_CLEAVAGES_SITES));
			final double maxPrecursorMass = Double
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MAX_PRECURSOR_MASS));
			final double minPrecursorMass = Double
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MIN_PRECURSOR_MASS));
			final boolean useIndex = Boolean
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.USE_INDEX));
			final String enzymeNocutResidues = String
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.ENZYME_NOCUT_RESIDUES));

			final String enzymeResidues = String
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.ENZYME_RESIDUES));

			final int enzymeOffset = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.ENZYME_OFFSET));

			final boolean useMono = Boolean
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MASS_TYPE_PARENT));
			// SET TO FALSE TO ALLOW CROSSLINKER ANALYSES
			final boolean isH2OPlusProtonAdded = false;
			final int massGroupFactor = Integer
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.MASS_GROUP_FACTOR));

			final char[] mandatoryInternalAAs = PropertiesReader.getProperties()
					.getProperty(PropertiesReader.MANDATORY_INTERNAL_AAs).toCharArray();
			final String uniprotVersion = String
					.valueOf(PropertiesReader.getProperties().getProperty(PropertiesReader.UNIPROT_VERSION));
			final DBIndexSearchParamsImpl ret = new DBIndexSearchParamsImpl(indexType, inMemoryIndex, indexFactor,
					dataBaseName, maxMissedCleavages, maxPrecursorMass, minPrecursorMass, useIndex, enzymeNocutResidues,
					enzymeResidues, enzymeOffset, useMono, isH2OPlusProtonAdded, massGroupFactor, mandatoryInternalAAs,
					false, null, uniprotVersion, null);
			return ret;
		} catch (final NumberFormatException e) {
			e.printStackTrace();
		} catch (final Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * @return the indexer
	 */
	public DBIndexer getIndexer() {
		return indexer;
	}

	@Override
	public IndexedProtein getIndexedProteinById(int proteinId) {
		final String proteinDef = indexer.getProteinCache().getProteinDef(proteinId);
		if (proteinDef != null) {
			final IndexedProtein ip = new IndexedProtein(proteinDef, proteinId);
			return ip;
		}
		return null;
	}

	@Override
	public String getProteinSequenceById(int proteinId) {
		return indexer.getProteinCache().getProteinSequence(proteinId);
	}

}
