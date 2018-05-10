package edu.scripps.yates.dbindex;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.dbindex.model.DBIndexSearchParams;
import edu.scripps.yates.dbindex.model.ResidueInfo;
import edu.scripps.yates.dbindex.util.IndexUtil;
import edu.scripps.yates.utilities.bytes.DynByteBuffer;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;

/**
 * Simple table implementation with: blazmass_sequences integer mass (PRI KEY,
 * ASC and DESC index for quick range lookup) byte data : encoded peptide info
 * (prot id, offset, length)
 *
 * Meant to be used by DBIndexStoreSQLiteMult only, internal private class
 *
 * Caches sequences in intermediate structure and commits to db every 200k
 * sequences.
 *
 */
public class DBIndexStoreSQLiteByte extends DBIndexStoreSQLiteAbstract {

	private PreparedStatement addSeqStatement;
	protected PreparedStatement updateSeqStatement;
	private PreparedStatement updateAppendSeqStatement;
	protected PreparedStatement getSeqStatement;
	protected PreparedStatement[] getSeqMassRanges2_Statements;
	private PreparedStatement getSeqDataStatement;
	private PreparedStatement getSeqExistsStatement;
	private final String blazmass_sequences_pri_key = "precursor_mass_key";
	// changed by Salva 11Nov2014, using the value on IndexUtil
	// private final DynByteBuffer data[] = new DynByteBuffer[MAX_MASS];
	// private final DynByteBuffer data[] = new DynByteBuffer[MAX_MASS
	// * IndexUtil.getMassGroupFactor()];
	protected final TIntObjectHashMap<DynByteBuffer> dataMap = new TIntObjectHashMap<DynByteBuffer>();

	// we really only need MAX_MASS / NUM_BUCKETS, but need to shift

	private final int bucketId;
	// single sequence commit mode coutners
	private int curCommit = 0;
	private static final int SEQ_UPDATE_INTERVAL = 20000;
	// seq count since last commit of entire cache
	protected int seqCount = 0;
	// total sequence count indexed
	protected int totalSeqCount = 0;
	private final TIntHashSet massKeys = new TIntHashSet();
	// full cache commit interval
	// NOTE: optimized for 7GB heap space, adjsut this setting accordingly
	// private static final int FULL_CACHE_COMMIT_INTERVAL = (100 * 1000 * 1000)
	// / Constants.NUM_BUCKETS;
	protected static final int FULL_CACHE_COMMIT_INTERVAL = (1000000) / Constants.NUM_BUCKETS;

	protected DBIndexStoreSQLiteByte(int bucketId, ProteinCache proteinCache) {
		super(proteinCache);
		this.bucketId = bucketId;
	}

	protected DBIndexStoreSQLiteByte(boolean inMemory, int bucketId, ProteinCache proteinCache) {
		super(SearchParams.getInstance(), inMemory, proteinCache);
		this.bucketId = bucketId;
	}

	protected DBIndexStoreSQLiteByte(DBIndexSearchParams searchParams, boolean inMemory, int bucketId,
			ProteinCache proteinCache) {
		super(searchParams, inMemory, proteinCache);
		this.bucketId = bucketId;
	}

	@Override
	protected void initStatementsCommon() throws DBIndexStoreException {
	}

	@Override
	protected void initStatements() throws DBIndexStoreException {
		try {

			addSeqStatement = con.prepareStatement(
					"INSERT INTO " + getIndexTableName() + " (precursor_mass_key, data) " + "VALUES (?, ?);");

			updateSeqStatement = con.prepareStatement(
					"UPDATE " + getIndexTableName() + " SET data = ? " + "WHERE precursor_mass_key = ?;");

			updateAppendSeqStatement = con.prepareStatement(
					"UPDATE " + getIndexTableName() + " SET data = (data || ?) " + "WHERE precursor_mass_key = ?;");

			getSeqStatement = con.prepareStatement("SELECT precursor_mass_key, data " + "FROM " + getIndexTableName()
					+ " " + "WHERE precursor_mass_key BETWEEN ? AND ?;");

			getSeqMassRanges2_Statements = new PreparedStatement[Constants.MAX_MASS_RANGES];

			final StringBuilder stSb = new StringBuilder();
			stSb.append("SELECT DISTINCT precursor_mass_key, data ");
			stSb.append("FROM ").append(getIndexTableName()).append(" WHERE ");
			for (int st = 0; st < Constants.MAX_MASS_RANGES; ++st) {
				if (st > 0) {
					stSb.append(" OR ");
				}
				stSb.append("precursor_mass_key BETWEEN ? AND ? ");
				getSeqMassRanges2_Statements[st] = con.prepareStatement(stSb.toString());
			}

			getSeqDataStatement = con.prepareStatement(
					"SELECT data " + "FROM " + getIndexTableName() + " " + "WHERE precursor_mass_key = ?;");

			getSeqExistsStatement = con.prepareStatement("SELECT precursor_mass_key " + "FROM " + getIndexTableName()
					+ " " + "WHERE precursor_mass_key = ? LIMIT 1;");

		} catch (final SQLException e) {
			logger.error("Error initializing statements in db, path: " + dbPath, e);
			throw new DBIndexStoreException("Error initializing statements in db, path: " + dbPath, e);
		}
	}

	@Override
	protected String getSequencesTablePriKey() {
		return blazmass_sequences_pri_key;
	}

	/**
	 * Get sequence data for long representation of precursor mass or null if
	 * not present
	 *
	 * @param massKey
	 * @return encoded sequence data string for multiple peptide, or null
	 * @throws SQLException
	 */
	private byte[] getSequenceData(int massKey) throws SQLException {
		if (massKey == 1952)
			System.out.println("a ver");
		byte[] ret = null;
		getSeqDataStatement.setInt(1, massKey);
		final ResultSet rs = getSeqDataStatement.executeQuery();
		if (rs.next()) {
			ret = rs.getBytes(1);
		}
		rs.close();

		return ret;
	}

	/**
	 * check if sequence exists by mass
	 *
	 * @param massKey
	 * @return
	 * @throws SQLException
	 */
	private boolean getSequenceExists(int massKey) throws SQLException {
		// instead of accessing the DB, keep a HashSet in memory with the
		// massKeys
		return massKeys.contains(massKey);
		/*
		 * boolean ret = false; getSeqExistsStatement.setInt(1, massKey); final
		 * ResultSet rs = getSeqExistsStatement.executeQuery(); if (rs.next()) {
		 * ret = true; } rs.close();
		 * 
		 * return ret;
		 */
	}

	/**
	 * Add the new sequence to in memory cached pre-commit sequence hash
	 *
	 * @param precMass
	 * @param seqOffset
	 * @param seqLength
	 * @param proteinId
	 */
	protected void updateCachedData(double precMass, int seqOffset, int seqLength, int proteinId) {
		// changed by Salva 11Nov2014, using the value on IndexUtil
		final int rowId = (int) (precMass * sparam.getMassGroupFactor());

		// change by Salva 21Nov2014
		// DynByteBuffer byteBuffer = data[rowId];
		DynByteBuffer byteBuffer = dataMap.get(rowId);
		if (byteBuffer == null) {
			byteBuffer = new DynByteBuffer();
			// change by Salva 21Nov2014
			// data[rowId] = byteBuffer;
			dataMap.put(rowId, byteBuffer);
		}

		// store in Little-Endian order
		final byte[] seqMassB = DynByteBuffer.toByteArray(precMass);
		byteBuffer.add(seqMassB);

		final byte[] seqOffsetB = DynByteBuffer.toByteArray(seqOffset);
		byteBuffer.add(seqOffsetB);

		final byte[] seqLengthB = DynByteBuffer.toByteArray(seqLength);
		byteBuffer.add(seqLengthB);

		final byte[] proteinIdB = DynByteBuffer.toByteArray(proteinId);
		byteBuffer.add(proteinIdB);

		// commit this mass after COMMIT_SEQUENCES sequences
		if (true) {
			if (byteBuffer.getSize() > Constants.COMMIT_SEQUENCES * Constants.BYTE_PER_SEQUENCE) {
				try {
					insertSequence(rowId, byteBuffer);
					// clear since we wrote it to db
					byteBuffer.clear();
				} catch (final SQLException e) {
					logger.error("Error commiting mass buffer");
				}
			}
		}
	}

	/**
	 * updates db rows with the sequence, optionally commits it
	 *
	 * @param massKey
	 * @param buf
	 * @throws SQLException
	 */
	protected void insertSequence(int massKey, DynByteBuffer buf) throws SQLException {

		// logger.info( bucketId +
		// ": Commiting cached sequence data for rowId: " + rowId);
		if (getSequenceExists(massKey)) {
			// merge, bytes in db with the bytes in cache
			// build up updated string

			updateAppendSeqStatement.setBytes(1, buf.getData());
			updateAppendSeqStatement.setInt(2, massKey);
			updateAppendSeqStatement.executeUpdate();
		} else {
			// write a brand new row for the mass

			// insert
			addSeqStatement.setInt(1, massKey);
			addSeqStatement.setBytes(2, buf.getData());
			addSeqStatement.executeUpdate();
			// keep mass key
			massKeys.add(massKey);
		}

		if (curCommit++ == SEQ_UPDATE_INTERVAL) {
			curCommit = 0;
			// con.commit(); //do not db commit every sequence for now, let
			// the big cache commit handle it every
			// FULL_CACHE_COMMIT_INTERVAL
		}

	}

	@Override
	protected final void commitCachedData() throws SQLException {
		// called at the end of indexing to commit all buffers
		logger.debug("Bucket " + bucketId + ": Commiting cached data");
		// change by Salva 21Nov2014
		// for (int massKey = 0; massKey < MAX_MASS; ++massKey) {
		// for (int massKey = 0; massKey < data.length; ++massKey) {
		for (final Integer massKey : dataMap.keys()) {

			// see if already in index
			// change by Salva 21Nov2014
			// DynByteBuffer cached = data[massKey];
			final DynByteBuffer cached = dataMap.get(massKey);
			if (cached == null || cached.getSize() == 0) {
				continue;
			}
			if (getSequenceExists(massKey)) {
				// merge, bytes in db with the bytes in cache
				// build up updated string

				updateAppendSeqStatement.setBytes(1, cached.getData());
				updateAppendSeqStatement.setInt(2, massKey);
				updateAppendSeqStatement.executeUpdate();
			} else {
				// write a brand new row for the mass

				// insert
				addSeqStatement.setInt(1, massKey);
				addSeqStatement.setBytes(2, cached.getData());
				addSeqStatement.executeUpdate();
				// keep mass key
				massKeys.add(massKey);
			}
			cached.clear();
		}

		// clear cached data
		// salva change 21Nov2014
		// Arrays.fill(data, null);
		dataMap.clear();

		// logger.info( "Commit start");
		con.commit();
		// logger.info( "Commit end");

		// logger.info( bucketId +
		// ": Done commiting cached data");

	}

	@Override
	public FilterResult filterSequence(double precMass, String sequence) {
		// no filtering, we add every sequence
		return FilterResult.INCLUDE;
	}

	@Override
	public void addSequence(double precMass, int seqOffset, int seqLength, String sequence, String resLeft,
			String resRight, long proteinId) throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}
		if ("YFDRDDVALKNF".equals(sequence)) {
			System.out.println(sequence + "\t" + precMass);
		}
		totalSeqCount++;

		updateCachedData(precMass, (short) seqOffset, (short) seqLength, (int) proteinId);

		if (true) {
			// despite commiting each sequence, commit and clear all after
			// 100mil
			// to free up memory, optimize buffers, and ensure a commit
			if (seqCount++ == FULL_CACHE_COMMIT_INTERVAL) {
				// logger.info( "Commiting all");
				seqCount = 0;
				try {
					commitCachedData();
				} catch (final SQLException ex) {
					logger.error("Error commiting cached data", ex);
				}
			}
		}

	}

	@Override
	public ResidueInfo getResidues(IndexedSequence peptideSequence, IndexedProtein protein)
			throws DBIndexStoreException {
		// implemented by outter class
		return null;

	}

	@Override
	public List<IndexedSequence> getSequences(double precMass, double tolerance) throws DBIndexStoreException {

		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		// logger.log(Level.FINE, "Starting peptide sequences query");
		// long start = System.currentTimeMillis();

		final double toleranceFactor = tolerance; // precMass * tolerance;
		double minMassF = precMass - toleranceFactor;
		if (minMassF < 0f) {
			minMassF = 0;
		}
		final double maxMassF = precMass + toleranceFactor;
		// changed by Salva 11Nov2014, using the value on IndexUtil
		final int precMassInt = (int) (precMass * sparam.getMassGroupFactor());
		// long toleranceInt = (long) (tolerance * MASS_STORE_MULT);
		// changed by Salva 11Nov2014, using the value on IndexUtil
		final int toleranceInt = (int) (toleranceFactor * sparam.getMassGroupFactor()); // Robin
		// fixed
		// the
		// ppm
		// calculation

		int minMass = precMassInt - toleranceInt;
		// change by Salva 21Nov2014
		minMass = (int) (minMassF * sparam.getMassGroupFactor());
		// System.out.println(precMassInt + " " + toleranceInt);

		if (minMass < 0) {
			minMass = 0;
		}
		int maxMass = precMassInt + toleranceInt;
		// change by Salva 21Nov2014
		maxMass = (int) (maxMassF * sparam.getMassGroupFactor());
		final List<IndexedSequence> ret = new ArrayList<IndexedSequence>();

		ResultSet rs = null;
		try {
			// changed by Salva 11Nov2014, using the value on IndexUtil
			getSeqStatement.setInt(1, minMass - IndexUtil.getNumRowsToLookup());// PPM_PER_ENTRY);
			// changed by Salva 11Nov2014, using the value on IndexUtil
			getSeqStatement.setInt(2, maxMass + IndexUtil.getNumRowsToLookup());// PPM_PER_ENTRY);
			rs = getSeqStatement.executeQuery();

			while (rs.next()) {
				// final int massKey = rs.getInt(1);
				final byte[] seqData = rs.getBytes(2);

				parseAddPeptideInfo(seqData, ret, minMassF, maxMassF);

			}

		} catch (final SQLException e) {
			final String msg = "Error getting peptides ";
			logger.error(msg, e);
			throw new DBIndexStoreException(msg, e);
		} finally {
			if (rs != null) {
				try {
					rs.close();
				} catch (final SQLException ex) {
					logger.error(null, ex);
				}
			}
		}

		// long end = System.currentTimeMillis();
		// logger.info( "Peptide sequences query took: " + (end -
		// start) + " milisecs, results: " + ret.size());

		return ret;
	}

	@Override
	public List<IndexedSequence> getSequences(List<MassRange> ranges) throws DBIndexStoreException {
		if (ranges.size() == 1) {
			final MassRange range = ranges.get(0);
			return getSequences(range.getPrecMass(), range.getTolerance());
		}

		throw new DBIndexStoreException("Multiple mass ranges not implemented yet!");
	}

	/**
	 * Parse the encoded sequences and return them in toInsert
	 *
	 * @param data
	 * @param toInsert
	 * @param minMass
	 * @param maxMass
	 */
	protected void parseAddPeptideInfo(byte[] data, List<IndexedSequence> toInsert, double minMass, double maxMass) {

		// to collapse multiple sequences into single one, with multproteins
		final Map<String, List<IndexedSeqInternal>> temp = new THashMap<String, List<IndexedSeqInternal>>();

		final int dataLength = data.length;

		// if (dataLength % 4 != 0) {
		// throw new RuntimeException("Unexpected number of peptide items: "
		// + dataLength);
		// }

		for (int i = 0; i < dataLength; i += Constants.BYTE_PER_SEQUENCE) {
			final int first = i + 8; // double mass
			final int second = i + 12;
			final int third = i + 16;
			final int fourth = i + 20;

			byte[] slice = Arrays.copyOfRange(data, i, first);
			final double seqMass = DynByteBuffer.toDouble(slice);

			if (seqMass < minMass || seqMass > maxMass) {
				// skip the sequence, it qualified the bucket, but not the
				// actual mass
				continue;
			}

			slice = Arrays.copyOfRange(data, first, second);
			final int offset = DynByteBuffer.toInt(slice);
			slice = Arrays.copyOfRange(data, second, third);
			final int length = DynByteBuffer.toInt(slice);
			slice = Arrays.copyOfRange(data, third, fourth);
			final int proteinId = DynByteBuffer.toInt(slice);

			String peptideSequence = null;
			// System.out.print("prot id: " + proteinId + " offset: " +
			// offset + " len: " + length);

			peptideSequence = proteinCache.getPeptideSequence(proteinId, offset, length);
			// System.out.println("pep: " + peptideSequence);

			// we are cheating and supplying protein id instead of peptide
			// id
			// to set it temporarily, before we merge them into a list
			final IndexedSeqInternal tempSequence = new IndexedSeqInternal(seqMass, offset, length, proteinId,
					peptideSequence);

			List<IndexedSeqInternal> sequences = temp.get(peptideSequence);
			if (sequences == null) {
				sequences = new ArrayList<IndexedSeqInternal>();
				temp.put(peptideSequence, sequences);
			}
			sequences.add(tempSequence);

		}

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
			final String protSequence = proteinCache.getProteinSequence(firstSeq.proteinId);
			final ResidueInfo residues = Util.getResidues(null, firstSeq.offset, firstSeq.length, protSequence);
			mergedSequence.setResidues(residues);
			toInsert.add(mergedSequence);
		}

	}

	@Override
	public boolean supportsProteinCache() {
		return true;
	}

	@Override
	public void setProteinCache(ProteinCache proteinCache) {
		this.proteinCache = proteinCache;
	}

	@Override
	public long addProteinDef(long num, String definition, String proteinSequence) throws DBIndexStoreException {
		// implemented by outter class
		return num;

	}

	@Override
	public List<IndexedProtein> getProteins(IndexedSequence sequence) throws DBIndexStoreException {
		// implemented by outter class
		return null;

	}

	@Override
	protected void createTempIndex() {
		// not needed
	}

	@Override
	protected void deleteTempIndex() {
		// not needed
	}

	protected String getIndexTableName() {
		return "blazmass_sequences";
	}

	@Override
	protected void createTables() throws DBIndexStoreException {
		try {

			executeStatement("CREATE TABLE IF NOT EXISTS " + getIndexTableName() + " "
					+ "(precursor_mass_key INTEGER PRIMARY KEY, " + "data BINARY" + ");");

		} catch (final SQLException e) {
			logger.error("Error creating tables, ", e);
			throw new DBIndexStoreException("Error creating tables. ", e);
		}
	}

	@Override
	protected void createIndex() throws DBIndexStoreException {
		try {
			// for fast lookup by mass range
			// logger.info(
			// "Creating precursor_mass_key_index_asc ");
			// executeStatement("CREATE INDEX IF NOT EXISTS
			// precursor_mass_key_index_asc ON blazmass_sequences
			// (precursor_mass_key ASC);");
			// dsc index speeds up the range queries another 10x more
			// logger.info(
			// "Creating precursor_mass_key_index_dsc ");
			executeStatement("CREATE INDEX IF NOT EXISTS precursor_mass_key_index_dsc ON " + getIndexTableName()
					+ " (precursor_mass_key DESC);");
		} catch (final SQLException e) {
			e.printStackTrace();
			logger.error("Error creating index, ", e);
			throw new DBIndexStoreException("Error creating index, ", e);
		}
	}

	@Override
	protected void closeConnection() {
		try {

			if (addSeqStatement != null) {
				addSeqStatement.close();
				addSeqStatement = null;
			}

			if (updateSeqStatement != null) {
				updateSeqStatement.close();
				updateSeqStatement = null;
			}

			if (updateAppendSeqStatement != null) {
				updateAppendSeqStatement.close();
				updateAppendSeqStatement = null;
			}

			if (getSeqStatement != null) {
				getSeqStatement.close();
				getSeqStatement = null;
			}

			for (int i = 0; i < Constants.MAX_MASS_RANGES; ++i) {
				if (getSeqMassRanges2_Statements[i] != null) {
					getSeqMassRanges2_Statements[i].close();
					getSeqMassRanges2_Statements[i] = null;
				}
			}

			if (getSeqDataStatement != null) {
				getSeqDataStatement.close();
				getSeqDataStatement = null;
			}

			if (getSeqExistsStatement != null) {
				getSeqExistsStatement.close();
				getSeqExistsStatement = null;
			}

			if (con != null) {
				con.close();
				con = null;
			}
		} catch (final SQLException ex) {
			logger.debug("Error closing SQLite connection", ex);
		}
	}

	@Override
	public long getNumberSequences() throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		ResultSet rs = null;
		try {
			// we don't really now, return number of rows for now
			rs = executeQuery("SELECT COUNT(*) FROM " + getIndexTableName());
			return rs.getLong(1);
		} catch (final SQLException ex) {
			logger.error("Error executing query to get number of sequences.", ex);
			throw new DBIndexStoreException("Error executing query to get number of sequences.", ex);
		} finally {
			if (rs != null) {
				try {
					rs.close();
				} catch (final SQLException ex) {
					logger.error(null, ex);
				}
			}
		}

	}

	@Override
	public List<Integer> getEntryKeys() throws DBIndexStoreException {
		if (!inited) {
			throw new DBIndexStoreException("Indexer is not initialized");
		}

		ResultSet rs = null;
		try {
			// we don't really now, return number of rows for now
			rs = executeQuery("SELECT precursor_mass_key FROM " + getIndexTableName());
			final List<Integer> ret = new ArrayList<Integer>();
			while (rs.next()) {
				ret.add(rs.getInt(1));
			}
			return ret;
		} catch (final SQLException ex) {
			logger.error("Error executing query to get number of sequences.", ex);
			throw new DBIndexStoreException("Error executing query to get number of sequences.", ex);
		} finally {
			if (rs != null) {
				try {
					rs.close();
				} catch (final SQLException ex) {
					logger.error(null, ex);
				}
			}
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
