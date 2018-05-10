package edu.scripps.yates.dbindex.util;

import java.io.File;

import edu.scripps.yates.dbindex.DBIndexInterface;
import edu.scripps.yates.dbindex.io.DBIndexSearchParamsImpl;
import edu.scripps.yates.dbindex.model.DBIndexSearchParams;
import edu.scripps.yates.dbindex.model.PeptideFilter;

public class FastaDigestionConfiguration extends DigestionConfiguration {
	private final File fasta;
	private final boolean ignorePeptidesNotFoundInDB;

	public FastaDigestionConfiguration(File fasta, char[] enzymeArray, int numMisscleavages, boolean semiCleavage,
			String peptideFilterString, boolean ignorePeptidesNotFoundInDB) {
		super(enzymeArray, numMisscleavages, semiCleavage, peptideFilterString);
		this.fasta = fasta;
		this.ignorePeptidesNotFoundInDB = ignorePeptidesNotFoundInDB;
	}

	/**
	 * @return the fasta
	 */

	public File getFasta() {
		return fasta;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return super.toString() + " fasta: " + fasta.getAbsolutePath() + " ignorePeptidesNotFoundInDB: "
				+ ignorePeptidesNotFoundInDB;
	}

	public static DBIndexInterface getFastaDBIndex(FastaDigestionConfiguration fastadigestionConfiguration) {
		return getFastaDBIndex(fastadigestionConfiguration.getFasta(), fastadigestionConfiguration.getEnzymeArray(),
				fastadigestionConfiguration.getNumMisscleavages(), fastadigestionConfiguration.isSemiCleavage(),
				fastadigestionConfiguration.getPeptideFilter());
	}

	public static DBIndexInterface getFastaDBIndex(File fastaFile, char[] enzymeArray, int missedCleavages,
			boolean semicleavage, PeptideFilter peptideFilter) {
		if (fastaFile != null) {

			DBIndexSearchParams defaultDBIndexParams = DBIndexInterface.getDefaultDBIndexParams(fastaFile);
			String fastaIndexKey = IndexUtil.createFullIndexFileName(defaultDBIndexParams);

			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeArr(enzymeArray, missedCleavages, semicleavage);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeOffset(0);
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setEnzymeNocutResidues("");
			((DBIndexSearchParamsImpl) defaultDBIndexParams).setH2OPlusProtonAdded(true);

			if (peptideFilter != null) {
				((DBIndexSearchParamsImpl) defaultDBIndexParams).setPeptideFilter(peptideFilter);
			}
			DBIndexInterface dbIndex = new DBIndexInterface(defaultDBIndexParams);

			return dbIndex;
		}
		return null;
	}

	public boolean isIgnorePeptidesNotFoundInDB() {
		return ignorePeptidesNotFoundInDB;
	}

}
