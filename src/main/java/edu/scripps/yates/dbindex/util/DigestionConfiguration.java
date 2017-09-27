package edu.scripps.yates.dbindex.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DigestionConfiguration {
	private final char[] enzymeArray;
	private final int numMisscleavages;
	private final boolean semiCleavage;
	private final String peptideFilterString;

	public DigestionConfiguration(char[] enzymeArray, int numMisscleavages, boolean semiCleavage,
			String peptideFilterString) {

		this.enzymeArray = enzymeArray;
		this.numMisscleavages = numMisscleavages;
		this.semiCleavage = semiCleavage;
		this.peptideFilterString = peptideFilterString;
	}

	/**
	 * @return the enzymeArray
	 */
	public char[] getEnzymeArray() {
		return enzymeArray;
	}

	/**
	 * @return the numMisscleavages
	 */
	public int getNumMisscleavages() {
		return numMisscleavages;
	}

	/**
	 * @return the semiCleavage
	 */
	public boolean isSemiCleavage() {
		return semiCleavage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "[enzymeArray=" + Arrays.toString(enzymeArray) + ", numMisscleavages=" + numMisscleavages
				+ ", semiCleavage=" + semiCleavage + ", peptideFilterString=" + peptideFilterString + "]";
	}

	public PeptideFilter getPeptideFilter() {
		if (peptideFilterString != null) {
			String aa = String.valueOf(peptideFilterString.charAt(0));
			int numMax = Integer.valueOf(String.valueOf(peptideFilterString.charAt(1)));
			return new PeptideFilterByMaxOccurrencies(aa, numMax);
		}
		return null;
	}

	public List<String> digestProtein(String proteinSequence) {
		com.compomics.util.experiment.biology.Enzyme enzyme = new com.compomics.util.experiment.biology.Enzyme(1,
				"enzyme name", getEnzymeArrayString(), "", "", "");
		ArrayList<String> digest = enzyme.digest(proteinSequence, numMisscleavages, 6, Integer.MAX_VALUE);
		PeptideFilter peptideFilter = getPeptideFilter();
		if (peptideFilter == null) {
			return digest;
		}
		List<String> ret = new ArrayList<String>();
		for (String peptideSequence : digest) {
			if (peptideFilter.isValid(peptideSequence)) {
				ret.add(peptideSequence);
			}
		}
		return ret;
	}

	private String getEnzymeArrayString() {
		StringBuilder sb = new StringBuilder();
		for (char c : enzymeArray) {
			sb.append(c);
		}
		return sb.toString();
	}
}
