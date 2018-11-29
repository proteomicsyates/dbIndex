package edu.scripps.yates.dbindex.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.compomics.util.protein.AASequenceImpl;
import com.compomics.util.protein.Enzyme;
import com.compomics.util.protein.Protein;

import edu.scripps.yates.utilities.fasta.dbindex.PeptideFilter;

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
		if (peptideFilterString != null && !"".equals(peptideFilterString)) {
			final String aa = String.valueOf(peptideFilterString.charAt(0));
			final int numMax = Integer.valueOf(String.valueOf(peptideFilterString.charAt(1)));
			return new PeptideFilterByMaxOccurrencies(aa, numMax);
		}
		return null;
	}

	public List<String> digestProtein(String proteinSequence) {
		final Enzyme enzyme = new Enzyme("enzyme name", getEnzymeArrayString(), "", "", numMisscleavages);

		final Protein protein = new Protein(new AASequenceImpl(proteinSequence));
		final Protein[] digest = enzyme.cleave(protein, 6, Integer.MAX_VALUE);
		final List<String> ret = Arrays.asList(digest).stream().map(p -> p.getSequence().getSequence())
				.collect(Collectors.toList());
		final PeptideFilter peptideFilter = getPeptideFilter();
		if (peptideFilter == null) {
			return ret;
		}
		final List<String> ret2 = new ArrayList<String>();
		for (final String peptideSequence : ret) {
			if (peptideFilter.isValid(peptideSequence)) {
				ret2.add(peptideSequence);
			}
		}
		return ret2;
	}

	private String getEnzymeArrayString() {
		final StringBuilder sb = new StringBuilder();
		for (final char c : enzymeArray) {
			sb.append(c);
		}
		return sb.toString();
	}
}
