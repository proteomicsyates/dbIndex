package edu.scripps.yates.dbindex.util;

/**
 * A {@link PeptideFilter} implementation that discard peptides that contains a
 * AA more than X times
 *
 * @author Salva
 *
 */
public class PeptideFilterByMaxOccurrencies extends PeptideFilter {
	private final char aa;
	private final int numMaxOccurrencies;

	public PeptideFilterByMaxOccurrencies(String aa, int numMax) {
		this.aa = aa.toCharArray()[0];
		numMaxOccurrencies = numMax;
	}

	@Override
	public boolean isValid(String peptideSequence) {
		int numMatches = 0;
		for (Character c : peptideSequence.toCharArray()) {
			if (c == aa) {
				numMatches++;
				if (numMatches > numMaxOccurrencies) {
					return false;
				}
			}
		}

		return true;
	}

	@Override
	public String toString() {
		return String.valueOf(aa) + numMaxOccurrencies;
	}
}
