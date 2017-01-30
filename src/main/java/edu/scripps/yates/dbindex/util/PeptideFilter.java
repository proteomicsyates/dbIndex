package edu.scripps.yates.dbindex.util;

public abstract class PeptideFilter {
	public abstract boolean isValid(String peptideSequence);

	@Override
	public abstract String toString();
}
