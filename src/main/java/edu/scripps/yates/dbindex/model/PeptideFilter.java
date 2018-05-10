package edu.scripps.yates.dbindex.model;

public abstract class PeptideFilter {
	public abstract boolean isValid(String peptideSequence);

	@Override
	public abstract String toString();
}
