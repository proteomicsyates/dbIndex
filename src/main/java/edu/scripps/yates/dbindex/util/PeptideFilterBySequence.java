package edu.scripps.yates.dbindex.util;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A {@link PeptideFilter} implementation that discard peptides that match a
 * regular expresion
 *
 * @author Salva
 *
 */
public class PeptideFilterBySequence extends PeptideFilter {
	private final Pattern pattern;

	public PeptideFilterBySequence(String regexp) {
		pattern = Pattern.compile(regexp);
	}

	@Override
	public boolean isValid(String peptideSequence) {
		final Matcher matcher = pattern.matcher(peptideSequence);
		if (matcher.find()) {
			return false;
		}
		return true;
	}

	@Override
	public String toString() {
		return pattern.pattern();
	}
}
