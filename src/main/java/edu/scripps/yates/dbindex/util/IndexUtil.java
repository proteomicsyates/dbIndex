package edu.scripps.yates.dbindex.util;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.RandomAccessFile;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.annotations.uniprot.fasta.FastaReaderFromUniprot;
import edu.scripps.yates.dbindex.Constants;
import edu.scripps.yates.dbindex.SearchParams;
import edu.scripps.yates.dbindex.Util;
import edu.scripps.yates.utilities.fasta.FastaReader;
import edu.scripps.yates.utilities.fasta.dbindex.DBIndexSearchParams;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;
import edu.scripps.yates.utilities.fasta.dbindex.PeptideFilter;
import edu.scripps.yates.utilities.fasta.dbindex.ResidueInfo;
import edu.scripps.yates.utilities.masses.AssignMass;

/**
 * Indexing util
 * 
 * @author Adam
 */
public class IndexUtil {

	private static final Logger log = Logger.getLogger(IndexUtil.class);

	public static String getMd5(String in) {
		try {
			final MessageDigest md5 = MessageDigest.getInstance("MD5");
			md5.update(in.getBytes(), 0, in.length());

			return new BigInteger(1, md5.digest()).toString(16);

		} catch (final NoSuchAlgorithmException e) {
			log.log(Level.INFO, "Could not calculate md5sum ", e);
			return "";
		}
	}

	public static byte[] getMd5Bytes(String in) {
		try {
			final MessageDigest md5 = MessageDigest.getInstance("MD5");
			md5.update(in.getBytes(), 0, in.length());
			return md5.digest();

		} catch (final NoSuchAlgorithmException e) {
			log.log(Level.INFO, "Could not calculate md5sum ", e);
			return new byte[1];
		}
	}

	/**
	 * get last N lines from a file (fast)
	 * 
	 * @param file
	 * @param lines
	 * @return
	 */
	public static String tail(File file, int lines) {
		RandomAccessFile fileHandler = null;
		try {
			fileHandler = new java.io.RandomAccessFile(file, "r");
			final long fileLength = file.length() - 1;
			final StringBuilder sb = new StringBuilder();
			int line = 0;

			for (long filePointer = fileLength; filePointer != -1; filePointer--) {
				fileHandler.seek(filePointer);
				final int readByte = fileHandler.readByte();

				if (readByte == 0xA) {
					if (line == lines) {
						if (filePointer == fileLength) {
							continue;
						} else {
							break;
						}
					}
				} else if (readByte == 0xD) {
					line = line + 1;
					if (line == lines) {
						if (filePointer == fileLength - 1) {
							continue;
						} else {
							break;
						}
					}
				}
				sb.append((char) readByte);
			}

			sb.deleteCharAt(sb.length() - 1);
			final String lastLine = sb.reverse().toString();
			return lastLine;
		} catch (final java.io.FileNotFoundException e) {
			e.printStackTrace();
			return null;
		} catch (final java.io.IOException e) {
			e.printStackTrace();
			return null;
		} finally {
			if (fileHandler != null) {
				try {
					fileHandler.close();
				} catch (final IOException ex) {
					Logger.getLogger(IndexUtil.class).log(Level.ERROR, null, ex);
				}
			}
		}
	}

	/**
	 * Get number of lines in a file (fast)
	 * 
	 * @param filename
	 * @return
	 * @throws IOException
	 */
	public static int countLines(String filename) throws IOException {
		final LineNumberReader reader = new LineNumberReader(new FileReader(filename));
		int cnt = 0;
		String lineRead = "";
		while ((lineRead = reader.readLine()) != null) {
		}

		cnt = reader.getLineNumber();
		reader.close();
		return cnt;
	}

	public static ResidueInfo getResidues(IndexedSequence peptideSequence, int seqOffset, int seqLen,
			String proteinSequence) {
		final int protLen = proteinSequence.length();

		final int resLeftI = seqOffset >= Constants.MAX_INDEX_RESIDUE_LEN ? seqOffset - Constants.MAX_INDEX_RESIDUE_LEN
				: 0;
		final int resLeftLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, seqOffset);
		final StringBuilder sbLeft = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
		sbLeft.append(proteinSequence.substring(resLeftI, resLeftI + resLeftLen));

		final int end = seqOffset + seqLen;
		final int resRightI = end; // + 1;
		final int resRightLen = Math.min(Constants.MAX_INDEX_RESIDUE_LEN, protLen - end - 1);
		final StringBuilder sbRight = new StringBuilder(Constants.MAX_INDEX_RESIDUE_LEN);
		if (resRightI < protLen) {
			sbRight.append(proteinSequence.substring(resRightI, resRightI + resRightLen));
		}

		// add -- markers to fill Constants.MAX_INDEX_RESIDUE_LEN length
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

		return new ResidueInfo(resLeft, resRight);
	}

	/**
	 * Calculates the mass of a sequence, by summing all the masses of all the AAs
	 * and then summing the mass of a water molecule H2O and the proton mass, as
	 * well as the mass of the cTerm and nTerm if configured in {@link AssignMass}
	 *
	 * @param sequence
	 * @param isH2OPlusProtonAdded
	 * @return
	 */
	public static double calculateMass(String sequence) {
		return calculateMass(sequence, true);
	}

	/**
	 * Calculates the mass of a sequence, by summing all the masses of all the AAs
	 * and then summing the mass of a water molecule H2O and the proton mass (is the
	 * parameter is true), as well as the mass of the cTerm and nTerm if configured
	 * in {@link AssignMass}
	 *
	 * @param sequence
	 * @param isH2OPlusProtonAdded
	 * @return
	 */
	public static double calculateMass(String sequence, boolean isH2OPlusProtonAdded) {
		double mass = 0;
		if (isH2OPlusProtonAdded)
			mass += AssignMass.H2O_PROTON;
		mass += AssignMass.getcTerm();
		mass += AssignMass.getnTerm();
		for (int i = 0; i < sequence.length(); i++) {
			// System.out.println(AssignMass.getMass(sequence.charAt(i)));
			mass += AssignMass.getMass(sequence.charAt(i));
		}
		return mass;
	}

	/**
	 * This is the number of rows to lookup in the index from the key it is
	 * calculated from the mass to search.<br>
	 * So, if someone looks to a mass of 1340.386, it will be multiplied by the
	 * massGroupFactor, and then, then integer part will be taken as the key (if
	 * massGroupFactor is 1, then 1340 is the key). Then, the index will look into
	 * 1340-1, 1340 and 1340+1 entries in the SQLLite database, being the 1, in this
	 * case, the value of numRowsToLookup.<br>
	 * Remain this on 0, since from the last change by Salva 24Nov2014, setting
	 * massGroupFactor to 10000 this is not needed.
	 */
	private static int numRowsToLookup = 0;// 1;

	public static void setNumRowsToLookup(int numRowsToLookup) {
		IndexUtil.numRowsToLookup = numRowsToLookup;
	}

	public static int getNumRowsToLookup() {
		return numRowsToLookup;
	}

	/**
	 * Gets the tolerance in Da referred to an actualMass and the ppmTolerance
	 *
	 * @param actualMass
	 * @param ppmTolerance
	 * @return
	 */
	public static double getToleranceInDalton(double actualMass, double ppmTolerance) {
		return actualMass * (1 - 1 / (ppmTolerance / Constants.ONE_MILLION + 1));
	}

	public static String createFullIndexFileName(DBIndexSearchParams params) {
		return createFullIndexFileName(params, null, null, false, null, false, null);
	}

	public static String createFullIndexFileName(DBIndexSearchParams params, String sufix,
			Integer maxVariationsPerPeptide, boolean useUniprot, String uniprotVersion, boolean usePhophoSite,
			String phosphoSiteSpecies) {
		return createFullIndexFileName(params.getDatabaseName(), sufix, params.getEnzymeOffset(),
				params.getEnzymeResidues(), params.getEnzymeNocutResidues(), params.getMaxMissedCleavages(),
				params.getMinPrecursorMass(), params.getMaxPrecursorMass(), SearchParams.getStaticParams(),
				params.isSemiCleavage(), params.getPeptideFilter(), params.isH2OPlusProtonAdded(),
				params.getMassGroupFactor(), params.getMandatoryInternalAAs(), params.isLookProteoforms(),
				maxVariationsPerPeptide, useUniprot, uniprotVersion, usePhophoSite, phosphoSiteSpecies,
				params.getDiscardDecoyRegexp(), params.isUseMonoParent());
	}

	public static String createFullIndexFileName(SearchParams params, String sufix, Integer maxVariationsPerPeptide,
			boolean useUniprot, String uniprotVersion, boolean usePhophoSite, String phosphoSiteSpecies) {
		return createFullIndexFileName(params.getDatabaseName(), sufix, params.getEnzymeOffset(),
				params.getEnzymeResidues(), params.getEnzymeNocutResidues(), params.getMaxMissedCleavages(),
				params.getMinPrecursorMass(), params.getMaxPrecursorMass(), SearchParams.getStaticParams(),
				params.isSemiCleavage(), params.getPeptideFilter(), params.isH2OPlusProtonAdded(),
				params.getMassGroupFactor(), params.getMandatoryInternalAAs(), params.isLookProteoforms(),
				maxVariationsPerPeptide, useUniprot, uniprotVersion, usePhophoSite, phosphoSiteSpecies,
				params.getDiscardDecoyRegexp(), params.isUseMonoParent());

	}

	public static String createFullIndexFileName(String databaseName, String sufix, int enzymeOffset,
			String enzymeResidues, String enzymeNocutResidues, int maxMissedCleavages, double minPrecursorMass,
			double maxPrecursorMass, StringBuffer staticParams, boolean semiCleavage, PeptideFilter peptideFilter,
			boolean H2OPlusProtonAdded, int massGroupFactor, char[] mandatoryInternalAAs, Boolean isLookProteoforms,
			Integer maxVariationsPerPeptide, boolean useUniprot, String uniprotVersion, boolean usePhophoSite,
			String phosphoSiteSpecies, String discardDecoys, boolean mono) {
		final String uniqueIndexName = databaseName + "_";

		// generate a unique string based on current params that affect the
		// index
		final StringBuilder uniqueParams = new StringBuilder();
		// uniqueParams.append(sparam.getEnzyme().toString());
		// uniqueParams.append(sparam.getEnzymeNumber());
		uniqueParams.append("enzymeOffset:" + enzymeOffset);
		uniqueParams.append(", enzymeResidues:").append(enzymeResidues);
		uniqueParams.append(", enzymeNoCutResidues:").append(enzymeNocutResidues);

		// uniqueParams.append(" ").append(getIndexFactor() ) ;

		uniqueParams.append(", maxCleavages:").append(maxMissedCleavages);
		uniqueParams.append(", minPrecursorMass:" + minPrecursorMass);
		uniqueParams.append(", maxPrecursorMass:" + maxPrecursorMass);
		uniqueParams.append(", static:").append(staticParams);
		uniqueParams.append(", semiCleave:" + semiCleavage);

		if (peptideFilter != null) {
			uniqueParams.append(", pepFilter:").append(peptideFilter.toString());
		}

		uniqueParams.append(", isH2OPlusProtonAdded:" + H2OPlusProtonAdded);

		uniqueParams.append(", massGroupFactor:" + massGroupFactor);
		uniqueParams.append(", massType:" + mono);
		if (mandatoryInternalAAs != null) {
			uniqueParams.append(", mandatoryInternalAAs:" + mandatoryInternalAAs.toString());
		}
		// added by Salva 29March 2018

		uniqueParams.append(", proteoForms:" + (isLookProteoforms != null && isLookProteoforms));

		uniqueParams.append(", maxVariationsPerPeptide:" + maxVariationsPerPeptide);

		uniqueParams.append(", useUniprot:" + useUniprot);
		uniqueParams.append(", uniprotVersion:" + (useUniprot ? uniprotVersion : null));
		uniqueParams.append(", usePhosphoSite:" + usePhophoSite);
		uniqueParams.append(", phosphoSiteSpecies:" + (usePhophoSite ? phosphoSiteSpecies : null));
		uniqueParams.append(", sufix:" + sufix);
		uniqueParams.append(", discardDecoys:" + discardDecoys);

		final String uniqueParamsStr = uniqueParams.toString();
		final String uniqueParamsStrHash = Util.getMd5(uniqueParamsStr);
		log.info("Index Unique String:" + uniqueParamsStr);
		log.info("Index Unique hashkey:" + uniqueParamsStrHash);
		return uniqueIndexName + uniqueParamsStrHash;
	}

	public static FastaReader getFastaReader(DBIndexSearchParams params) {
		if (params.isLookProteoforms() != null && params.isLookProteoforms()) {

//			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(
//					params.getUniprotReleasesFolder(), true);
//			final UniprotProteoformRetriever proteoFormRetriever = new UniprotProteoformRetrieverFromXML(uplr,
//					params.getUniprotVersion());
//			return new ProteoFormFastaReader(params.getDatabaseName(), proteoFormRetriever);

			final File uniprotReleasesFolder = params.getUniprotReleasesFolder();
			final boolean useIndex = params.isUseIndex();
			final UniprotProteinLocalRetriever uplr = new UniprotProteinLocalRetriever(uniprotReleasesFolder, useIndex);
			uplr.setRetrieveFastaIsoformsFromMainForms(true);
			uplr.setRetrieveFastaIsoforms(true);
			final String fastaFileName = params.getDatabaseName();
			return new FastaReaderFromUniprot(fastaFileName, params.getUniprotVersion(), uplr, true);
		} else {
			return new FastaReader(params.getDatabaseName());
		}
	}
}
