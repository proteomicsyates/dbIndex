package edu.scripps.yates.dbindex;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.fasta.dbindex.IndexedProtein;
import edu.scripps.yates.utilities.fasta.dbindex.IndexedSequence;

/**
 * Protein sequence cache
 * 
 * Supports addition of proteins and constant access of protein sequence by id
 * Abstracts out the underlying storage.
 * 
 * Singleton, so we do not waste memory in multi-threaded setup Thread-safe, can
 * be shared between multiple reader threads.
 * 
 * @author Adam
 */
public class ProteinCache {
	private static final Logger log = Logger.getLogger(ProteinCache.class);
	private final ArrayList<String> sequences;
	private final ArrayList<String> defs;
	public static final int EST_SIZE = 200000;
	private static ProteinCache instance = null;

	// /////////////////
	// Salva 15 Apr 2015
	// Disabled because it cause problems to use the same protein cache for
	// multiple indexes
	// /////////////////
	// public synchronized static ProteinCache getInstance() {
	// if (instance == null) {
	// instance = new ProteinCache();
	// }
	// return instance;
	// }
	// /////////////////

	/**
	 * Create a new empty cache
	 */
	public ProteinCache() {
		sequences = new ArrayList<String>(EST_SIZE);
		defs = new ArrayList<String>(EST_SIZE);
	}

	public int getNumberProteins() {
		return sequences.size();
	}

	/**
	 * get protein by id, protein id start at 1
	 * 
	 * @param proteinId protein id, 1 based
	 * @return protein sequence string
	 */
	public String getProteinSequence(long proteinId) {
		if (sequences.size() > proteinId) {
			return sequences.get(Long.valueOf(proteinId).intValue());
		}
		return null;
	}

	/**
	 * get protein by id, protein id start at 1
	 * 
	 * @param proteinId protein id, 1 based
	 * @return protein sequence string
	 */
	String getProteinDef(long proteinId) {
		return defs.get(Long.valueOf(proteinId).intValue());
	}

	/**
	 * add protein to the cache Assumes the protein is added in the protein
	 * id-ordering
	 * 
	 * @param def     prot. def to add
	 * @param protein sequence of protein to add to the cache
	 */
	public int addProtein(String def, String protein) {
		// make sure that there is no TAB in def, otherwise it will throw an
		// error trying to get the 3rd element separated by TABs
		if (def.contains("\t")) {
			def = def.replace("\t", " ");
		}
		defs.add(def);
		if (protein != null) {
			sequences.add(protein);
		}
		return defs.size() - 1;
	}

	/*
	 * @harshil Shah
	 */
	public int addProtein(String def) {
		return addProtein(def, null);
	}

	/**
	 * Get peptide sequence as protein subsequence, does not check bounds
	 * 
	 * @param protId
	 * @param seqOffset
	 * @param seqLen
	 * @return
	 */
	public String getPeptideSequence(int protId, int seqOffset, int seqLen) {
		final String protSeq = getProteinSequence(protId);
		String peptideSeq = null;

		try {

			peptideSeq = protSeq.substring(seqOffset, seqOffset + seqLen);
		} catch (final Exception e) {
			e.printStackTrace();
			log.error("Error tryin to get substring from " + protSeq);
			log.error("Using beginIndex:" + seqOffset + " and length: " + seqLen);
			log.error(e.getMessage());
		}

		return peptideSeq;
	}

	/**
	 * clear the protein cache
	 */
	void clear() {
		sequences.clear();
		defs.clear();
	}

	/**
	 * check if the cache has been populated
	 * 
	 * @return true if populated / not empty
	 */
	protected boolean isPopulated() {
		return !sequences.isEmpty();
	}

	public List<IndexedProtein> getProteins(IndexedSequence sequence) {
		final List<IndexedProtein> ret = new ArrayList<IndexedProtein>();

		for (final Integer proteinId : sequence.getProteinIds()) {
			final IndexedProtein ip = new IndexedProtein(getProteinDef(proteinId), proteinId);
			ret.add(ip);
		}
		return ret;
	}

	public void clearSequences() {
		sequences.clear();
	}

	public void clearDefs() {
		defs.clear();
	}

	public void addSequence(String seq) {
		sequences.add(seq);
	}

	public int getIndexOfDef(String def) {
		return defs.indexOf(def);
	}

	public String getDef(Integer index) {
		return defs.get(index);
	}

	public String getSequence(Integer index) {
		return sequences.get(index);
	}
}
