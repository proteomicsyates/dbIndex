package edu.scripps.yates.dbindex;

import edu.scripps.yates.utilities.fasta.FastaParser;

/**
 * A representation of an indexed protein (a fasta header)
 * 
 * @author Adam
 */
public class IndexedProtein {
	private long id;
	private String accession;
	private String fastaDefLine;

	IndexedProtein(long id) {
		this.id = id;
	}

	IndexedProtein(String fastaDefLine, long id) {
		// this.accession = Fasta.getAccession(fastaDefLine);
		// accession = Fasta.getSequestLikeAccession(fastaDefLine);
		accession = FastaParser.getACC(fastaDefLine).getFirstelement();
		setFastaDefLine(fastaDefLine);
		this.id = id;
	}

	public long getId() {
		return id;
	}

	public void setId(long id) {
		this.id = id;
	}

	@Override
	public String toString() {
		return "IndexedProtein{" + "accession=" + accession + ", id=" + id + '}';
	}

	public String getAccession() {
		return accession;
	}

	public void setAccession(String accession) {
		this.accession = accession;
	}

	/**
	 * @return the fastaDefLine
	 */
	public String getFastaDefLine() {
		return fastaDefLine;
	}

	/**
	 * @param fastaDefLine
	 *            the fastaDefLine to set
	 */
	public void setFastaDefLine(String fastaDefLine) {
		this.fastaDefLine = fastaDefLine;
	}

}
