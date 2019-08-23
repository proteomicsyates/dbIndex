/*
 * To change this template, choose Tools | Templates and open the template in
 * the editor.
 */
package edu.scripps.yates.dbindex.io;

/**
 *
 * @author rpark
 */
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.io.FilenameUtils;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import edu.scripps.yates.dbindex.SearchParams;
import edu.scripps.yates.dbindex.model.AssignMassToStaticParam;
import edu.scripps.yates.dbindex.model.DiffModification;
import edu.scripps.yates.dbindex.model.ModResidue;
import edu.scripps.yates.utilities.fasta.dbindex.IndexType;
import edu.scripps.yates.utilities.masses.AssignMass;
import gnu.trove.set.hash.THashSet;

/**
 * @author Robin Park
 * @version $Id: SearchParamReader.java,v 1.1 2014/04/25 22:33:09 linhe Exp $
 */
public class SearchParamReader {

	private String path;
	private String fileName;
	private SearchParams param;
	// private boolean isModification;
	private Hashtable<String, String> ht = new Hashtable<String, String>();
	private final Logger logger = Logger.getLogger(SearchParamReader.class.getName());
	private final char[] MOD_SYMBOL = { '*', '#', '@', };
	// private TIntDoubleHashMap symbolMap = new TIntDoubleHashMap(10);
	private final double[] massShiftArr = new double[3];
	public static final String DEFAULT_PARAM_FILE_NAME = "blazmass.params";

	public static void main(String args[]) throws Exception {

		if (args.length < 2) {
			System.out.println("Usage: java SearchParamReader path param_filename");
			return;
		}

		final SearchParamReader p = new SearchParamReader(args[0], args[1]);
		// else
		// p = new SearchParamReader(".", DEFAULT_PARAM_FILE_NAME);

		final SearchParams param = p.getParam();

		// System.out.println("===" + param.getFullIndexFileName());
		// System.out.println("dbname===============" +
		// param.getFullIndexFileName());
		// System.out.println("dbname===============" +
		// param.getIndexFileNameOnly());

		/*
		 * ParamReader p = new ParamReader(args[0], args[1], Integer.parseInt(args[2]));
		 * Hashtable<String, String> h = p.getHashtable();
		 * System.out.println(h.get("diff_search_options")); System.out.println(
		 * p.isModification() ); List l = p.getResidueList(); for(int
		 * i=0;i<l.size();i++) { ModResidue residue = (ModResidue)l.get(i);
		 * System.out.println(residue.getMassShift() + "\t" + residue.getResidue()); }
		 * double[] d = p.getMassShiftArr(); for(int i=0;i<d.length;i++) {
		 * System.out.println(d[i]); }
		 */
	}

	public SearchParamReader(String path, String fileName) throws IOException {
		this.path = path;
		this.fileName = fileName;
		init();
	}

	public SearchParamReader(File file) throws IOException {
		fileName = FilenameUtils.getName(file.getAbsolutePath());
		path = FilenameUtils.getFullPathNoEndSeparator(file.getAbsolutePath());
		init();
	}

	public void init() throws IOException {
		path = path + File.separator + fileName;
		if (!new File(path).exists()) {
			throw new IOException("Could not locate params file at path: " + path);
		}
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(path));
			String eachLine;

			// Hashtable<String> ht = new Hashtable<String>();

			param = SearchParams.getInstance();
			// StringBuffer sb = new StringBuffer();
			// br.skip(2000);
			while ((eachLine = br.readLine()) != null) {
				// sb.append(eachLine);
				// sb.append("\n");

				// System.out.println("1----------" + eachLine);

				if (eachLine.startsWith("#")) // ||
				// (strArr = eachLine.split("=")).length < 2)
				{
					continue;
				}

				// System.out.println("-----------" + eachLine);
				/*
				 * if (eachLine.startsWith("[BLAZMASS_ENZYME_INFO]") ||
				 * eachLine.startsWith("[SEQUEST_ENZYME_INFO]")) { int enzymeNumber =
				 * trimValueAsInt(getParam("enzyme_number")); int currentEnz; while ((eachLine =
				 * br.readLine()) != null) { String[] arr = eachLine.split("\\."); if
				 * (arr.length < 2) { continue; } currentEnz = Integer.parseInt(arr[0]);
				 * //dTempMass: where is it used? if (currentEnz == enzymeNumber) { String[]
				 * tarr = eachLine.split("\\s+"); //szEnzymeName = tarr[1].trim();
				 * param.setEnzymeOffset(trimValueAsInt(tarr[2]));
				 * param.setEnzymeBreakAA(tarr[3].trim());
				 * param.setEnzymeNoBreakAA(tarr[4].trim()); } } if (enzymeNumber == 0) {
				 * param.setUseEnzyme(false); } else { param.setUseEnzyme(true); } }
				 */

				if (null == eachLine) {
					break;
				}
				final String[] strArr = eachLine.split("=");
				if (strArr.length < 2) {
					continue;
				}

				ht.put(strArr[0].trim(), strArr[1].split(";")[0].trim());

			}

			param.setDatabaseName(trimValue(getParam("database_name")));
			param.setIndexDatabaseName(trimValue(getParam("index_database_name")));

			final int useIndex = trimValueAsInt(getParam("use_index"), 1);
			param.setUseIndex(useIndex == 1);

			final int indexType = trimValueAsInt(getParam("index_type"), 1);
			param.setIndexType(indexType == 1 ? IndexType.INDEX_NORMAL : IndexType.INDEX_LARGE);

			final int inMemoryIndexI = trimValueAsInt(getParam("index_inmemory"), 1);
			final boolean inMemoryIndex = inMemoryIndexI == 1;
			param.setInMemoryIndex(inMemoryIndex);

			final int indexFactor = trimValueAsInt(getParam("index_factor"), 6);
			param.setIndexFactor(indexFactor);

			param.setParametersFile(path);
			// param.setParameters(sb.toString());

			// String pepTolerance = getParam("peptide_mass_tolerance");
			// if (null == pepTolerance) {
			final String pepTolerance = getParam("ppm_peptide_mass_tolerance");
			if (null == pepTolerance)
				System.out.println("Error: ppm_peptide_mass_tolerance is missing");
			// }

			param.setPeptideMassTolerance(trimValueAsDouble(pepTolerance));

			param.setFragmentIonTolerance(trimValueAsDouble(getParam("ppm_fragment_ion_tolerance")));
			param.setFragmentIonToleranceInt((int) trimValueAsDouble((getParam("ppm_fragment_ion_tolerance"))));

			final String highRes = getParam("ppm_fragment_ion_tolerance_high");
			if (null != highRes && "1".equals(highRes))
				param.setHighResolution(true);
			else
				param.setHighResolution(false);

			// String mongodbServer = getParam("mongodb_server");
			//
			// if (null != mongodbServer)
			// param.setMongodbServer(mongodbServer);

			final String mongoParam = getParam("use_mongodb");
			final String seqDBParam = getParam("use_SeqDB");
			final String protDBParam = getParam("use_ProtDB");
			if (mongoParam == null || Integer.parseInt(mongoParam) == 0) {
				param.setUsingMongoDB(false);
//				System.out.println("Must use mongodb");
//				System.exit(1);
			} else {
				// force dbSize
				param.setUsingMongoDB(true);
				final String mongoDBURI = getParam("mongoDB_URI");
				if (mongoDBURI != null) {
					param.setMongoDBURI(mongoDBURI);
				} else {
					final String massDBURI = getParam("MassDB_URI");
					if (massDBURI != null) {
						param.setMongoDBURI(massDBURI);
					} else {
						final String seqDBURI = getParam("SeqDB_URI");
						if (seqDBURI != null) {
							param.setMongoDBURI(seqDBURI);
						} else {
							final String protDBURI = getParam("ProtDB_URI");
							if (protDBURI != null) {
								param.setMongoDBURI(protDBURI);
							}
						}
					}
				}
				param.setMassDBName(getParam("MassDB_dbname"));
				param.setMassDBCollection(getParam("MassDB_collection"));
				param.setDatabaseName(param.getMassDBName() + "." + param.getMassDBCollection() + " [MongoDB]");

				if (seqDBParam != null && Integer.parseInt(seqDBParam) == 1) {
					// use SeqDB
					param.setUsingSeqDB(true);
					final String seqDBName = getParam("SeqDB_dbname");
					final String seqDBCollection = getParam("SeqDB_collection");
					param.setSeqDBName(seqDBName);
					param.setSeqDBCollection(seqDBCollection);

					if (protDBParam != null && Integer.parseInt(protDBParam) == 1) {
						// use ProtDB
						param.setUsingProtDB(true);
						param.setProtDBName(getParam("ProtDB_dbname"));
						param.setProtDBCollection(getParam("ProtDB_collection"));
					} else {
						param.setUsingProtDB(false);
						System.out.println("Not using ProtDB MongoDB database -- SQT output will be incomplete");
					}
				} else {
					param.setUsingSeqDB(false);
					System.out.println("Not using SeqDB MongoDB database -- SQT output will be incomplete");
				}
			}
			final String massTypeParent = getParam("mass_type_parent");

			if ("0".equals(massTypeParent)) {
				param.setUseMonoParent(false);
			} else if ("1".equals(massTypeParent)) {
				param.setUseMonoParent(true);
			} else {
				param.setUseMonoParent(true);
			}
			param.setMassTypeParent(trimValueAsInt(massTypeParent)); // ;
																		// 0=average
																		// masses,
																		// 1=monoisotopic
																		// masses

			final String massTypeFrag = getParam("mass_type_fragment");

			// System.out.println("=====---=====" + massTypeParent);
			// System.out.println("=====---=====" + massTypeFrag);
			if ("1".equals(massTypeFrag)) {
				param.setUseMonoFragment(true);
			} else {
				param.setUseMonoFragment(false);
			}
			param.setMassTypeFragment(trimValueAsInt(massTypeFrag)); // ;
																		// 0=average
																		// masses,
																		// 1=monoisotopic
																		// masses

			param.setNumPeptideOutputLnes(trimValueAsInt(getParam("num_output_lines")));
			param.setRemovePrecursorPeak(trimValueAsInt(getParam("remove_precursor_peak")));
			final String ionSeries = getParam("ion_series");
			if (ionSeries != null) {
				final String arr[] = ionSeries.split(" ");

				// System.out.println("===========" + arr[0] + arr[1] + " " +
				// arr[2]);
				param.setNeutralLossAions(Integer.parseInt(arr[0]));
				param.setNeutralLossBions(Integer.parseInt(arr[1]));
				param.setNeutralLossYions(Integer.parseInt(arr[2]));

				final int[] ions = new int[9];
				final double[] weightArr = new double[12];

				for (int i = 0; i < 9; i++)
					weightArr[i] = Double.parseDouble(arr[i + 3]);

				for (int i = 0; i < 3; i++)
					weightArr[i + 9] = Double.parseDouble(arr[i]);

				param.setWeightArr(weightArr);

				ions[0] = (int) (10 * Double.parseDouble(arr[3])); // a
				ions[1] = (int) (10 * Double.parseDouble(arr[4])); // b
				ions[2] = (int) (10 * Double.parseDouble(arr[5])); // c
				ions[3] = (int) (10 * Double.parseDouble(arr[6]));
				ions[4] = (int) (10 * Double.parseDouble(arr[7]));
				ions[5] = (int) (10 * Double.parseDouble(arr[8]));
				ions[6] = (int) (10 * Double.parseDouble(arr[9])); // x
				ions[7] = (int) (10 * Double.parseDouble(arr[10])); // y
				ions[8] = (int) (10 * Double.parseDouble(arr[11])); // z

				int numIonUsed = 0;
				for (final int eachIon : ions) {
					if (eachIon > 0)
						numIonUsed++;
				}

				param.setNumIonSeriesUsed(numIonUsed);
				param.setIonSeries(ions);
			}

			param.setMaxNumDiffMod(trimValueAsInt(getParam("max_num_differential_AA_per_mod")));

			final Object obj = getParam("peptide_mass_tolerance");
			if (null != obj) {
				param.setPeptideMassTolerance(trimValueAsDouble(obj.toString()));
			}

			final String varMassTol = getParam("var_peptide_mass_tolerance");
			if (null != varMassTol) {
				param.setVariableTolerance(true);
				param.setVariablePeptideMassTolerance(trimValueAsDouble(varMassTol));
			}

			final String amuMassTol = getParam("amu_peptide_mass_tolerance");
			if (null != amuMassTol) {
				param.setPeptideMassTolerance(trimValueAsDouble(amuMassTol));
			}

			final String ppmMassTol = getParam("ppm_peptide_mass_tolerance");
			if (null != ppmMassTol) {
				param.setUsePPM(true);
				param.setRelativePeptideMassTolerance(trimValueAsDouble(ppmMassTol) / 1000f);
			}

			param.setIsotopes(trimValueAsInt(getParam("isotopes")));
			final String n15enrich = getParam("n15_enrichment");
			if (null != n15enrich) {
				param.setN15Enrichment(trimValueAsDouble(n15enrich));
			}

			param.setMatchPeakTolerance(trimValueAsDouble(getParam("match_peak_tolerance")));
			param.setNumPeptideOutputLnes(trimValueAsInt(getParam("num_output_lines")));
			param.setRemovePrecursorPeak(trimValueAsInt(getParam("remove_precursor_peak")));
			// param.setAddCterminus(trimValueAsDouble(getParam("add_C_terminus")));
			// param.setAddNterminus(trimValueAsDouble(getParam("add_N_terminus")));
			AssignMass.setcTerm(trimValueAsDouble(getParam("add_C_terminus")));
			if (AssignMass.getcTerm() > 0)
				SearchParams.addStaticParam("cterm", AssignMass.getcTerm());

			AssignMass.setnTerm(trimValueAsDouble(getParam("add_N_terminus")));
			if (AssignMass.getnTerm() > 0)
				SearchParams.addStaticParam("nterm", AssignMass.getnTerm());

			final AssignMass amassPar = new AssignMass(param.isUseMonoParent());
			param.setHplusparent(AssignMass.getHplus());
			param.setHparent(AssignMass.getH());
			param.setoHparent(AssignMass.getOh());

			final String minPrecursor = getParam("min_precursor_mass");
			if (null != minPrecursor)
				param.setMinPrecursorMass(Double.parseDouble(minPrecursor));

			final String maxPrecursor = getParam("max_precursor_mass");
			if (null != maxPrecursor)
				param.setMaxPrecursorMass(Double.parseDouble(maxPrecursor));
			final String minFragPeakNum = getParam("min_frag_num");
			if (null != minFragPeakNum)
				param.setMinFragPeakNum(Integer.parseInt(minFragPeakNum));

			final AssignMass amassFrag = new AssignMass(param.isUseMonoFragment());

			System.out.println("============" + param.isUseMonoParent() + " " + param.isUseMonoFragment());

			// param.setYionfragment(param.getAddCterminus() + amassPar.getOh()
			// + amassPar.getH() + amassPar.getH());
			// param.setBionfragment(param.getAddNterminus() + amassPar.getH());
			AssignMass.setBionfragment(AssignMass.getnTerm() + AssignMass.getH());
			AssignMass.setYionfragment(
					AssignMass.getcTerm() + AssignMass.getOh() + AssignMass.getH() + AssignMass.getH());

			// System.out.println("============" + (amassPar.getOh() +
			// amassPar.getH() + amassPar.getH()));
			// System.out.println(AssignMass.getcTerm() + "============" +
			// (amassPar.getOh() + "\t" + amassPar.getH() + "\t" +
			// amassPar.getH()));

			param.setBinWidth(AssignMass.getBinWidth());

			// amassPar.addMass('G', 1.1f);
			double f = Double.parseDouble(getParam("add_C_terminus", true));
			AssignMass.setcTerm(f);
			// amassFrag.setcTerm(f);

			f = Double.parseDouble(getParam("add_N_terminus", true));
			AssignMass.setnTerm(f);
			// amassFrag.setnTerm(f);

			f = Double.parseDouble(getParam("add_G_Glycine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('G', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('G', f);

			f = Double.parseDouble(getParam("add_A_Alanine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('A', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('A', f);

			f = Double.parseDouble(getParam("add_S_Serine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('S', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('S', f);

			f = Double.parseDouble(getParam("add_P_Proline", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('P', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('P', f);
			// amassFrag.addMass('P', f );

			f = Double.parseDouble(getParam("add_V_Valine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('V', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('V', f);
			// amassFrag.addMass('V', f );

			f = Double.parseDouble(getParam("add_T_Threonine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('T', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('T', f);
			// amassFrag.addMass('T', f );

			// System.out.println("====" + AssignMass.getMass('C') + " " +
			// param.getN15Enrichment());
			// System.out.println("====>>" + f);
			f = Double.parseDouble(getParam("add_C_Cysteine", true));
			if (param.getN15Enrichment() > 0)
				f = f + param.getN15Enrichment() * (f - 57.02146f);
			// amassPar.addMass('C', f*param.getN15Enrichment());

			AssignMassToStaticParam.addMassAndStaticParam('C', f);

			/*
			 * if (f > 56.9f) { System.out.println("===" + param.getN15Enrichment() + " " +
			 * f); f = f + param.getN15Enrichment() * (f - 57.02146f);
			 * System.out.println("===" + param.getN15Enrichment() + " " + f); }
			 */

			// amassFrag.addMass('C', f );

			f = Double.parseDouble(getParam("add_L_Leucine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('L', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('L', f);
			// amassFrag.addMass('L', f );

			f = Double.parseDouble(getParam("add_I_Isoleucine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('I', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('I', f);
			// amassFrag.addMass('I', f );

			f = Double.parseDouble(getParam("add_X_LorI", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('X', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('X', f);
			// amassFrag.addMass('X', f );

			f = Double.parseDouble(getParam("add_N_Asparagine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('N', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('N', f);
			// amassFrag.addMass('N', f );

			f = Double.parseDouble(getParam("add_O_Ornithine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('O', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('O', f);
			// amassFrag.addMass('O', f );

			f = Double.parseDouble(getParam("add_B_avg_NandD", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('B', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('B', f);
			// amassFrag.addMass('B', f );

			f = Double.parseDouble(getParam("add_D_Aspartic_Acid", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('D', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('D', f);
			// amassFrag.addMass('D', f );

			f = Double.parseDouble(getParam("add_Q_Glutamine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('Q', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('Q', f);
			// amassFrag.addMass('Q', f );

			f = Double.parseDouble(getParam("add_K_Lysine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('K', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('K', f);
			// amassFrag.addMass('K', f );

			f = Double.parseDouble(getParam("add_Z_avg_QandE", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('Z', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('Z', f);
			// amassFrag.addMass('Z', f );

			f = Double.parseDouble(getParam("add_E_Glutamic_Acid", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('E', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('E', f);
			// amassFrag.addMass('E', f );

			f = Double.parseDouble(getParam("add_M_Methionine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('M', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('M', f);
			// amassFrag.addMass('M', f );

			f = Double.parseDouble(getParam("add_H_Histidine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('H', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('H', f);
			// amassFrag.addMass('H', f );

			f = Double.parseDouble(getParam("add_F_Phenyalanine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('F', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('F', f);
			// amassFrag.addMass('F', f );

			f = Double.parseDouble(getParam("add_R_Arginine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('R', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('R', f);
			// amassFrag.addMass('R', f );

			f = Double.parseDouble(getParam("add_Y_Tyrosine", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('Y', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('Y', f);
			// amassFrag.addMass('Y', f );

			f = Double.parseDouble(getParam("add_W_Tryptophan", true));
			if (param.getN15Enrichment() > 0)
				AssignMassToStaticParam.addMassAndStaticParam('W', f * param.getN15Enrichment());
			else
				AssignMassToStaticParam.addMassAndStaticParam('W', f);
			// amassFrag.addMass('W', f );

			param.setPdAAMassParent(amassPar.getPdAAMass());
			param.setPdAAMassFragment(amassFrag.getPdAAMass());

			// param.setMaxInternalCleavageSites(Integer
			// .parseInt(getParam("max_num_internal_cleavage_sites")));
			final String numOutputLinesString = getParam("num_output_lines");
			if (numOutputLinesString != null && !"".equals(numOutputLinesString)) {
				param.setNumPeptideOutputLnes(Integer.parseInt(numOutputLinesString));
			}
			final String maxNumInternalCleavageSites = getParam("max_num_internal_cleavage_sites");
			if (maxNumInternalCleavageSites != null && !"".equals(maxNumInternalCleavageSites)) {
				param.setMaxMissedCleavages(Integer.parseInt(maxNumInternalCleavageSites));
			}
			final boolean parseBoolean = Boolean.parseBoolean(getParam("miscleavage"));
			param.setSemiCleavage(parseBoolean);
			param.setEnzymeName(getParam("enzyme_name"));
			param.setEnzymeResidues(getParam("enzyme_residues"));
			param.setEnzymeBreakAA(getParam("enzyme_residues"), param.getMaxMissedCleavages(), param.isSemiCleavage());
			param.setEnzymeCut(getParam("enzyme_cut"));
			param.setEnzymeNocutResidues(getParam("enzyme_nocut_residues"));
			param.setEnzymeNoBreakAA(getParam("enzyme_nocut_residues"));

			if ("c".equals(getParam("enzyme_cut")))
				param.setEnzymeOffset(1);
			else
				// if("n".equals( getParam("enzyme_cut") )
				param.setEnzymeOffset(0);

			/*
			 * dN15Enrichment // Load ion series to consider, useA, useB, useY are for
			 * neutral losses iNumIonSeriesUsed = 0; for (int i = 0; i < 9; i++) {
			 * iIonVal10[i] = 10 * iIonVal[i]; iIonVal25[i] = 25 * iIonVal[i]; iIonVal50[i]
			 * = 50 * iIonVal[i]; if (iIonVal[i] > 0) {
			 * piSelectedIonSeries[iNumIonSeriesUsed++] = i; } } // differential mod. search
			 * for AAs listed in szDiffChar if (sParam.getDiffMass1() != 0.0 &&
			 * sParam.getDiffChar1().length() > 0) { sParam.isDiffSearch() = true; } if
			 * (sParam.getDiffMass2() != 0.0 && sParam.getDiffChar2().length() > 0) {
			 * sParam.isDiffSearch() = true; } if (sParam.getDiffMass3() != 0.0 &&
			 * sParam.getDiffChar3().length() > 0) { sParam.isDiffSearch() = true; } // if
			 * (sParam.isDiffSearch()) // iSizepiDiffSearchSites = MAX_PEPTIDE_SIZE*4;
			 * sParam.getPdAAMassParent()[AssignMass.S] += dN15Enrichment * addSmass;
			 * pdAAMassFragment[AssignMass.S] += dN15Enrichment * addSmass; } if (addPmass
			 * != 0.0) { sParam.getPdAAMassParent()[AssignMass.P] += dN15Enrichment *
			 * addPmass; pdAAMassFragment[AssignMass.P] += dN15Enrichment * addPmass; }
			 */

			if (null != getParam("diff_search_options")) {
				final String[] modArr = getParam("diff_search_options").split(" ");

				// Read each set of residue
				// e.g. diff_search_options = 80.0 ST -18.0 ST 0.0 X
				// *, #, @

				final Set<Double> modShiftSet = new THashSet<Double>();

				if (null != getParam("diff_search_options") && !"".equals(getParam("diff_search_options")))
					for (int i = 0; i < modArr.length; i += 2) {

						final double massShift = Double.parseDouble(modArr[i]);

						if (massShift != 0) {
							// System.out.println("fdsafsda");
							param.setDiffSearch(true);
							// symbolMap.put( MOD_SYMBOL[i/2], massShift);
							// this.isModification = true;

							for (int j = 0; j < modArr[i + 1].length(); j++) {

								final ModResidue residue = new ModResidue(modArr[i + 1].charAt(j), massShift);
								// param.addModHt(modArr[i + 1].charAt(j),
								// massShift);
								// DiffModification
								param.addModResidue(residue);

								DiffModification.setDiffModMass(residue.getResidue(), massShift);
								// System.out.println("fsdafsda" + " " +
								// residue.getResidue() + " " + massShift);
								// System.out.println("fsdafsda" + " " +
								// DiffModification.getDiffModMass('T'));
								modShiftSet.add(massShift);
							}
						}
					}

				// param.setModShiftSet(modShiftSet);
				// double[] modShiftArr = modShiftSet.toArray(new Double[0]);
				final ICombinatoricsVector<Double> initialVector = Factory.createVector(modShiftSet);

				for (int i = 1; i <= param.getMaxNumDiffMod(); i++) {
					final Generator<Double> gen = Factory.createMultiCombinationGenerator(initialVector, i);

					for (final ICombinatoricsVector<Double> combination : gen) {
						final List<Double> ml = combination.getVector();
						final List<Double> fl = new ArrayList<Double>();

						for (final Iterator<Double> mitr = ml.iterator(); mitr.hasNext();) {
							// System.out.println("----" + mitr.next());
							fl.add(mitr.next());
						}

						param.addModGroupList(fl);
					}
				}

				/*
				 * Set<Double> modShiftSet = sParam.getModShiftSet();
				 * //modShiftSet.toArray(pArr); ICombinatoricsVector<String> initialVector2 =
				 * Factory.createVector( new Integer[] { 79, 16} ); // Create a
				 * multi-combination generator to generate 3-combinations of // the initial
				 * vector gen = Factory.createMultiCombinationGenerator(initialVector, 2); //
				 * Print all possible combinations for (ICombinatoricsVector<String> combination
				 * : gen) { System.out.println(combination); } gen =
				 * Factory.createMultiCombinationGenerator(initialVector, 1); // Print all
				 * possible combinations for (ICombinatoricsVector<String> combination : gen) {
				 * System.out.println(combination); }
				 */

			}

			// added by Salva 24nov2014
			param.setMassGroupFactor(10000);
			param.setH2OPlusProtonAdded(true);
			//

		} catch (final IOException e) {
			System.out.println("Error reading param file " + e);
			throw new IOException(e.toString());
		} finally {
			if (null != br)
				;
			br.close();
		}
	}

	public double[] getMassShiftArr() {
		return massShiftArr;
	}

	public void pepProbeInit() {
	}

	public void gutentagInit() {
	}

	public SearchParams getSearchParams() {
		return param;
	}

	public int trimValueAsInt(String str, int defaultVal) {
		if (null == str) {
			return defaultVal;
		}
		int ret = defaultVal;
		try {
			ret = Integer.valueOf(trimValue(str));
		} catch (final NumberFormatException e) {
			logger.log(Level.SEVERE, "Blazmass parameter invalid, expecting int, got: " + str);
		}
		return ret;
	}

	public int trimValueAsInt(String str) {
		return trimValueAsInt(str, 0);
	}

	// public double trimValueAsDouble(String str) {
	// double ret = 0;
	// try {
	// ret = double.valueOf(trimValue(str));
	// } catch (NumberFormatException e) {
	// logger.log(Level.SEVERE,
	// "Blazmass parameter invalid, expecting double, got: " + str);
	// }
	// return ret;
	// }

	public double trimValueAsDouble(String str) {
		double ret = 0;
		try {
			ret = Double.valueOf(trimValue(str));
		} catch (final NumberFormatException e) {
			logger.log(Level.SEVERE, "Blazmass parameter invalid, expecting double, got: " + str);
		} catch (final NullPointerException e2) {
			logger.log(Level.SEVERE, "Blazmass parameter invalid, expecting double, got NULL (empty)");
		}
		return ret;
	}

	public String trimValue(String str) {
		if (null == str) {
			return null;
		}

		final int index = str.indexOf(';');
		if (index > 0) {
			str = str.substring(0, index);
		}

		return str.trim();
	}

	public Hashtable<String, String> getHashtable() {
		return ht;
	}

	public String getFileName() {
		return fileName;
	}

	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	public Hashtable<String, String> getHt() {
		return ht;
	}

	public void setHt(Hashtable<String, String> ht) {
		this.ht = ht;
	}

	public SearchParams getParam() {
		return param;
	}

	public void setParam(SearchParams param) {
		this.param = param;
	}

	public String getPath() {
		return path;
	}

	public void setPath(String path) {
		this.path = path;
	}

	/**
	 * Get a param value from the internal map, log warning if not present
	 *
	 * @param paramName name of the param to get
	 * @return the param value or empty string if not present
	 */
	private String getParam(String paramName) {
		// System.out.println("== " + paramName + "\t" + ht.get(paramName) + " "
		// + ht.contains(paramName));

		String paramValue = ht.get(paramName);
		if (paramValue == null) {
			System.out.println("warning: missing parameter: " + paramName);
		} else {
			paramValue = paramValue.trim();
		}

		return paramValue;
	}

	/**
	 * Get a param value from the internal map, log warning if not present
	 *
	 * @param paramName name of the param to get
	 * @return the param value or empty string if not present
	 */
	private String getParam(String paramName, boolean assureNotNullNumber) {
		// System.out.println("== " + paramName + "\t" + ht.get(paramName) + " "
		// + ht.contains(paramName));

		final String paramValue = getParam(paramName);
		if (assureNotNullNumber && paramValue == null) {
			return "0.0";
		}
		return paramValue;
	}

}
