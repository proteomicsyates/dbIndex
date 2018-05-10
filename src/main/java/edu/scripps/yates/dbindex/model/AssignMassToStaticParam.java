package edu.scripps.yates.dbindex.model;

import edu.scripps.yates.dbindex.SearchParams;

public class AssignMassToStaticParam {

	public static void addMassAndStaticParam(int i, double mass) {

		if (mass <= 0)
			return;
		SearchParams.addStaticParam((char) i, mass);

		edu.scripps.yates.utilities.masses.AssignMass.addMass(i, mass);

	}

}
