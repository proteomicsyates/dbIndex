package edu.scripps.yates.dbindex.model;

/**
 * Different operational modes supported
 */
public enum IndexType {

	INDEX_NORMAL {
		@Override
		public String toString() {
			return "Normal index (best for small and medium db)";
		}

	},
	INDEX_LARGE {
		@Override
		public String toString() {
			return "Large index (best for large db)";
		}

	},
};
