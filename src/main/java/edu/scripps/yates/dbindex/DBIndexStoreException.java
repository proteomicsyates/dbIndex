package edu.scripps.yates.dbindex;

/**
 * Exception thrown when error accessing the db store occurred
 */
public class DBIndexStoreException extends Throwable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 750954811615264144L;

	public DBIndexStoreException(String message) {
		super(message);
	}

	public DBIndexStoreException(String message, Throwable cause) {
		super(message, cause);
	}

	public DBIndexStoreException(Throwable cause) {
		super(cause);
	}
}
