package edu.scripps.yates.dbindex;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;

import org.apache.log4j.Logger;

/**
 * Dynamic byte buffer utility
 * 
 * @author Adam
 */
public class DynByteBuffer {
	private static int DEFAULT_CAPACITY = 4 * 4; // 1 sequence = 4 * 4 bytes
	private int curSize = 0;
	private int curCapacity = 16;
	private byte[] data;
	private static Logger log = Logger.getLogger(DynByteBuffer.class);
	private static int maxCapacity = 0;

	public DynByteBuffer() {
		data = new byte[curCapacity];
	}

	public void add(byte toAdd) {
		final byte[] array = new byte[1];
		array[0] = toAdd;
		add(array);
	}

	public void add(byte[] toAdd) {
		final int toAddLen = toAdd.length;
		if (curCapacity < curSize + toAddLen) {
			final int newCapacity = (curCapacity + toAddLen) * 3;
			if (newCapacity > maxCapacity) {
				log.debug("newCapacity= " + newCapacity);
				maxCapacity = newCapacity;
			}
			final byte[] newData = new byte[newCapacity];
			System.arraycopy(data, 0, newData, 0, curSize);
			curCapacity = newCapacity;
			data = newData;
		}

		// insert new data
		System.arraycopy(toAdd, 0, data, curSize, toAddLen);
		curSize += toAddLen;
	}

	public void clear() {
		curSize = 0;
		curCapacity = DEFAULT_CAPACITY * 4;
		data = new byte[curCapacity];

	}

	public int getSize() {
		return curSize;
	}

	public byte[] getData() {
		final byte[] copyOfRange = Arrays.copyOfRange(data, 0, curSize);
		return copyOfRange;
	}

	public static byte[] toByteArray(int myInteger) {
		return ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN).putInt(myInteger).array();
	}

	public static int toInt(byte[] byteArray) {
		return ByteBuffer.wrap(byteArray).order(ByteOrder.LITTLE_ENDIAN).getInt();
	}

	public static byte[] toByteArray(char myChar) {
		return ByteBuffer.allocate(2).order(ByteOrder.LITTLE_ENDIAN).putChar(myChar).array();
	}

	public static char toChar(byte[] byteArray) {
		return ByteBuffer.wrap(byteArray).order(ByteOrder.LITTLE_ENDIAN).getChar();
	}

	public static byte[] toByteArray(short myShort) {
		return ByteBuffer.allocate(2).order(ByteOrder.LITTLE_ENDIAN).putShort(myShort).array();
	}

	public static short toShort(byte[] byteBarray) {
		return ByteBuffer.wrap(byteBarray).order(ByteOrder.LITTLE_ENDIAN).getShort();
	}

	// public static byte[] toByteArray(float theFloat) {
	// return ByteBuffer.allocate(4).order(ByteOrder.LITTLE_ENDIAN)
	// .putFloat(theFloat).array();
	// }

	public static byte[] toByteArray(double theDouble) {
		return ByteBuffer.allocate(8).order(ByteOrder.LITTLE_ENDIAN).putDouble(theDouble).array();
	}

	// public static float toFloat(byte[] byteBarray) {
	// return ByteBuffer.wrap(byteBarray).order(ByteOrder.LITTLE_ENDIAN)
	// .getFloat();
	// }

	public static double toDouble(byte[] byteBarray) {
		return ByteBuffer.wrap(byteBarray).order(ByteOrder.LITTLE_ENDIAN).getDouble();
	}

	public static int toInt(byte myByte) {
		return myByte & (0xff);
	}
}
