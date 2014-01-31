package org.rcsb.codec;

import java.util.Arrays;

public class BitEncoder {
	private static final int BYTE7_MIN_VALUE  =   -64;
	private static final int BYTE7_MAX_VALUE  =    63;
	private static final int BYTE12_MIN_VALUE = -2048;
	private static final int BYTE12_MAX_VALUE =  2047;	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {		
//		int distance = 1502;
//		int deltaX = 1018;
//		int deltaY = -916;
//		int deltaZ = -616;
		
		int distance = 1325;
		int deltaX = 901;
		int deltaY = 894;
		int deltaZ = -398;
		
		
		int[] outVal = new int[1];
		BitEncoder.toInt(distance, deltaX, deltaY, deltaZ, outVal);
		System.out.println("new val: " + outVal[0]);
		byte[] buffer = new byte[4];
		int[] out = new int[4];
		fromInt(outVal[0], distance, buffer, out);
		System.out.println("out: " + Arrays.toString(out));
	}
	
	public static boolean toInt(int distance, int x, int y, int z, int[] out) {
		if (x < BYTE12_MIN_VALUE || x > BYTE12_MAX_VALUE) {
			return false;
		}
		if (y < BYTE12_MIN_VALUE || y > BYTE12_MAX_VALUE) {
			return false;
		}

		int delta = (int)Math.round(Math.sqrt(distance*distance - x*x - y*y)) - Math.abs(z);
		if (delta < BYTE7_MIN_VALUE || delta > BYTE7_MAX_VALUE) {
			return false;
		}
		
		int sign = 0; // -1 -> 1b; 0 -> 0b;
		if (z < 0) {
			sign = -1;
		}
		out[0] = (int)compressSigned4(x, y, delta, sign);
		
		return true;
	}
	
	public static void fromInt(int value, int distance, byte[]buffer, int[] out) {
		longTo4Bytes(value, buffer);
		from4Bytes(buffer, distance, out);
	}
	
	private static void from4Bytes(byte[] bytes, int distance, int[] out) {
		long c4 = bytes4ToLong(bytes);
		decompressSigned4(c4, out);

		int z = (int)Math.round(Math.sqrt(distance*distance - out[0]*out[0] - out[1]*out[1])) - out[2];
		if (out[3] == -1) {
			z *= -1;
		}
		out[2] = z;
	}
	
	private static long bytes4ToLong(byte[] b) {
		long v = b[0] & 0xFF;
		v |= ((long)(b[1] & 0xFF) << 8);
		v |= ((long)(b[2] & 0xFF) << 16);
		v |= ((long)(b[3] & 0xFF) << 24);
		return v;		
	}
	
	private static long compressSigned4(int i1, int i2, int i3, int i4) {
		return compress4(toUnsignedInt(i1),toUnsignedInt(i2),toUnsignedInt(i3),toUnsignedInt(i4));
	}
	
	private static long compress4(int i1, int i2, int i3, int i4) {
		return (long)i1 | (((long)i2) << 12) | (((long)i3) << 24) | (((long)i4) << 31);
	}
	
	private static void decompressSigned4(long value, int[] out) {
		decompress4Special(value, out);
		out[0] = toSignedInt(out[0]);
		out[1] = toSignedInt(out[1]);
		out[2] = toSignedInt(out[2]);
		out[3] = toSignedInt(out[3]);
	}
	
	private static void decompress4Special(long value, int[] out) {
		out[0] = (int)((value & 0xFFF));
		out[1] = (int)((value >> 12) & 0xFFF);
		out[2] = (int)((value >> 24) & 0x7F);
		out[3] = (int)((value >> 31));
	}

	private static void longTo4Bytes(long value, byte[] bytes) {
		bytes[0] = (byte)(value >>  0); 
		bytes[1] = (byte)(value >>  8);
		bytes[2] = (byte)(value >> 16);
		bytes[3] = (byte)(value >> 24); 	
	}
	
	private static int toUnsignedInt(final int n) {
	    return (n << 1) ^ (n >> 31);
	 }

	private static int toSignedInt(final int n) {
	    return (n >>> 1) ^ -(n & 1);
    }
	
}
