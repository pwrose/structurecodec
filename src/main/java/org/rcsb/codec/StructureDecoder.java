/**
 * 
 */
package org.rcsb.codec;

import java.io.DataInputStream;
import java.io.IOException;

/**
 * @author Peter
 *
 */
public abstract class StructureDecoder {
	/**
	 * Returns a structure decoder for supported versions and compression levels
	 * @param majorVersion
	 * @param minorVersion
	 * @param compressionLevel
	 * @param dataInputStream
	 * @param inflator
	 * @return structure decoder
	 * @throws IOException
	 */
	public static StructureDecoder getDecoder(byte majorVersion, byte minorVersion, int compressionLevel, DataInputStream dataInputStream, StructureInflatorInterface inflator) throws IOException {
        if (majorVersion == 0 && minorVersion == 0 && compressionLevel == 1) {
			return new StructureDecoderImpl1(dataInputStream, inflator);
		}
        throw new IOException("StructureDecoder: invalid version or compression level: major version: " + 
		majorVersion + " minor version: " + minorVersion + " compression level: " + compressionLevel);
	}
	
	public abstract void decode() throws IOException;
}
