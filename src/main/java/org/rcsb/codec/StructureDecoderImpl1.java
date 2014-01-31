/**
 * 
 */
package org.rcsb.codec;

import static org.rcsb.codec.CodecConstants.*;


import java.io.DataInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * StructureToBinary writes/reads the "Atom" section of a PDB/mmCIF file to/from a binary files.
 * 
 * Format:
 * magic number "PDBb" (4 bytes)
 * version string (4 bytes)
 * 
 * The following records have the general format:
 * 
 * Record id, record length in number of bytes, data ..
 *    
 * Byte number  0   1   2   3   4   5   6   7   8   9
 *            +---+---+---+---+---+---+---+---+---+---+--
 *            |id | record length | data ..
 *            +---+---+---+---+---+---+---+---+---+---+--
 * 
 *            | N |             4 | model count (int = 4 bytes)
 *            
 *            | M |             4 | chains count in model (int = 4 bytes)
 *            
 *            | C |             8 | chain id (4 bytes), group count (int = 4 bytes)
 *            
 *            | G |            12 | short group record: group name (3 bytes), group number (int = 4 bytes), insertion code (1 byte), atom count (int = 4 bytes)
 *            
 *            | G |            24 | long group record: group name (3 bytes), group number (4 bytes), insertion code (1 byte), atom count (4 bytes), 
 *                                  xOffset(4 bytes), yOffset (4 bytes), zOffset (4 bytes)
 *            
 *            | A |            17 | short atom record: atom name with PDB spacing (4 bytes), element string (2 bytes), 
 *                                  alternative location indicator (1 byte), deltaX (short = 2 bytes), deltaY (short = 2 bytes), deltaZ (short = 2 bytes), 
 *                                  occupancy (short = 2 bytes), temperature factor (short = 2 bytes)
 *            
 *            | A |            27 | long atom record: atom name with PDB spacing (4 bytes), element string (2 bytes), 
 *                                  alternative location indicator (1 byte), x (float = 4 bytes), y (4 bytes), z (4 bytes), 
 *                                  occupancy (4 bytes), temperature factor (4 bytes)
 *                                  
 * for short atom record:  x = deltaX * 0.001 + xOffset; deltaY * 0.001 + yOffset; deltaZ * 0.001 + zOffset; 
 *                         o = occupancy * 0.001, b = temperature factor * 0.001;
 *            
 * Record order in file:
 * 
 *  magic number
 *  version
 *  N
 *  M model 1
 *  C chain 1
 *    G long group record 1    OR     G short group record 1
 *      A short atom record 1           A long atom record 1
 *      A short atom record 2           A long atom record 1
 *      ...
 *    G group 2
 *      ...
 *  C chain 2
 *  ..
 *  M model 2
 *  ..
 * 
 *
 * @author Peter Rose
 *
 */
public class StructureDecoderImpl1 extends StructureDecoder {
	protected DataInputStream inStream = null;

	// byte arrays for temporary data
	private byte[] b1 = new byte[1];
	private byte[] b2 = new byte[2];
	private byte[] b3 = new byte[3];
	private byte[] b4 = new byte[4];
	private int[] i4 = new int[4];

	private int modelCount = 0;
	private int sequenceIndex = 0;
	private int chainCount = 0;
	private int groupNumber = 0;
	private int intType = 4;
	private int byteOffset = 0;
	
	private boolean homogeneousModel = false;
	private List<String> sequences = new ArrayList<String>();
	private List<Integer> sequenceIndices = new ArrayList<Integer>();
	private List<String> chainIds = null;
	private List<Integer> chainCounts = new ArrayList<Integer>();
	private List<Integer> groupCounts = new ArrayList<Integer>();
	private List<Integer> groupIndices = new ArrayList<Integer>();
	private List<Integer> groupNumbers = new ArrayList<Integer>();
	private List<String[]> groupInfo = new ArrayList<String[]>();
	private List<int[]> bondInfo = new ArrayList<int[]>();
	private List<Byte> flagInfo = new ArrayList<Byte>();
	private List<Integer> bFactors = new ArrayList<Integer>();
	private int[] occupancy = null;
	
	private boolean useBfactor = false;
	private boolean useOccupancy = false;

	private StructureInflatorInterface deflator;

	public StructureDecoderImpl1(DataInputStream dataInputStream, StructureInflatorInterface deflator) {
		this.inStream = dataInputStream;
		this.deflator = deflator;
	}
	
	public void decode() throws IOException {	
		modelCount = 0;
		chainCount = 0;
		homogeneousModel = false;
		
		byte recordId = 0;
		
		while ((recordId = (byte) inStream.readByte()) != END) {

			switch (recordId) {
			case STRUCTURE:
				readStructureRecord();
				break;
			case MODEL:
				readModelRecord();
				break;
			case SEQUENCE:
				readSequenceRecord();
				break;
			case CHAIN:
				readChainRecord();
				break;
			case GROUP:
				readGroupRecord();
				break;
			case GINFO:
				readGInfoRecord();
				break;
			case COORD:
				readCoordRecord();
				break;
			case BFACTOR:
				readBFactorRecord();
				break;
			case OCCUPANCY:
				readOccupancyRecord();
				break;
//			default:
//				if (recordId <= ZRecord) {
//					inStream.skipBytes(inStream.readByte());
//				} else {
//					inStream.skipBytes(inStream.readInt());
//				}
//				System.out.println("recordId: " + recordId + " skipping");
//				break;
			}
		}
	}
	
	private void readSequenceRecord() throws IOException {
		int len = inStream.readInt();
		String sequence = readFixedLengthString(len);
		sequences.add(sequence);
	}

	private void readStructureRecord() throws IOException {
	//	System.out.println("readStructureRecord");
		inStream.skipBytes(1);
		modelCount = inStream.readInt();	
		homogeneousModel = inStream.readBoolean();
	//	System.out.println("modelCount: " + modelCount);
	//	System.out.println("homogeneous: " + homogeneousModel);
		chainIds = new ArrayList<String>(chainCount);
		groupCounts = new ArrayList<Integer>(chainCount);
	}
	
	private void readModelRecord() throws IOException {
//		System.out.println("readModelRecord");
		// TODO this info needs to go into an array (size model number if not homogeneous)
		inStream.skipBytes(1);
		chainCounts.add(inStream.readInt());
//		System.out.println("chainCount: " + chainCounts.get(chainCounts.size()-1));
	}
	
	private void readChainRecord() throws IOException {
//		System.out.println("readChainRecord");
		inStream.skipBytes(1);
		sequenceIndex = inStream.readInt();
		sequenceIndices.add(sequenceIndex);
		String chainId = readFixedLengthString(4);
//		System.out.println("chainId: " + chainId);
		chainIds.add(chainId);
		int groupCount = inStream.readInt();	
		groupCounts.add(groupCount);
//		System.out.println("groupCount: " + groupCount);
		groupNumber = 0;	
	}
	
	private void readGroupRecord() throws IOException {
//		System.out.println("readGroupRecord");
		int len = inStream.readByte();
		int groupIndex = inStream.readInt();
		groupIndices.add(groupIndex);
//		System.out.println("groupIndex: " + groupIndex);
		if (len == 8) {
			groupNumber = inStream.readInt();
			groupNumbers.add(groupNumber);
//			System.out.println("groupNumber read: " + groupNumber);
		} else {
			groupNumbers.add(++groupNumber);
//			System.out.println("groupNumber calc: " + groupNumber);
		}
	}
	
	private void readGInfoRecord() throws IOException {
//		System.out.println("readGInfoRecord");
		int len = inStream.readInt();
		int aCount = inStream.readShort();
		len -= 4*aCount + 3;
		byte flags = inStream.readByte();	
		String[] info = readGroupInfo(len);
		int[] bondList = new int[aCount*2];
		for (int bl = 0; bl < aCount*2; bl++) {
			bondList[bl] = (int)inStream.readShort();
		}
		groupInfo.add(info);
		bondInfo.add(bondList);
		flagInfo.add(flags);
	}
	
	private void readBFactorRecord() throws IOException {
//		System.out.println("readBFactorRecord");
		int len = inStream.readInt();
		byteOffset = 0;
		intType = 4;
		while (byteOffset < len) {
			int val = readNextInt();
			bFactors.add(val);
		}
		useBfactor = true;
	}
	
	private void readOccupancyRecord() throws IOException {
	//	System.out.println("readOccupancyRecord");
		int len = inStream.readInt();	
		int n = len/2;
		occupancy = new int[n];
		for (int i = 0; i < n; i++) {
			occupancy[i] = inStream.readShort();
		}
		useOccupancy = true;
	}

	private void readCoordRecord() throws IOException {
		//		System.out.println("readCoordRecord");
	//	int len = inStream.readInt();
		inStream.skipBytes(4);
		int[] x = new int[1024];
		int[] y = new int[1024];
		int[] z = new int[1024];
		int[] b = new int[1024];

		deflator.setModelCount(modelCount);
		
		byteOffset = 0;
		intType = 4;

		int chainIndex = 0;
		int groupIndex = 0;
		int atomIndex = 0;

		for (int m = 0; m < modelCount; m++) {
			if (homogeneousModel) {
				chainIndex = 0;
				groupIndex = 0;
				chainCount = chainCounts.get(0);
			} else {
			     chainCount = chainCounts.get(m);
			}
			
			deflator.setModelInfo(m, chainCount);

			for (int i = 0; i < chainCount; i++) {			
				int xOffset = 0;
				int yOffset = 0;
				int zOffset = 0;
				int bOffset = 0;
				
				int xLink = 0;
				int yLink = 0;
				int zLink = 0;
				int bLink = 0;
				boolean hasTail = false;
						
				String chainId = chainIds.get(chainIndex);
				int seqIndex = sequenceIndices.get(chainIndex);
				int groupCount = groupCounts.get(chainIndex);
				
				deflator.setChainInfo(chainId, groupCount);
				chainIndex++;

				for (int j = 0; j < groupCount; j++) {
					groupNumber = groupNumbers.get(groupIndex);
					int gIndex = groupIndices.get(groupIndex);
					groupIndex++;

					String[] info = groupInfo.get(gIndex);
					int[] bondList = bondInfo.get(gIndex);
					int atomCount = bondList.length/2;
					System.out.println("atomCount: " + atomCount);
					byte flags = flagInfo.get(gIndex);

					boolean isAminoAcid = (flags & AMINO_ACID) != 0;
					boolean isNucleotide = (flags & NUCLEOTIDE) != 0;
				
					boolean hasHead =  (flags & HEAD) != 0;
					if (! hasTail) {
						xLink = 0;
						yLink = 0;
						zLink = 0;
						bLink = 0;
					}

					int index = 0;
					String groupName = info[index++];
					char insertionCode = info[index++].charAt(0);

					int polymerType = 0;
					if (isAminoAcid) {
						polymerType = 1;
					} else if (isNucleotide) {
						polymerType = 2;
					} 
					
					deflator.setGroupInfo(groupName, groupNumber, insertionCode, polymerType, atomCount);

					// TODO does this zero out the previous coords?
					if (atomCount > 1024) {
						x = new int[atomCount];
						y = new int[atomCount];
						z = new int[atomCount];
						b = new int[atomCount];
					}

					for (int k = 0; k < atomCount; k++) {
						String atomName = info[index++];
						String atomNameTrimmed = atomName.trim();
						String element = info[index++];
						char altLoc = info[index++].charAt(0);

	//					boolean isLink = xLink != 0 || yLink != 0 || zLink != 0 || bLink != 0;
						int distance = 0;
						if (bondList[k] >=0) {
							distance = bondList[atomCount+k];
						} else if (k == 0 && isAminoAcid && hasHead && hasTail) {
							distance = PEPTIDE_BOND_LENGTH;
						} else if (k == 0 && isNucleotide && hasHead && hasTail) {
							distance = NUCLEOTIDE_BOND_LENGTH;
						}

						int[] xyz = decodeCoords(distance);
						
//						if (xyz[2] == -14821) {
//							System.out.println(Arrays.toString(xyz));
//							System.out.println("distance: " + distance);
//						}
				//		System.out.println("decode: " + atomData[0] + "," + atomData[1] + "," + atomData[2]);
						if (bondList[k] >=0) {
							// use coordinates from a previous atom in this group
							int ind = bondList[k];
							xOffset = x[ind];
							yOffset = y[ind];
							zOffset = z[ind];
							bOffset = b[ind];
						} else if (k == 0 && hasTail && hasHead) {
							// TODO need to check if the first atom is really a link atom (could be missing, i.e., CA model)
							// take difference to previous linked residue
							// need to check if amino or nucleic
							xOffset = xLink;
							yOffset = yLink;
							zOffset = zLink;
							bOffset = bLink;
						}

						xOffset += xyz[0];
						yOffset += xyz[1];
						zOffset += xyz[2];
						if (useBfactor) {
							bOffset += bFactors.get(atomIndex);
						}
//						if (zOffset == 45887) {
//							System.out.println("isLink: " + isLink);
//							System.out.println("distance: " + distance);
//							System.out.println("link: " + xLink+","+yLink+","+zLink+","+bLink);
//							System.out.println("delta: " + xyz[0]+","+xyz[1]+","+xyz[2]+","+bFactors.get(atomIndex));
//							System.out.println("offset: " + xOffset+","+yOffset+","+zOffset+","+bOffset);
//							System.out.println("atomIndex: " + atomIndex);
//							System.out.println("bOffset: " + bOffset);
//						}
						int occ = B_PRECISION;
						if (useOccupancy) {
							occ = occupancy[atomIndex];
						}
						atomIndex++;

						x[k] = xOffset;
						y[k] = yOffset;
						z[k] = zOffset;
						b[k] = bOffset;

						deflator.setAtomInfo(atomName, 0, altLoc, xOffset*XYZ_SCALE, yOffset*XYZ_SCALE, zOffset*XYZ_SCALE, occ*B_SCALE, bOffset*B_SCALE, element);

						if ((isAminoAcid && atomNameTrimmed.equals("C")) || (isNucleotide && atomNameTrimmed.equals("O3'"))) {
							xLink = xOffset;
							yLink = yOffset;
							zLink = zOffset;
							bLink = bOffset;
						}	
					
					}
					
					hasTail =  (flags & TAIL) != 0;
				}
			}
		}
	}

	private void skipToNextRecord() throws IOException {
		int len = inStream.readInt();
		inStream.skipBytes(len);
	}
	
	private int readNextInt() throws IOException {
		int v = 0;
		switch (intType) {
		case 2:
			v = inStream.readShort();
			byteOffset += 2;
			
			switch (v) {
			case BYTE2_MARKER4: 
				intType = 4;
				return readNextInt();
			case BYTE2_MARKER5: 
				intType = 5;
				return readNextInt();
			}
			break;

		case 4:
			v = inStream.readInt();
			byteOffset += 4;

			switch (v) {
			case BYTE4_MARKER2: 
				intType = 2;
				return readNextInt();
			case BYTE4_MARKER5: 
				intType = 5;
				return readNextInt();
			}
			break;
			
		case 5:
			v = inStream.readInt();
			byteOffset += 4;

			switch (v) {
			case BYTE4_MARKER2: 
				intType = 2;
				return readNextInt();
			case BYTE4_MARKER4: 
				intType = 4;
				return readNextInt();
			}
			break;
		}

		return v;
	}
	
	private int[] decodeCoords(int distance) throws IOException {
		int v = readNextInt();
		if (intType == 5) {
			BitEncoder.fromInt(v, distance, b4, i4);
		} else {
			i4[0] = v;
			i4[1] = readNextInt();
			i4[2] = readNextInt();
		}
		return i4;
	}

	private String[] readGroupInfo(int length) throws IOException {
		int len = 2 + 3 * (length - 4)/7;
		String[] info = new String[len];

		inStream.read(b3); 
		info[0] = new String(b3); // group name
		inStream.read(b1); 
		info[1] = new String(b1); // insertion code
		for (int i = 2; i < len; i += 3) {
			inStream.read(b4);
			info[i] = new String(b4); // atom name
			inStream.read(b2);
			info[i+1] = new String(b2).trim(); // element name
			inStream.read(b1);
			info[i+2] = new String(b1); // alternative location
		}
		return info;
	}

	private String readFixedLengthString(int length) throws IOException {
		byte[] bytes = new byte[length];
		inStream.read(bytes);
		return new String(bytes);
	}
}

