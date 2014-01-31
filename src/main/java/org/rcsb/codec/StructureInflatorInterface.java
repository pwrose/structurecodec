package org.rcsb.codec;

public interface StructureInflatorInterface {

	void setModelCount(int modelCount);
	void setModelInfo(int modelNumber, int chainCount);
	void setChainInfo(String chainId, int groupCount);
	void setGroupInfo(String groupName, int groupNumber, char insertionCode, int polymerType, int atomCount);
	void setAtomInfo(String atomName, int serialNumber, char alternativeLocationId, 
			float x, float y, float z, float occupancy, float temperatureFactor, String element);
	
}
