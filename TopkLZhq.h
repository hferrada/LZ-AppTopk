/*
 * TopkLZhq.h
 *
 *  Created on: 10-09-2014
 *      Author: hector
 */

#ifndef TOPKLZHQ_H_
#define TOPKLZHQ_H_

#include "LZ78Tries64.h"
#include <RangeMMTree64.h>

using namespace sdsl;
using namespace std;
using namespace drf64;

class TopkLZhq {
private:
	ulong lgD;				// Binary logarithm (ceiling) of nDocs
	char boundSymbol;		// original symbol delimiter of documents when we read all documents in 1 file.
	ulong* docLength;		// store the length for each document
	ulong* charTf;			// store the tf for each symbol. This is no include in the final structure

	ulong nod;				// number of phrase nodes nod. Here nod = n'+1, n' of LZIndex
	uint lgNod;				// ceiling binary log of nod
	ulong nF;				// number of fictitious nodes in RevTrie, nF only include one fictitious node for each unary path of fictitious nodes
	ulong nEF;				// this is the difference between the total fictitious nodes in the original RevTrie an the RevTrie represented where we join the fictitious nodes in unary path
	uint nEFHead;			// number of fictitious nodes are head of a unary path
	ulong nodRev;			// nod + nF
	ulong *PLz;				// LZTrie's DFUDS sequence representation in 2*nod bits
	ulong nLz;				// length of PLz DFUDS sequences, these have exactly nLz = 2*nod nodes
	ulong nRev;				// length of PRev DFUDS sequences, these have exactly nRev = 2(nod+nF) nodes
	ulong *PRev;			// RevTrie's DFUDS sequence representation in 2(nod+nF) bits

	// Labels...
	uchar *LbRev;			// RevTrie's labels, for the n' phrase nodes, in DFUDS representation
	bit_vector fRev_b;		// this marks the fictitious nodes in LbRev
	rrr_vector<127> fRev_rrr;
	rrr_vector<127>::rank_1_type fRev_rank;
	uchar *LbRevF;			// RevTrie's labels for fictitious nodes with unary paths
	bit_vector fictU_b;		// this marks the fictitious nodes with unary path in FRev
	rrr_vector<127> fictU_rrr;
	rrr_vector<127>::rank_1_type fictU_rank;
	bit_vector listF_b;		// Delimiter for each unary sequence in LbRevF
	rrr_vector<127> listF_rrr;
	rrr_vector<127>::select_1_type lisF_sel;

	RangeMMTree64 *treeLz;	// range min-max tree representations for PLz
	RangeMMTree64 *treeRev;	// range min-max tree representations for PRev

	ulong *DocLZ;			// Document array to phrases with length nod-1
	ulong *Node;			// This array give the preorder value in LZTrie, that is: Range[j] = i <-> the node with preorder j in RevTrie corresponds to the node with preorder i in LZTrie

	// MARKS NODES...
	ulong nTopk;			// number of nodes with topk list
	bit_vector bTopk_b;		// marks all node that stores a topk List, length: [1..n']. This has nTopk 1's and z 0's
	rrr_vector<127> bTopk_rrr;
	rrr_vector<127>::rank_1_type bTopk_rank;

	bit_vector bDL_b;		// for the bTopk nodes marked, it marks the nodes that store all its different documents

	// TOP-K LIST...
	ulong *docTopk;			// this are the nTopk lists of id documents, with log(D) bits per item, maximum k id per list
	bit_vector limDTk_b;	// this marks the boundaries for each list in docTopk
	rrr_vector<127> limDTk_rrr;
	rrr_vector<127>::rank_1_type limDTk_rank;
	rrr_vector<127>::select_1_type limDTk_sel;

	// sort the list in increasing mode by frequency
	void sortTopkDocsByFq(uint* fNod, uint k, uint* dNodK, uint* fNodK);

	// set boundaries between each sublist of docId stored in docTopk
	void marksBoundaries(LZ78Tries64::RevNode* node, ulong *pre, ulong *posD);

	// It stores documents and generates boundaries of documents for each node marked in RevTrie
	// posDHI is the total number of increasing list
	//void storeDocuments_old(ulong* fNod, LZ78Tries64::RevNode* nod, ulong *pre, ulong *posD, ulong *nodTopk, ulong* dNodK, ulong* fNodK);

	void create_sNode(LZ78Tries64::RevNode* nod, ulong *pre, ulong *aux_sNode, ulong *posNod);

	// search the current pattern (given the cross of the RevTrie) in order to construct the real top-k answer
	bool searchPatternInFMI(uchar* pat, uint m, uint *fqCSA);

	// it counts and stores in 'posBit' the total length frequencies (in bits)
	void countLengthFrequencies(uint *fNod, LZ78Tries64::RevNode* nod, ulong *pre, ulong *posD, ulong *nodTopk, uint* dNodK, uint* fNodK, ulong *docTopkAux);
	//void countLengthFrequencies_old(ulong* fNod, LZ78Tries64::RevNode* nod, ulong *pre, ulong *posD, ulong *nodTopk);

	// determine the doc Id for the text position
	uint searchDocPosition(ulong pos);

	// count the frequencies for all occurrence in this node (node for a true phrase or for a fictitious node)
	ulong countMyFrequencies(LZ78Tries64::RevNode* nod, ulong* fNod);
	ulong countMyFrequencies(LZ78Tries64::RevNode* nod);

	// count the length of its document segments that includes this node and also the number of different documents
	ulong nodeWithTopkList(LZ78Tries64::RevNode *nod);

	// marks the nodes that store a topk list in bTopk[]. Set the bits in bBlockLZ[], and filled the array fatherLZ[1..n']
	void marksTopkNodes(LZ78Tries64::RevNode *node, ulong *pre, uint childFather);

	// create the document array for LZTrie's phrases
	void genDocArray(LZ78Tries64::LZNode* nod, ulong *DocLZ, ulong *ini);

	//	create 	LbRevF	:	RevTrie's labels for fictitious nodes with unary paths
	// 			ListF	:	Delimiter for each unary sequence in LbRevF
	// 			fRev	:	this marks the fictitious nodes in LbRev (RavTrie's labels)
	// 			FictU	:	this marks the fictitious nodes with unary path in FRev
	void genSeqExtraFictRevTrie(LZ78Tries64::RevNode *node, ulong *posExtF, ulong *pfRev, ulong *posFU);

	//	create 	PRev, LbRev and LbRevF structures
	void genDFUDSSeqRevTrie(LZ78Tries64::RevNode *node, ulong *pos, ulong *posRev);

	//	create 	PLZ	:	LZTrie's DFUDS representation without Labels
	void genDFUDSSeqLZTrieNotLbl(LZ78Tries64::LZNode *node, ulong *pos);

	void createStructure(char dirSaveLoad[300]);
	void readListFiles(char *inputFile, bool bitLowercase);
	void readListFilesTrec(char *inputFile, ulong max_size);			// ELIMINAR
	void readUniqueFile(char *inputFile, bool bitLowercase);
	void makePriorities(ulong* cPosTb);

	// the doc with key dictD[k] go up until its correct position in the maximum queue
	void swimInMaxQ(uint k);

	// set the new item at top of the queue heap[1..k]
	void setTopMaxQ(uint k);

public:
	static bool TRACE;		// true: print all details for console
	static bool TEST;		// true: print all details for console

	ulong sizeDS;			// size in bytes for all data structure

	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	uint nDocs;				// Number of distinct documents D
	char cutDoc;			// symbol to separate documents
	uint g;					// minimum length of document segments to store a topk list in a node
	LZ78Tries64 *tries;		// structure to store the original LZTrie and RevTrie.

	ulong* EndDocs;			// this stores the final phrase number for each document. This is no include in the final structure
	ulong *EndDocsPos;		// this stores the final text position for each document. This is no include in the final structure

	// BRUTE FORCE FOR COUNT AND SORT FREQUENCIES...
	// ** IMPORTAT: dictD only can store keys of length of uint (16 bits)
	ulong *dictD;			// dictionary of keys, length: [1..D] with log(k*) bits each id-Doc
	ulong *keyDoc;			// document of each key, length: [1..k*] with lgD bits
	ulong *keyFq;			// frequencies for each key, length: [1..k*] with 8 bits
	uint *heap;				// maximum priority queue for keyFq, length: [1..k*]

	bool *dictDst;			// documents stored before of brute force search.

	char dirStore[200];		// directory to save/load the data structure (files *.tk)


	// this bit_vector is only for test search of patterns and not belong to the index!!
	bit_vector separator_b;		// This marks the beginning of each LZ-phrase
	rrr_vector<127> sep_rrr;
	rrr_vector<127>::rank_1_type sep_rank;
	rrr_vector<127>::select_1_type sep_sel;
	sdsl::csa_wt<wt_huff<rrr_vector<127> >, 4, 4> fmi;

	TopkLZhq(uint gVal, char *inputFile, bool filesInList, char cutDocCode, bool bitLowercase, char dirSaveLoad[300], char boundS);
	TopkLZhq(uint gVal, char *inputFile, bool fileType, char cutDocCode, bool bitLowercase, char dirSaveLoad[300], char boundS, ulong max_size);  // ELIMINAR
	TopkLZhq(char dirSaveLoad[200], bool showSize);
	virtual ~TopkLZhq();

	// It stores in docList[1..k] and frqList[1..k] the list and frequencies of the top-k most frequent documents for the pattern p[1..m]
	ulong topKDocument(uchar* pat, uint m, uint* docList, uint k);

	// save the Data Structure in folder 'dirStore'
	void saveDS(bool showSize);

	// load the Data Structure from the folder 'fileName'
	void loadDS(bool showSize);

	// determine the doc Id for the phrase idNode
	uint searchDocument(ulong idNode);

	// return true if the pattern is in the RevTrie and in '*pv' the respective preorder value in the tree,
	// in '*x' stores its DFUDS bit position. Return false if the pattern does not exist
	bool searchPattern_Rev(uchar* pat, uint m, ulong *x, ulong *pv);
};

#endif /* TOPKLZHQ_H_ */
