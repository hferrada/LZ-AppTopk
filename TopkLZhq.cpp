/*
 * TopkLZhq.cpp
 *
 *  Created on: 10-09-2014
 *      Author: hector
 */

#include "TopkLZhq.h"
bool TopkLZhq::TRACE = false;
bool TopkLZhq::TEST = false;
const uint MAX_DEP = 100000;
const uint MAX_LEN_DOCS = 40; // it is 40 bpc
uchar auxPath[MAX_DEP];	// this is used to store the labels of unary path at construction

ulong auxLenD;

TopkLZhq::TopkLZhq(char dirSaveLoad[200], bool showSize) {
	strcpy(dirStore, dirSaveLoad);
	loadDS(showSize);

	ulong size = nDocs+1;
	this->dictD = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	if (showSize) cout << " ** size of dictD[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->dictDst = new bool[size];
	size *= sizeof(bool);
	sizeDS += size;
	if (showSize) cout << " ** size of dictDst[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyDoc = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	if (showSize) cout << " ** size of keyDoc[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyFq = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	if (showSize) cout << " ** size of keyFq[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->heap = new uint[size];
	size *= sizeof(uint);
	sizeDS += size;
	if (showSize) cout << " ** size of heap[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

}

// ===========================================================================================
// only for TREC files !!
// ===========================================================================================
TopkLZhq::TopkLZhq(uint gVal, char *inputFile, bool fileType, char cutDocCode, bool bitLowercase, char dirSaveLoad[300], char boundS, ulong max_size) {
	ulong i;

	// size for variables
	this->sizeDS = 9*sizeof(ulong) + 4*sizeof(uint) + sizeof(char);

	this->g = gVal;
	this->cutDoc = cutDocCode;
	this->boundSymbol = boundS;
	strcpy(dirStore, dirSaveLoad);
	this->nDocs = this->n = this->nEFHead = 0;

	charTf = new ulong[LZ78Tries64::SIGMA];
	for(i=0; i<LZ78Tries64::SIGMA; i++)
		charTf[i] = 0;
	if (fileType){
		readListFilesTrec(inputFile, max_size); // 1GB = 1073741824 bytes			// ELIMINAR solo para trec
	}else
		readUniqueFile(inputFile, bitLowercase);

	// create the FMI, needed at time construction, to compute the real top-k list for market nodes.
	char fileName[300];
	cout << "____________________________________________________" << endl;
	cout << " Make fmi (for top-k brute force)..." << endl;
	cout << " Reading... " << inputFile << endl;
	strcpy(fileName, "");
	strcpy(fileName, inputFile);
	strcat(fileName, "_copy.txt");
	construct(fmi, fileName, 1); // generate index
	cout << " **  FMI size " << size_in_bytes(fmi) << " bytes = " << (float)size_in_bytes(fmi)/(float)n << "|T|" << endl;
	cout << " **  FMI length " << fmi.size() << endl;
	if (fmi.size() != n+1){
		cout << "ERROR. FMI length != n = " << n << endl;
		exit(1);
	}
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fmi.test");
	store_to_file(fmi, fileName);

	ulong* cPosTb = new ulong[LZ78Tries64::SIGMA]; // hasTable position for each character
	makePriorities(cPosTb);

	this->EndDocs = new ulong[nDocs];
	this->EndDocsPos = new ulong[nDocs];

	//cout << "  Create Original LZTrie and RevTrie..." << endl;
	tries = new LZ78Tries64(nDocs, cutDoc, cPosTb);
	this->separator_b = bit_vector(n, 0);

	tries->createTriesFromText(seq, n, &Node, EndDocs, EndDocsPos, &separator_b);
	cout << "  LZTrie with " << tries->nPhra << " nodes (phrases)" << endl;
	sizeDS += tries->sizeNode;
	cout << " ** Size of Node array " << tries->sizeNode << " = " << tries->sizeNode*8.0/(float)n << " bpc" << endl;

	sep_rrr = rrr_vector<127>(separator_b);
	sep_rank = rrr_vector<127>::rank_1_type(&sep_rrr);
	sep_sel = rrr_vector<127>::select_1_type(&sep_rrr);

	tries->countFictNodRevTrie(tries->revTrie);	// this set tries->nFictU
	cout << "  RevTrie with " << tries->nPhra << " nodes, " << tries->nFictU << " nFictU nodes and " << tries->nExtraFict << " nExtraFict nodes" << endl;

	if (TRACE){
		cout << endl << " LZTrie with n'= " << tries->nPhra << " nodes..." << endl;
		//tries->listLZTrie();
		cout << "====================================================" << endl;
		cout << endl << " RevTrie with n'= " << tries->nPhra << " phrase nodes, nFictU = " << tries->nFictU << " fictitious nodes, and " << tries->nExtraFict << " extra fictitious nodes (nExtraFict)... " << endl;
		cout << " Node[0.." << tries->nPhra-1 << "]..." << endl;

		ulong len = separator_b.bit_size();
		cout << " separator_b[0.." << len-1 << "]:" << endl;
		cout << separator_b[0];
		for(i=1; i<len; i++)
			cout << separator_b[i];
		cout << endl;
	}

	this->lgD = ceilingLog64(nDocs, 2);
	this->nod = tries->nPhra;
	this->lgNod = ceilingLog64(nod, 2);
	this->nF = tries->nFictU;
	this->nEF = tries->nExtraFict;
	this->nodRev = nod+nF;
	this->nRev = 2*nodRev;
	this->nLz = 2*nod;

	cout << "  DFUFS RevTrie seq with " << nRev << " bits" << endl;
	cout << "  DFUFS LZTrie seq with " << nLz << " bits" << endl;

	ulong size = nLz/W64;
	if (nLz%W64)
		size++;
	this->PLz = new ulong[size];

	size = nRev/W64;
	if (nRev%W64)
		size++;
	this->PRev = new ulong[size];

	this->LbRev = new uchar[nodRev];
	cout << " length of LbRev  " << nodRev << endl;
	size = nodRev*sizeof(uchar);
	sizeDS += size;
	cout << " ** size of LbRev[]: " << size << " = " << size*8.0/(float)n << " bpc" << endl;
	cout << " LbRev lenght " << nodRev << endl;

	this->LbRevF = new uchar[nEF];
	//cout << "## len of LbRevF  " << nEF << endl;
	size = nEF*sizeof(uchar);
	sizeDS += size;
	cout << " ** size of LbRevF[]: " << size << " = " << size*8.0/(float)n << " bpc" << endl;
	cout << " LbRevF lenght " << nEF << endl;

	size = nDocs+1;
	this->dictD = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of dictD[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->dictDst = new bool[size];
	size *= sizeof(bool);
	sizeDS += size;
	cout << " ** size of dictDst[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyDoc = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of keyDoc[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyFq = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of keyFq[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->heap = new uint[size];
	size *= sizeof(uint);
	sizeDS += size;
	cout << " ** size of heap[1..nDocs]" << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nod*lgD;
	if(size%W64)
		size = size/W64 + 1;
	else
		size /= W64;
	this->DocLZ = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of DocLZ[1..n']: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;
	cout << "====================================================" << endl << endl;

	// set initial open parenthesis in DFUDS sequences
	setBit64(PLz, 0);
	setBit64(PRev, 0);
	createStructure(dirSaveLoad);
}

TopkLZhq::TopkLZhq(uint gVal, char *inputFile, bool filesInList, char cutDocCode, bool bitLowercase, char dirSaveLoad[300], char boundS) {
	ulong i;
	char fileName[300];

	// size for variables
	this->sizeDS = 9*sizeof(ulong) + 4*sizeof(uint) + sizeof(char);

	this->g = gVal;
	this->cutDoc = cutDocCode;
	this->boundSymbol = boundS;
	strcpy(dirStore, dirSaveLoad);
	this->nDocs = this->n = this->nEFHead = 0;

	charTf = new ulong[LZ78Tries64::SIGMA];
	for(i=0; i<LZ78Tries64::SIGMA; i++)
		charTf[i] = 0;
	if (filesInList){
		readListFiles(inputFile, bitLowercase);
	}else
		readUniqueFile(inputFile, bitLowercase);

	ulong* cPosTb = new ulong[LZ78Tries64::SIGMA]; // hasTable position for each character
	makePriorities(cPosTb);

	// create the FMI, needed at time construction, to compute the real top-k list for market nodes.
	cout << "____________________________________________________" << endl;
	cout << " Make fmi (for top-k brute force)..." << endl;
	cout << " Reading... " << inputFile << endl;
	strcpy(fileName, "");
	strcpy(fileName, inputFile);
	cout << " Reading... " << inputFile << endl;
	strcat(fileName, "_copy.txt");
	construct(fmi, fileName, 1); // generate index
	cout << " **  FMI size " << size_in_bytes(fmi) << " bytes = " << (float)size_in_bytes(fmi)/(float)n << "|T|" << endl;
	cout << " **  FMI length " << fmi.size() << endl;
	if (fmi.size() != n){
		cout << "ERROR. FMI length != n = " << n << endl;
		exit(1);
	}
	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fmi.test");
	store_to_file(fmi, fileName);

	this->EndDocs = new ulong[nDocs];
	this->EndDocsPos = new ulong[nDocs];

	//cout << "  Create Original LZTrie and RevTrie..." << endl;
	tries = new LZ78Tries64(nDocs, cutDoc, cPosTb);
	this->separator_b = bit_vector(n, 0);

	tries->createTriesFromText(seq, n, &Node, EndDocs, EndDocsPos, &separator_b);
	cout << "  LZTrie with " << tries->nPhra << " nodes (phrases)" << endl;
	sizeDS += tries->sizeNode;
	cout << " ** Size of Node array " << tries->sizeNode << " = " << tries->sizeNode*8.0/(float)n << " bpc" << endl;

	sep_rrr = rrr_vector<127>(separator_b);
	sep_rank = rrr_vector<127>::rank_1_type(&sep_rrr);
	sep_sel = rrr_vector<127>::select_1_type(&sep_rrr);

	tries->countFictNodRevTrie(tries->revTrie);	// this set tries->nFictU
	cout << "  RevTrie with " << tries->nPhra << " nodes, " << tries->nFictU << " nFictU nodes and " << tries->nExtraFict << " nExtraFict nodes" << endl;

	if (TRACE){
		cout << endl << " LZTrie with n'= " << tries->nPhra << " nodes..." << endl;
		//tries->listLZTrie();
		cout << "====================================================" << endl;
		cout << endl << " RevTrie with n'= " << tries->nPhra << " phrase nodes, nFictU = " << tries->nFictU << " fictitious nodes, and " << tries->nExtraFict << " extra fictitious nodes (nExtraFict)... " << endl;
		cout << " Node[0.." << tries->nPhra-1 << "]..." << endl;

		ulong len = separator_b.bit_size();
		cout << " separator_b[0.." << len-1 << "]:" << endl;
		cout << separator_b[0];
		for(i=1; i<len; i++)
			cout << separator_b[i];
		cout << endl;
	}

	this->lgD = ceilingLog64(nDocs, 2);
	this->nod = tries->nPhra;
	this->lgNod = ceilingLog64(nod, 2);
	this->nF = tries->nFictU;
	this->nEF = tries->nExtraFict;
	this->nodRev = nod+nF;
	this->nRev = 2*nodRev;
	this->nLz = 2*nod;

	cout << "  DFUFS RevTrie seq with " << nRev << " bits" << endl;
	cout << "  DFUFS LZTrie seq with " << nLz << " bits" << endl;

	ulong size = nLz/W64;
	if (nLz%W64)
		size++;
	this->PLz = new ulong[size];

	size = nRev/W64;
	if (nRev%W64)
		size++;
	this->PRev = new ulong[size];

	this->LbRev = new uchar[nodRev];
	cout << " length of LbRev  " << nodRev << endl;
	size = nodRev*sizeof(uchar);
	sizeDS += size;
	cout << " ** size of LbRev[]: " << size << " = " << size*8.0/(float)n << " bpc" << endl;
	cout << " LbRev lenght " << nodRev << endl;

	this->LbRevF = new uchar[nEF];
	//cout << "## len of LbRevF  " << nEF << endl;
	size = nEF*sizeof(uchar);
	sizeDS += size;
	cout << " ** size of LbRevF[]: " << size << " = " << size*8.0/(float)n << " bpc" << endl;
	cout << " LbRevF lenght " << nEF << endl;

	size = nDocs+1;
	this->dictD = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of dictD[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->dictDst = new bool[size];
	size *= sizeof(bool);
	sizeDS += size;
	cout << " ** size of dictDst[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyDoc = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of keyDoc[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->keyFq = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of keyFq[1..nDocs]: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nDocs+1;
	this->heap = new uint[size];
	size *= sizeof(uint);
	sizeDS += size;
	cout << " ** size of heap[1..nDocs]" << size << " = " << size*8.0/(float)this->n << " bpc" << endl;

	size = nod*lgD;
	if(size%W64)
		size = size/W64 + 1;
	else
		size /= W64;
	this->DocLZ = new ulong[size];
	size *= sizeof(ulong);
	sizeDS += size;
	cout << " ** size of DocLZ[1..n']: " << size << " = " << size*8.0/(float)this->n << " bpc" << endl;
	//cout << "====================================================" << endl << endl;

	// set initial open parenthesis in DFUDS sequences
	setBit64(PLz, 0);
	setBit64(PRev, 0);
	createStructure(dirSaveLoad);
}

void TopkLZhq::createStructure(char dirSaveLoad[300]){
	ulong i, j, pos, posRev, posF, posExtF, posFU;
	cout << "  Create LZ Top-k list Structure from RevTrie..." << endl;

	//[1] Create the DFUDS sequences of parentheses PRev, the bitstring FRev to mark the fictitious nodes, and the label vectors LbRev and LbRevF...
	fRev_b = bit_vector(nRev/2, 0);
	cout << "lenght of fRev_b (posF) = " << fRev_b.bit_size() << endl;
	fictU_b = bit_vector(nF, 0);
	cout << "lenght of fictU_b (posFU) = " << fictU_b.bit_size() << endl;
	listF_b = bit_vector(nEF, 0);
	cout << "lenght of listF_b (posExtF) = " << listF_b.bit_size() << endl;

	posF = posExtF = posFU = 0;
	genSeqExtraFictRevTrie(tries->revTrie, &posExtF, &posF, &posFU);
	cout << "posExtF " << posExtF << ", posF " << posF << ", posFU " << posFU << endl;
	cout << "posF = " << posF << " ? nRev/2 = " << nRev/2 << endl;
	cout << "genSeqExtraFictRevTrie OK" << endl;

	pos = 1; posRev = 0;
	genDFUDSSeqRevTrie(tries->revTrie, &pos, &posRev);
	cout << "pos " << pos << ", posRev (LbRev) " << posRev << endl;

	pos = 1;
	genDFUDSSeqLZTrieNotLbl(tries->lzTrie, &pos);
	cout << "pos (Plz) " << pos << " =? nLz = " << nLz <<endl;

	if(TEST){ // test for balanced sequences
		long int sum = 0;
		ulong r=0;

		for (; r<nLz; r++){
			if(readBit64(PLz, r))
				sum++;
			else sum--;
		}
		if(sum != 0){
			cout << " PLz DFUDS is not a balanced sequence of parentheses !! " << endl;
			exit(1);
		}
		//else cout << " Test for PLz. DFUDS is a well balanced sequence of parentheses !! " << endl;

		sum = r = 0;
		for (; r<nRev; r++){
			if(readBit64(PRev, r))
				sum++;
			else sum--;
		}
		if(sum != 0){
			cout << " PRev DFUDS is not a balanced sequence of parentheses !! " << endl;
			exit(1);
		}
		//else cout << " Test for PRev. DFUDS is a well balanced sequence of parentheses !! " << endl;
	}

	fRev_rrr = rrr_vector<127>(fRev_b);
	sizeDS += size_in_bytes(fRev_rrr);
	cout << " ** size of fRev_rrr: " << size_in_bytes(fRev_rrr) << " = " << size_in_bytes(fRev_rrr)*8.0/(float)n << " bpc" << endl;
	fRev_rank = rrr_vector<127>::rank_1_type(&fRev_rrr);

	fictU_rrr = rrr_vector<127>(fictU_b);
	sizeDS += size_in_bytes(fictU_rrr);
	cout << " ** size of fictU_rrr: " << size_in_bytes(fictU_rrr) << " = " << size_in_bytes(fictU_rrr)*8.0/(float)n << " bpc" << endl;
	fictU_rank = rrr_vector<127>::rank_1_type(&fictU_rrr);

	listF_rrr = rrr_vector<127>(listF_b);
	sizeDS += size_in_bytes(listF_rrr);
	cout << " ** size of listF_rrr: " << size_in_bytes(listF_rrr) << " = " << size_in_bytes(listF_rrr)*8.0/(float)n << " bpc" << endl;
	lisF_sel = rrr_vector<127>::select_1_type(&listF_rrr);

	if(TRACE){
		cout << "RevTrie with " << nod << " nodes, " << tries->nFictU << " nFictU nodes and " << tries->nExtraFict << " nExtraFict nodes" << endl;
		//tries->listRevTrie();

		cout << " DFUDFS PLz[0.." << nLz-1 << "]:" << endl;
		cout << readBit64(PLz, 0);
		for(i=1; i<nLz; i++){
			if(i%10 == 0)
				cout << "-";
			cout << readBit64(PLz, i);
		}
		cout << endl;

		cout << " DFUDFS PRev[0.." << nRev-1 << "]:" << endl;
		cout << readBit64(PRev, 0);
		for(i=1; i<nRev; i++){
			if(i%10 == 0)
				cout << "-";
			cout << readBit64(PRev, i);
		}
		cout << endl;
		cout << " fRev[0.." << fRev_b.bit_size()-1 << "]:" << endl;
		cout << fRev_b[0];
		for(i=1; i<fRev_b.bit_size(); i++){
			if(i%10 == 0)
				cout << "-";
			cout << fRev_b[i];
		}
		cout << endl;
		cout << " fictU[0.." << fictU_b.bit_size()-1 << "]:" << endl;
		cout << fictU_b[0];
		for(i=1; i<fictU_b.bit_size(); i++){
			if(i%10 == 0)
				cout << "-";
			cout << fictU_b[i];
		}
		cout << endl;
		cout << " LbRev[0.." << nodRev-1 << "]:" << endl;
		for(i=0; i<nodRev; i++)
			cout << LbRev[i];

		cout << endl;
		if (nEF){
			cout << " LbRevF[0.." << nEF-1 << "]:" << endl;
			for(i=0; i<nEF; i++)
				cout << LbRevF[i];
			cout << endl;
			cout << " listF[0.." << listF_b.bit_size()-1 << "]:" << endl;
			for(i=0; i<listF_b.bit_size(); i++)
				cout << listF_b[i];
			cout << endl;
		}else
			cout << " LbRevF does not have any extra labels (unary path of fictitious node with length > 1)" << endl;

		cout << " EndDocs[0.." << nDocs-1 << "]:" << endl;
		for(i=0; i<nDocs; i++)
			cout << EndDocs[i] << " ";
		cout << endl;

		cout << " Node array[0.." << nod-1 << "]:" << endl;
		cout << getNum64(Node, 0, lgNod) << " ";
		for(i=1, j=lgNod; i<nod; i++, j+=lgNod){
			if(i%10 == 0)
				cout << "- ";
			cout << getNum64(Node, j, lgNod) << " ";
		}
		cout << endl;
	}
	cout << "  Create the representation with Range min-max of PLz and PRev DFUDS..." << endl;
	treeLz = new RangeMMTree64(PLz, nLz, NULL, false);
	sizeDS += treeLz->sizeRMM;
	cout << " ** size of treeLz: " <<  treeLz->sizeRMM << " = " <<  treeLz->sizeRMM*8.0/(float)n << " bpc" << endl;

	treeRev = new RangeMMTree64(PRev, nRev, LbRev, false);
	sizeDS += treeRev->sizeRMM;
	cout << " ** size of treeRev: " <<  treeRev->sizeRMM << " = " <<  treeRev->sizeRMM*8.0/(float)n << " bpc" << endl;

	cout << "  Create DocLZ..." << endl;
	pos = 0;
	// making DocLZ in LZ78Tries
	genDocArray(tries->lzTrie, DocLZ, &pos);
	cout << "  Document array[0.." << nod-1 << "] and " << pos << " bits" << endl;
	if(TRACE){
		cout << "  Document array[0.." << nod-1 << "]:" << endl;
		for(i=j=0; i<nod; i++, j+=lgD){
			if(i%10 == 0) cout << "-";
			cout << getNum64(DocLZ, j, lgD);
		}
		cout << endl;
	}

	////////////////////////// TOP-K STRUCTURE //////////////////////////
	//cout << " Create the top-k list in each node..." << endl;
	bTopk_b = bit_vector(nRev/2, 0);
	cout << "lenght of bTopk_b = " << bTopk_b.bit_size() << endl;

	// marks the true nodes that store topk list
	nodeWithTopkList(tries->revTrie);

	pos = 0; // pos count only the n' true nodes
	marksTopkNodes(tries->revTrie, &pos, 0);

	// bTopk has nTopk 1's and z 0's, where n'=nTopk+z
	bTopk_rrr = rrr_vector<127>(bTopk_b);
	sizeDS += size_in_bytes(bTopk_rrr);
	cout << " ** size of bTopk_rrr: " << size_in_bytes(bTopk_rrr) << " = " << size_in_bytes(bTopk_rrr)*8.0/(float)n << " bpc" << endl;

	bTopk_rank = rrr_vector<127>::rank_1_type(&bTopk_rrr);
	nTopk = bTopk_rank.rank(bTopk_b.bit_size());
	cout << "   bTopk has " << nTopk << " nodes marked (" << ((float)nTopk/(float)nod)*100.0 << "% of total = " << bTopk_b.bit_size() << ")" << endl;
	bDL_b = bit_vector(nTopk, 0);
	sizeDS += size_in_bytes(bDL_b);
	cout << " ** size of bDL_b: " << size_in_bytes(bDL_b) << " = " << size_in_bytes(bDL_b)*8.0/(float)n << " bpc" << endl;

	ulong len;
	if(TRACE){
		len = bTopk_b.bit_size();
		cout << " bTopk[0.." << len-1 << "]:" << endl;
		cout << bTopk_b[0];
		for(i=1; i<len; i++){
			if(i%10 == 0)
				cout << "-";
			cout << bTopk_b[i];
		}
		cout << endl;
	}

	//////////////////////////////////////////////////////// docTopk /////////////////////////////////////////////////////////////
	uint* dNodK = new uint[nDocs];		// array of K most frequent id's
	uint* fNodK = new uint[nDocs];		// array of K frequencies
	uint* fNod = new uint[nDocs];

	ulong lenD, nDTot;
	pos = nDTot = lenD = 0;
	len = MAX_LEN_DOCS*n;					// Maximum length allowed is 24 bpc = 3|Text|
	if(len%W64)
		len = len/W64 + 1;
	else
		len /= W64;
	ulong *docTopkAux = new ulong[len];
	//cout << "largo(W) de docTopkAux = " << len << endl;
	countLengthFrequencies(fNod, tries->revTrie, &pos, &lenD, &nDTot, dNodK, fNodK, docTopkAux);
	cout << "   Number of documents that will be stored in top-k list: " << lenD << ", nodes marked: " << nDTot << endl;

	{	// save seq[] for testing and delete it
		char fileName[300];
		strcpy(fileName, "");
		strcpy(fileName, dirSaveLoad);
		strcat(fileName, "sequence.test");
		ofstream os (fileName, ios::binary);
		os.write((const char*)seq, n*sizeof(uchar));
		if(TRACE) cout << " .- Seq[] " << n*sizeof(uchar) << " Bytes" << endl;
		os.close();
		delete []seq;
	}

	auxLenD = lenD;
	if (nDTot != nTopk){
		cout << "ERROR. nDTot = " << nDTot << " != nTopk = " << nTopk << endl;
		exit(1);
	}
	delete [] dNodK;
	delete [] fNodK;
	delete [] fNod;

	// for documents in topK list...
	limDTk_b = bit_vector(lenD, 0);
	cout << "??? limDTk_b's length = " << limDTk_b.bit_size() << ", lgD = " << (uint)lgD << endl;
	len = lenD * lgD;
	if(len%W64)
		len = len/W64 + 1;
	else
		len /= W64;
	docTopk = new ulong[len];
	cout << "largo(W) de docTopk = " << len << endl;
	// to store the documents
	for(i=0; i<len; i++)
		docTopk[i] = docTopkAux[i];
	delete [] docTopkAux;

	len *= sizeof(ulong);
	sizeDS += len;
	cout << " ** size of docTopk: " << len << " = " << len*8.0/(float)n << " bpc" << endl;
	pos = lenD = 0;
	marksBoundaries(tries->revTrie, &pos, &lenD);
	cout << "     Number of documents stored in top-k list: " << lenD << " == " << limDTk_b.bit_size() << endl;

	limDTk_rrr = rrr_vector<127>(limDTk_b);
	limDTk_rank = rrr_vector<127>::rank_1_type(&limDTk_rrr);
	limDTk_sel = rrr_vector<127>::select_1_type(&limDTk_rrr);
	sizeDS += size_in_bytes(limDTk_rrr);

	cout << " ** size of limDTk_rrr: " << size_in_bytes(limDTk_rrr) << " = " << size_in_bytes(limDTk_rrr)*8.0/(float)n << " bpc" << endl;
	cout << " limDTk has " << limDTk_rank.rank(limDTk_b.bit_size()) << " 1's " << endl;
	if(TRACE){
		cout << " bDL[0.." << nTopk << "]:" << endl;
		cout << bDL_b[0];
		for(i=1; i<nTopk; i++){
			if(i%10 == 0)
				cout << "-";
			cout << bDL_b[i];
		}
		cout << endl;

		cout << " limDTk[0.." << limDTk_b.bit_size()-1 << "]:" << endl;
		cout << limDTk_b[0];
		for(i=1; i<limDTk_b.bit_size(); i++){
			if(i%10 == 0)
				cout << "-";
			cout << limDTk_b[i];
		}
		cout << endl;

		ulong pv, r, l, end;
		cout << "documents stored as Topk-List:" << endl;
		end = nRev/2;
		for(uint pre=0; pre<end; pre++){
			if (bTopk_b[pre]){
				pv = bTopk_rank.rank(pre);

				l = limDTk_sel.select(pv);		// start position of this block of id documents
				if(pv<nTopk)
					r = limDTk_sel.select(pv+1);
				else
					r = limDTk_b.bit_size();

				for(;l<r; l++){
					cout << getNum64(docTopk, l*lgD, lgD) << " ";
				}
				cout << endl;
			}
		}
	}
	cout << " PERFECT !!" << endl;
}

// sort the list in increasing mode by frequency
void TopkLZhq::sortTopkDocsByFq(uint* fNod, uint k, uint* dNodK, uint* fNodK){
	uint i, freq, len, items;
	int j;

	for (i=0; i<nDocs; i++){
		fNodK[i] = 0;
		dNodK[i] = 0;
	}
	for (i=0; i<nDocs && fNod[i]==0; i++);
	fNodK[0] = fNod[i];	// insert first document 'i' with freq>0
	dNodK[0] = i;
	items = 1;

	for (i++; i<nDocs; i++){
		freq = fNod[i];
		len = docLength[i];
		if(freq && items < k){
			// insertion sort for this document...
			for (j=items; j>0; j--){
				if (freq > fNodK[j-1]){
					fNodK[j] = fNodK[j-1];
					dNodK[j] = dNodK[j-1];
				}else{
					if (freq == fNodK[j-1] && len > docLength[dNodK[j-1]]){
						fNodK[j] = fNodK[j-1];
						dNodK[j] = dNodK[j-1];
					}else
						break;
				}
			}
			fNodK[j] = freq;
			dNodK[j] = i;
			items++;
		}else{
			if (freq && freq >= fNodK[k-1]){
				// insertion sort for this document...
				for (j=k-2; j>=0 && freq > fNodK[j]; j--){
					fNodK[j+1] = fNodK[j];
					dNodK[j+1] = dNodK[j];
				}
				for (; j>=0 && freq == fNodK[j] && len > docLength[dNodK[j]]; j--){
					fNodK[j+1] = fNodK[j];
					dNodK[j+1] = dNodK[j];
				}
				fNodK[j+1] = freq;
				dNodK[j+1] = i;
			}
		}
	}

	if (TRACE && false){
		cout << "sorted top-k list..." << endl;
		for (i=0; i<k; i++)
			cout << "doc: " << dNodK[i] << ", freq: " << fNodK[i] << endl;
	}
}

// set boundaries between each sublist of docId stored in docTopk
void TopkLZhq::marksBoundaries(LZ78Tries64::RevNode* node, ulong *pre, ulong *posD){
	LZ78Tries64::RevNode* p = node->fChildRev;

	if (bTopk_b[*pre]){
		limDTk_b[*posD] = 1;		// set initial document in this list
		(*posD) += node->kStar;
	}
	(*pre)++;

	for(uint i=0; i<node->nChildren; i++){
		marksBoundaries(p, pre, posD);
		p = p->nextSibRev;
	}
}

void TopkLZhq::create_sNode(LZ78Tries64::RevNode* nod, ulong *pre, ulong *aux_sNode, ulong *posNod){
	LZ78Tries64::RevNode* p = nod->fChildRev;
	ulong i;

	if (bTopk_b[*pre] == false && nod->fict == false){
		aux_sNode[*posNod] = nod->nodeLZ->preorder;		// set initial document in this list
		(*posNod)++;
	}
	(*pre)++;

	for(i=0; i<nod->nChildren; i++){
		create_sNode(p, pre, aux_sNode, posNod);
		p = p->nextSibRev;
	}
}

bool TopkLZhq::searchPatternInFMI(uchar* pat, uint m, uint *fqCSA){
	uint doc;
	string query = string((char *)pat);
	size_t occs = sdsl::count(fmi, query.begin(), query.begin()+m);
	auto locations = locate(fmi, query.begin(), query.begin()+m);
	//cout << "Total occurrences found with FMI : " << occs << endl;

	//cout << "locations..." << endl;
	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << endl;
		doc = searchDocPosition(locations[i]);
		if (doc >= nDocs){
			cout << "ERROR. doc = " << doc << " >= nDocs = " << nDocs << endl;
			exit(1);
		}
		(fqCSA[doc])++;

	}
	//cout << endl;


	if (occs == 0)
		return false;
	else{
		decltype(locations) empty;
		locations.swap(empty);
	}

	return true;
}

// it counts and stores in 'posBit' the total length frequencies (in bits)
void TopkLZhq::countLengthFrequencies(uint *fNod, LZ78Tries64::RevNode* node, ulong *pre, ulong *posD, ulong *nodTopk, uint* dNodK, uint* fNodK, ulong *docTopkAux){
	LZ78Tries64::RevNode* p = node->fChildRev;
	uint i;

	if (bTopk_b[*pre]){
		LZ78Tries64::RevNode* q = p;
		uint m, extraLb = 0;

		if(node->fict){
			extraLb = 1+node->lenUPath;
			LZ78Tries64::RevNode* fat = node;
			for(i=0; i<node->nChildren && q->fict; i++)
				q = q->nextSibRev;
			while(!q){
				fat = fat->fChildRev;
				q = fat->fChildRev;
				extraLb += 1+fat->lenUPath;
				for(i=0; i<fat->nChildren && q->fict; i++)
					q = q->nextSibRev;
			}
		}else
			q = node;
		ulong j = sep_sel.select(q->idNode)+extraLb;
		if((ulong)q->idNode < nod)
			m = sep_sel.select(q->idNode+1) - j;
		else
			m = n - j;
		for(i=0; i<m; i++)
			auxPath[i] = seq[i+j];
		auxPath[m+1] = '\0';

		for (i=0; i<nDocs; i++)
			fNod[i] = 0;
		if (searchPatternInFMI(auxPath, m, fNod) == false){
			cout << "ERROR, pattern not found in the FMI. pattern[0.." << m-1 << "] = [" << auxPath << "]" << endl;
			cout << "nod->idNode = " << node->idNode << endl;
			exit(1);
		}
		for (i=j=0; i<nDocs; i++){
			if (fNod[i])
				j++;
		}
		node->kStar = log(node->nOccT1/g) / log(2);
		node->kStar = pow(2, node->kStar );
		if(g*node->kStar  <= node->nOccT1)
			(node->kStar) <<= 1;
		if (node->kStar > nDocs)
			node->kStar = nDocs;
		if (j <= node->kStar){
			bDL_b[*nodTopk] = 1;
			node->kStar = j;
		}
		if((*posD+node->kStar)*lgD >= MAX_LEN_DOCS*n){
			cout << "ERROR.... docTopk requires more than " << MAX_LEN_DOCS << " bpc !!" << endl;
			exit(0);
		}
		sortTopkDocsByFq(fNod, node->kStar , dNodK, fNodK);
		for (i=0, j=(*posD)*lgD; i<node->kStar; i++, j+=lgD){
			//if (TRACE)	cout << dNodK[i] << " ";
			setNum64(docTopkAux, j, lgD, dNodK[i]);
		}

		(*posD) += node->kStar;
		(*nodTopk)++;
	}
	(*pre)++;

	for(i=0; i<node->nChildren; i++){
		countLengthFrequencies(fNod, p, pre, posD, nodTopk, dNodK, fNodK, docTopkAux);
		p = p->nextSibRev;
	}
}

// count the frequencies for all occurrence in this node (node for a true phrase or for a fictitious node)
ulong TopkLZhq::countMyFrequencies(LZ78Tries64::RevNode* nod, ulong* fNod){
	ulong nOcc = 0;
	if (nod->fict == false){
		ulong doc, rLz = tries->getLastLZLeaf(nod->nodeLZ);
		ulong lLz = nod->nodeLZ->preorder;
		nOcc = rLz-lLz+1;
		for(; lLz<=rLz; lLz++){
			doc = getNum64(DocLZ, lLz*lgD, lgD);
			(fNod[doc])++;
		}
	}

	LZ78Tries64::RevNode* p = nod->fChildRev;
	for(uint i=0; i<nod->nChildren; i++){
		nOcc += countMyFrequencies(p, fNod);
		p = p->nextSibRev;
	}

	return nOcc;
}

// determine the doc Id for the text position
uint TopkLZhq::searchDocPosition(ulong pos){
	uint m, ini=0, end=nDocs-1;

	// binary search in the interval endDocs[0..nPhra-1]
	while(ini < end){
		m = ini+(end-ini)/2;
		if (pos > EndDocsPos[m])
			ini = m+1;
		else{
			if(m){
				if (pos > EndDocsPos[m-1])
					return m;
				else
					end = m-1;
			}else
				return 0;
		}
	}

	return ini;
}

// count the frequencies for all occurrence in this node (node for a true phrase or for a fictitious node)
ulong TopkLZhq::countMyFrequencies(LZ78Tries64::RevNode* nod){
	ulong nOcc = 0;
	if (nod->fict == false)
		nOcc = tries->getLastLZLeaf(nod->nodeLZ)-nod->nodeLZ->preorder+1;

	LZ78Tries64::RevNode* p = nod->fChildRev;
	for(uint i=0; i<nod->nChildren; i++){
		nOcc += countMyFrequencies(p);
		p = p->nextSibRev;
	}

	return nOcc;
}

// count the length of its document segments that includes this node and also the number of different documents
ulong TopkLZhq::nodeWithTopkList(LZ78Tries64::RevNode *nod){
	nod->nOccT1 = 0;

	if (nod->idNode > 0)
		(nod->nOccT1) += tries->getLastLZLeaf(nod->nodeLZ) - nod->nodeLZ->preorder + 1;

	LZ78Tries64::RevNode* p = nod->fChildRev;
	for(uint i=0; i<nod->nChildren; i++){
		(nod->nOccT1) += nodeWithTopkList(p);
		p = p->nextSibRev;
	}
	return nod->nOccT1;
}

// marks the nodes that store a topk list in bTopk[]. Set the bits in bBlockLZ[], and filled the array fatherLZ[1..n']
void TopkLZhq::marksTopkNodes(LZ78Tries64::RevNode *node, ulong *pre, uint childFather){
	if(node->idNode && node->symbol != cutDoc){	// isn't it the root node ?
		if (node->fict==false || (node->nChildren > 1)){
			if (node->nOccT1 > 2*g)
				bTopk_b[*pre] = 1;
		}
	}
	(*pre)++;

	if (node->nChildren){
		LZ78Tries64::RevNode *child = node->fChildRev;
		childFather = node->nChildren;
		for(uint i=0; i<node->nChildren; i++){
			/*if(child == NULL){
				cout << "marksTopkNodes EEEEEEEEEEEEEEERRRRRRRRRRRRR child == NULL" << endl;
				cout << "i = " << i << ", node->nChildren = " << node->nChildren << ", id = " << node->idNode << endl;
				exit(1);
			}*/
			marksTopkNodes(child, pre, childFather);
			child = child->nextSibRev;
		}
	}
}

// determine the doc Id for the phrase idNode
uint TopkLZhq::searchDocument(ulong idNode){
	uint m, ini=0, end=nDocs-1;

	// binary search in the interval endDocs[0..nPhra-1]
	while(ini < end){
		m = ini+(end-ini)/2;
		if (idNode > EndDocs[m])
			ini = m+1;
		else{
			if(m){
				if (idNode > EndDocs[m-1])
					return m;
				else
					end = m-1;
			}else
				return 0;
		}
	}

	return ini;
}

// create the document array for LZTrie's phrases
void TopkLZhq::genDocArray(LZ78Tries64::LZNode* nod, ulong *DocLZ, ulong *ini){
	LZ78Tries64::LZNode* p = nod->fChildLZ;
	uint doc = searchDocument(nod->idNode);

	setNum64(DocLZ, *ini, lgD, doc);
	(*ini) += lgD;
	for(uint i=0; i<nod->nChildren; i++){
		genDocArray(p, DocLZ, ini);
		p = p->nextSibLZ;
	}
}

// set extra fictitious labels...
void TopkLZhq::genSeqExtraFictRevTrie(LZ78Tries64::RevNode *node, ulong *posExtF, ulong *pfRev, ulong *posFU){
	LZ78Tries64::RevNode *child;
	ulong i;

	if (node->fict){
		fRev_b[*pfRev] = 1;
		(*posFU)++;
		if(node->nChildren == 1){
			child = node->fChildRev;
			if (child->fict){
				nEFHead++;
				fictU_b[*posFU-1] = 1;
				listF_b[*posExtF] = 1;
				LbRevF[*posExtF] = child->symbol;
				auxPath[0] = child->symbol;
				uint unary = 1;
				(*posExtF)++;

				// search the last node in this unary path...
				while (child->nChildren == 1 && child->fChildRev->fict){
					LbRevF[*posExtF] = child->fChildRev->symbol;
					auxPath[unary] = child->fChildRev->symbol;
					unary++;
					(*posExtF)++;
					child = child->fChildRev;
				}
				node->nChildren = child->nChildren;
				node->fChildRev = child->fChildRev;
				node->uPath = true;
				node->lenUPath = unary;
				node->uSymbols = new uchar[unary];
				if (unary >= MAX_DEP){
					cout << "ERROR unary = " << unary << " >= MAX_DEP = " << MAX_DEP << endl;
					exit(1);
				}

				for (i=0; i<unary; i++)
					node->uSymbols[i] = auxPath[i];

			}
		}
	}
	(*pfRev)++;
	/*if (*pfRev > fRev_b.bit_size()){
		cout << "ERRRRR: *pfRev = " << *pfRev << " > " << "fRev_b.bit_size() = " << fRev_b.bit_size() << endl;
		exit(1);
	}*/

	// recursive call for all children
	if (node->nChildren){
		child = node->fChildRev;
		for(i=0; i<node->nChildren; i++){
			/*if(!child){
				cout << "ERRRRR: NO child and node->nChildren = " << node->nChildren << ", for node = " << node->idNode << endl;
				exit(1);
			}*/
			genSeqExtraFictRevTrie(child, posExtF, pfRev, posFU);
			child = child->nextSibRev;
		}
	}
}

//	create 	PRev, LbRev and LbRevF structures
void TopkLZhq::genDFUDSSeqRevTrie(LZ78Tries64::RevNode *node, ulong *pos, ulong *posRev){
	LZ78Tries64::RevNode *child;
	ulong i;

	// set open parentheses
	for(i=0, child = node->fChildRev; i<node->nChildren; i++, (*pos)++, child=child->nextSibRev, (*posRev)++){
		LbRev[*posRev] = child->symbol;
		setBit64(PRev, *pos);
		//cout << "(";
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}

	// set close parenthesis
	cleanBit64(PRev, *pos);
	(*pos)++;
	//cout << ")";

	// recursive call for all children
	child = node->fChildRev;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqRevTrie(child, pos, posRev);
		child = child->nextSibRev;
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}
}

//	create 	PLZ	:	LZTrie's DFUDS representation without Labels
void TopkLZhq::genDFUDSSeqLZTrieNotLbl(LZ78Tries64::LZNode *node, ulong *pos){
	LZ78Tries64::LZNode *child;
	ulong i;

	// set open parentheses
	for(i=0, child = node->fChildLZ; i<node->nChildren; i++, (*pos)++, child=child->nextSibLZ)
		setBit64(PLz, *pos);

	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}

	// set close parenthesis
	cleanBit64(PLz, *pos);
	(*pos)++;

	// recursive call for all children
	child = node->fChildLZ;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqLZTrieNotLbl(child, pos);
		child = child->nextSibLZ;
	}
	if (child){
		cout << "ERRR.. nod " << node->idNode << " has " << node->nChildren << " children, but there is also: " << child->idNode << endl;
		exit(1);
	}
}

// ELIMINAR
void TopkLZhq::readListFilesTrec(char *inputFile, ulong max_size){
	ulong i, len, lenText;
	ulong sizeFiles = 0;
	char fileName[300];

	std::ifstream in(inputFile);
	string line;
	std::getline(in,line);
	while(in && (sizeFiles < max_size)){
	    strcpy(fileName,line.c_str());
	    //cout << "File: " << fileName << endl;
	    std::ifstream input(fileName);
		assert(input.good());
		input.seekg(0, ios_base::end);
		len = (size_t)input.tellg();
		if(len > 1){
			n += len;
			sizeFiles += len;
			nDocs++;
		}
		input.close();
		std::getline(in,line);
	}
	in.close();
	docLength = new ulong[nDocs];
	sizeDS += nDocs*sizeof(ulong);
	cout << " ** size of docLength[1..nDocs] " << nDocs*sizeof(ulong) << " = " << (nDocs*sizeof(ulong))*8.0/(float)this->n << " bpc" << endl;

	cout << "Length of generalize text(n): " << n << ", in " << nDocs << " Documents" << endl;
	// allocate to memory for text...
	seq = new uchar[n+1];
	char *aux;
	std::ifstream in2(inputFile);
	lenText = 0;

	for(ulong texts=0; texts < nDocs;){
		std::getline(in2,line);
		strcpy(fileName,line.c_str());
		std::ifstream input(fileName); 			// open file
		//cout << "... File: " << fileName << endl;
		assert(input.good()); 				// check open
		input.seekg(0, ios_base::end);			// move to the end
		len = (size_t)input.tellg();			// add the final symbol (smallest)
		if(len > 1){
			docLength[texts] = len;
			aux = new char[len];
			//if (TRACE)
			//cout << "To read " << fileName << " pos: " << lenText << "..." << lenText+len-1 << endl;

			input.seekg(0, ios_base::beg);		// move back to the beginning
			if (input.good()){
				input.read(aux, len);

				len--;
				aux[len] = cutDoc;
				//cout << aux << endl;
				for (i=0; i<len; i++, lenText++){
					if((uchar)(aux[i]) <= cutDoc)
						seq[lenText] = ' ';
					else
						seq[lenText] = (uchar)(aux[i]);

					(charTf[seq[lenText]])++;
				}
				seq[lenText] = cutDoc;
				lenText++;
				assert(!input.fail());
				input.close();
			}else{
				cout << "Can't to open the file: <" << fileName << ">";
				exit(1);
			}
			delete [] aux;
			texts++;
		}
	}
	seq[n] = '\0';
	in2.close();
	char fileCpy[300];
	strcpy(fileCpy,inputFile);
	strcat(fileCpy, "_copy.txt");
	ofstream myfile;
	myfile.open (fileCpy);
	myfile << seq;
	myfile.close();

	if(TEST){
		uint DD = 0;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				DD++;
		}
		if(nDocs != DD){
			cout << "ERROR with nDocs in Sequence !! " << endl;
			cout << "nDocs = " << nDocs << " != " << DD << endl;
			exit(1);
		}
	}
	if(TRACE){		// listing original sequence
		cout << endl << "T[0.." << n-1 << "]:" << endl;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				cout << "$";
			else
				cout << seq[i];
		}
		cout << endl;
		cout << "pdocLength[0.." << nDocs-1 << "]: ";
		for(i=0; i<nDocs; i++)
			cout << docLength[i] << " ";
		cout << endl;
		cout << "charTf[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++)
			if (charTf[i])
				cout << i << "(" << charTf[i] << ") ";
		cout << endl;
	}
}

void TopkLZhq::readListFiles(char *inputFile, bool lowercase){
	ulong i, len, lenText;
	char fileName[300];

	std::ifstream in(inputFile);
	string line;
	std::getline(in,line);
	while(in){
	    strcpy(fileName,line.c_str());
	    //cout << "File: " << fileName << endl;
	    std::ifstream input(fileName);
		assert(input.good());
		input.seekg(0, ios_base::end);
		len = (size_t)input.tellg();
		if(len > 1){
			n += len;
			nDocs++;
		}
		input.close();
		std::getline(in,line);
	}
	in.close();
	docLength = new ulong[nDocs];
	sizeDS += nDocs*sizeof(ulong);
	cout << " ** size of docLength[1..nDocs] " << nDocs*sizeof(ulong) << " = " << (nDocs*sizeof(ulong))*8.0/(float)this->n << " bpc" << endl;

	cout << "Length of generalize text(n): " << n << ", in " << nDocs << " Documents" << endl;
	// allocate to memory for text...
	seq = new uchar[n];
	char *aux;
	std::ifstream in2(inputFile);
	lenText = 0;

	bool *present = new bool[LZ78Tries64::SIGMA];
	for (i=0; i<LZ78Tries64::SIGMA; i++)
		present[i] = false;
	uint sigma = 0;

	for(ulong texts=0; texts < nDocs;){
		std::getline(in2,line);
		strcpy(fileName,line.c_str());
		std::ifstream input(fileName); 			// open file
		//cout << "... File: " << fileName << endl;
		assert(input.good()); 				// check open
		input.seekg(0, ios_base::end);			// move to the end
		len = (size_t)input.tellg();			// add the final symbol (smallest)
		if(len > 1){
			docLength[texts] = len;
			aux = new char[len];
			//if (TRACE)
			//cout << "To read " << fileName << " pos: " << lenText << "..." << lenText+len-1 << endl;

			input.seekg(0, ios_base::beg);		// move back to the beginning
			if (input.good()){
				input.read(aux, len);

				len--;
				aux[len] = cutDoc;
				//cout << aux << endl;
				for (i=0; i<len; i++, lenText++){

					if (present[(uint)aux[i]] == false){
						sigma++;
						present[(uint)aux[i]] = true;
					}
					if((uchar)(aux[i]) <= cutDoc)
						seq[lenText] = ' ';
					else{
						if (lowercase)
							seq[lenText] = (uchar)(tolower(aux[i]));
						else
							seq[lenText] = (uchar)(aux[i]);
					}
					(charTf[seq[lenText]])++;
				}
				seq[lenText] = cutDoc;
				lenText++;
				assert(!input.fail());
				input.close();
			}else{
				cout << "Can't to open the file: <" << fileName << ">";
				exit(1);
			}
			delete [] aux;
			texts++;
		}
	}

	cout << " ## sigma original = " << sigma << endl;

	seq[n-1] = '\0';
	in2.close();
	char fileCpy[300];
	strcpy(fileCpy,inputFile);
	strcat(fileCpy, "_copy.txt");
	ofstream myfile;
	myfile.open (fileCpy);
	myfile << seq;
	myfile.close();
	seq[n-1] = cutDoc;

	if(TEST){
		uint DD = 0;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				DD++;
		}
		if(nDocs != DD){
			cout << "ERROR with nDocs in Sequence !! " << endl;
			cout << "nDocs = " << nDocs << " != " << DD << endl;
			exit(1);
		}
	}
	if(TRACE){		// listing original sequence
		cout << endl << "T[0.." << n-1 << "]:" << endl;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				cout << "$";
			else
				cout << seq[i];
		}
		cout << endl;
		cout << "pdocLength[0.." << nDocs-1 << "]: ";
		for(i=0; i<nDocs; i++)
			cout << docLength[i] << " ";
		cout << endl;
		cout << "charTf[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++)
			if (charTf[i])
				cout << i << "(" << charTf[i] << ") ";
		cout << endl;
	}
}

void TopkLZhq::readUniqueFile(char *inputFile, bool lowercase){
	ulong i, j, len;
	n = nDocs = 0;
	ifstream input(inputFile);			// open file
	assert(input.good()); 				// check open
	input.seekg(0, ios_base::end);		// move to the end
	n = (size_t)input.tellg();
	seq = new uchar[n];
	char *aux = new char[n];

	input.seekg(0, ios_base::beg);		// move back to the beginning
	if (input.good()){
		input.read(aux, n-1);
		for (i=0; i<n-1; i++){
			if((uchar)(aux[i]) == boundSymbol){
				nDocs++;
				while((uchar)(aux[i+1]) == boundSymbol){
					seq[i] = '\n';
					(charTf['\n'])++;
					i++;
				}
				seq[i] = cutDoc;
			}else{
				if((uchar)(aux[i]) <= cutDoc)
					seq[i] = '\n';
				else{
					if (lowercase)
						seq[i] = (uchar)(tolower(aux[i]));
					else
						seq[i] = (uchar)(aux[i]);
				}
			}
			(charTf[seq[i]])++;
		}
		if(seq[n-2] == cutDoc){
			seq[n-2] = '\n';
			(charTf['\n'])++;
			nDocs--;
		}
		seq[n-1] = cutDoc;
		nDocs++;
		//cout << seq << endl;
		assert(!input.fail());
		input.close();
	}else{
		cout << "Can't to open the file: <" << inputFile << ">";
		exit(1);
	}
	input.close();
	delete [] aux;

	docLength = new ulong[nDocs];
	sizeDS += nDocs*sizeof(ulong);
	cout << " ** size of docLength[1..nDocs] " << nDocs*sizeof(ulong) << " = " << (nDocs*sizeof(ulong))*8.0/(float)this->n << " bpc" << endl;
	cout << "Length of generalize text(n): " << n << ", in " << nDocs << " Documents" << endl;

	for (i=j=len=0; i<n; i++){
		if(seq[i] == cutDoc){
			docLength[j] = len;
			len=0;
			j++;
		}else
			len++;
	}
	if (j != nDocs){
		cout << "Error cutDocs Symbols = " << j << " != nDocs = " << nDocs << endl;
		exit(1);
	}

	seq[n-1] = '\0';
	char fileCpy[300];
	strcpy(fileCpy,inputFile);
	strcat(fileCpy, "_copy.txt");
	ofstream myfile;
	myfile.open (fileCpy);
	myfile << seq;
	myfile.close();
	seq[n-1] = cutDoc;

	if(TEST){
		uint DD = 0;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				DD++;
		}
		if(nDocs != DD){
			cout << "ERROR with nDocs in Sequence !! " << endl;
			cout << "nDocs = " << nDocs << " != " << DD << endl;
			exit(1);
		}
	}

	if(TRACE){		// listing original sequence
		cout << "T[0.." << n-1 << "]:" << endl;
		for(i=0; i<n; i++){
			if (seq[i] == cutDoc)
				cout << "$";
			else
				cout << seq[i];
		}
		cout << endl;
		cout << "pdocLength[0.." << nDocs-1 << "]: ";
		for(i=0; i<nDocs; i++)
			cout << docLength[i] << " ";
		cout << endl;
		cout << "charTf[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++)
			if (charTf[i])
				cout << i << "(" << charTf[i] << ") ";
		cout << endl;
	}
}

// make priorities for hash table by probability
void TopkLZhq::makePriorities(ulong* cPosTb){
	uint i, j, sigma = 0;
	uchar *prior = new uchar[LZ78Tries64::SIGMA];
	prior[0] = 0; // --> charTf[0];

	bool tr = TRACE;
	TRACE = false;

	for (i=1; i<LZ78Tries64::SIGMA; i++){
		if (charTf[i] > 0)
			sigma++;
		for (j=i; j>0 && charTf[prior[j-1]] < charTf[i]; j--){
			prior[j] = prior[j-1];
		}
		prior[j] = i; // --> charTf[i];
	}

	cout << " ## sigma = " << sigma << endl;

	if(TRACE){
		cout << "prior[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++){
			cout << i << "(" << (uint)prior[i] << ") ";
		}
		cout << endl;
	}

	uint rows = LZ78Tries64::SIGMA / LZ78Tries64::LENHASH;
	if (LZ78Tries64::SIGMA % LZ78Tries64::LENHASH)
		rows++;
	uint code = 0;
	for(uint fil=0; fil<rows; fil++){
		for(i=0; i<LZ78Tries64::LENHASH && code<LZ78Tries64::SIGMA; i++){
			if(fil%2 == 0)
				cPosTb[prior[code]] = i;
			else
				cPosTb[prior[code]] = LZ78Tries64::LENHASH-1-i;
			code++;
		}
	}

	if(TRACE){
		cout << "charPosTb[0..255]: ";
		for(i=0; i<LZ78Tries64::SIGMA; i++){
			cout << i << "(" << cPosTb[i] << ") ";
		}
		cout << endl;
	}
	delete [] prior;
	TRACE = tr;
}

// return true if the pattern is in the RevTrie and in '*pv' the respective preorder value in the tree,
// in '*x' stores its DFUDS bit position. Return false if the pattern does not exist
bool TopkLZhq::searchPattern_Rev(uchar* pat, uint m, ulong *x, ulong *pv){
	ulong d, w1, w2, sig0, pre=0, pos=*x;

	while (m > 0){
		m--;
		if (fRev_rrr[pre]){						// is it a fictitious node ?
			w1 = fRev_rank.rank(pre);

			if (fictU_rrr[w1]){					// is it a fictitious node with unary path of other fictitious nodes ?
				w2 = fictU_rank.rank(w1+1);
				ulong r,l=lisF_sel.select(w2);
				if (w2 < nEFHead)
					r = lisF_sel.select(w2+1);
				else
					r = nEF;
				if(m){
					while (l<r && LbRevF[l]==pat[m]){
						m--; l++;
						if (m==0) break;
					}
					if (l<r){
						if (m==0){
							if (LbRevF[l]!=pat[0])return false;
						}else return false;
					}else{
						pre++;
						if (fRev_rrr[pre]){
							if (fictU_rrr[w1+1]){
								sig0 = treeRev->selectNext0(pos);
								if (sig0 > n) sig0 = n;
								else if (pos>=sig0) return false;

								if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
									if(w2 > 1){
										d=1;
										treeRev->fwd_search(sig0-w2+1, &d, &pos);
										pos++;
									}else
										pos = sig0+1;
								}else return false;

								pre = pos-treeRev->rank_1(pos-1);
							}else m++;
						}else m++;
					}
				}else{
					if (LbRevF[l]==pat[0]){
						*x = pos;
						*pv = pre;
						return true;
					}
					else return false;
				}
			}else{									// normal case
				sig0 = treeRev->selectNext0(pos);
				if (sig0 > nRev) sig0 = nRev;
				else if (pos>=sig0) return false;

				if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
					if(w2 > 1){
						d=1;
						treeRev->fwd_search(sig0-w2+1, &d, &pos);
						pos++;
					}else
						pos = sig0+1;
				}else return false;
				pre = pos-treeRev->rank_1(pos-1);
			}
		}else{										// normal case
			sig0 = treeRev->selectNext0(pos);
			if (sig0 > nRev) sig0 = nRev;
			else if (pos>=sig0) return false;

			if (treeRev->thereIsChild(pos, pat[m], &w2, sig0-pos)){
				if(w2 > 1){
					d=1;
					treeRev->fwd_search(sig0-w2+1, &d, &pos);
					pos++;
				}else
					pos = sig0+1;
			}else return false;
			pre = pos-treeRev->rank_1(pos-1);
		}
	}
	*x = pos;
	*pv = pre;
	return true;
}

// It stores in docList[1..k] and frqList[1..k] the list and frequencies of the top-k most frequent documents for the pattern p[1..m]
ulong TopkLZhq::topKDocument(uchar* pat, uint m, uint* docList, uint k){
	ulong pv, x=1;

	if (searchPattern_Rev(pat, m, &x, &pv)){
		ulong pos, l, r;
		uint st=0;
		bool uniq, stored;
		uniq = stored = false;

		//cout << "pv: " << pv << endl;
		if (bTopk_rrr[pv]){
			stored = true;
			//cout << "Stored" << endl;
		}else{
			if(fRev_rrr[pv]){
				if(treeRev->hasUniqueChild(x) && bTopk_rrr[pv+1]){
					//cout << "Stored as unique Child, pv = " << pv << endl;
					pv++;
					stored = uniq = true;
				}
			}
		}

		if (stored){
			// Answer stored in the node pv, from lowest to highest frequency
			ulong ra = bTopk_rank.rank(pv+1);
			l = limDTk_sel.select(ra);		// number of this block of id documents
			if(ra<nTopk)
				r = limDTk_sel.select(ra+1);
			else
				r = limDTk_b.bit_size();
			st=r-l;						// number of different id's
			if (st > k)
				st = k;

			//cout << "Stored with " << st << " different id's" << endl;
			//cout << "bDL->getBit(ra-1) " << bDL_b[ra-1] << endl;
			//cout << "Documents, from l*lgD, l = " << l << endl;
			l *= lgD;
			for(pos=0; pos<st; l+=lgD, pos++){
				docList[pos] = getNum64(docTopk, l, lgD);
				//cout << docList[pos] << " ";
			}
			//cout << docList[0] << endl;
			if ((st == k) || bDL_b[ra-1])
				return st;

			// NOW, THERE ARE st ID_DOCS STORE IN docList[0..st-1]
			for (pos=0; pos<st; pos++)
				dictDst[docList[pos]] = true;
		}

		// Computes frequencies by brute force
		uint key=1;
		ulong i, len, fictBef, doc;
		len = treeRev->subTreeSize(x)-uniq;
		fictBef = fRev_rank.rank(pv);
		len -= fRev_rank.rank(pv+len) - fictBef;
		pv -= fictBef;

		if(st){
			for(i=pv*lgNod; len > 0; len--, i+=lgNod){
				l = getNum64(Node, i, lgNod);
				x = treeLz->select_0(l)+1;
				r = l + treeLz->subTreeSize(x);
				// scan DocLZ[l..r]. A new document doc has a key = dictD[doc], and his frequencies is in fqKey[key]
				//cout << "Docs for l(" << l << ") :  ";

				for (pos=l*lgD; l<r; l++, pos+=lgD){
					doc = getNum64(DocLZ, pos, lgD);
					//cout << doc << " ";
					if(dictDst[doc]) continue;
					x = dictD[doc];	// key
					if(x)
						(keyFq[x])++;
					else{
						// x==0 then it is a new ID
						dictD[doc] = key;
						keyDoc[key] = doc;
						keyFq[key] = 1;
						//cout << doc << "[" << key << "](" << keyFq[key] << ") ";
						key++;
					}
				}
			}
		}else{
			for(i=pv*lgNod; len > 0; len--, i+=lgNod){
				l = getNum64(Node, i, lgNod);
				x = treeLz->select_0(l)+1;
				r = l + treeLz->subTreeSize(x);
				// scan DocLZ[l..r]. A new document doc has a key = dictD[doc], and his frequencies is in fqKey[key]
				//cout << "Docs for l(" << l << ") :  ";

				for (pos=l*lgD; l<r; l++, pos+=lgD){
					doc = getNum64(DocLZ, pos, lgD);
					//cout << doc << " ";
					x = dictD[doc];	// key
					if(x)
						(keyFq[x])++;
					else{
						// x==0 then it is a new ID
						dictD[doc] = key;
						keyDoc[key] = doc;
						keyFq[key] = 1;
						//cout << doc << "[" << key << "](" << keyFq[key] << ") ";
						key++;
					}
				}
			}
		}

		//cout << endl;
		//cout << "key " << key << endl;
		//cout << "Freq [1...]" << endl;
		//for (l=1; l<key; l++){
		//	cout << l << "(" << keyFq[l] << ") ";
		//}
		//cout << endl;

		for (l=1; l<key; l++){
			heap[l] = l;
			swimInMaxQ(l);
			//for (uint a=1; a<=l; a++)
			//	cout << " " << keyFq[heap[a]];
			//cout << endl;
		}
		//for (l=1; l<key; l++)
		//	cout << " f" << (uint)keyFq[heap[l]];
		//cout << endl;

		key--;
		uint kd = st+key;
		if (k < kd)
			kd = k;
		k = kd-st;
		for (l=0; l<k; l++){
			x = heap[1];
			setTopMaxQ(key-l);
			//cout << " F" << (uint)keyFq[x] << endl;
			docList[l+st] = keyDoc[x];
		}

		// CLEAN ... to unmark documents
		for (l=1; l<=key; l++){
			dictD[keyDoc[l]] = 0;
			keyDoc[l] = keyFq[l] = 0;
		}
		for (pos=0; pos<st; pos++)
			dictDst[docList[pos]] = false;

		return kd;
	}
	return 0;
}

//=================================================================================================
// the doc with key dictD[k] go up until its correct position in the maximum queue
void TopkLZhq::swimInMaxQ(uint k){
	uint sw=heap[k];
	ulong doc=keyDoc[k];
	float num=keyFq[sw];

	while(k > 1 && keyFq[heap[k>>1]] < num){
		heap[k] = heap[k>>1];
		k >>= 1;
	}
	while(k > 1 && keyFq[heap[k>>1]] == num && docLength[keyDoc[k>>1]] < docLength[doc]){
		heap[k] = heap[k>>1];
		k >>= 1;
	}
	heap[k] = sw;
}


// set the new item at top of the queue heap[1..k]
void TopkLZhq::setTopMaxQ(uint k){
	uint j, m=2, i=1;
	ulong doc=keyDoc[k];

	float num=keyFq[heap[k]];

	while(m<k){
		j=m+1;
		if (j < k){
			if (keyFq[heap[m]] < keyFq[heap[j]])
				m=j;
			else{
				if (keyFq[heap[m]] == keyFq[heap[j]] && docLength[keyDoc[m]] < docLength[keyDoc[j]])
					m=j;
			}
		}
		if(num <= keyFq[heap[m]]){
			if(num < keyFq[heap[m]] || docLength[doc] < docLength[keyDoc[m]]){
				heap[i] = heap[m];
				i = m;
				m <<= 1;
			}else
				break;
		}else
			break;
	}
	heap[i] = heap[k];
}


// save the Data Structure in file 'fileName'
void TopkLZhq::saveDS(bool showSize){
	char *fileName = new char[300];
	cout << "Save data structure in folder " << dirStore << endl;

	strcpy(fileName, dirStore);
	strcat(fileName, "dataStructures.tk");
	ofstream os (fileName, ios::binary);
	cout << "   Saving. Data structure size: " << sizeDS << endl;

	os.write((const char*)&n, sizeof(ulong));
	os.write((const char*)&nDocs, sizeof(uint));
	os.write((const char*)&lgD, sizeof(ulong));
	os.write((const char*)&g, sizeof(uint));
	os.write((const char*)&cutDoc, sizeof(char));
	os.write((const char*)&nod, sizeof(ulong));
	os.write((const char*)&lgNod, sizeof(uint));
	os.write((const char*)&nF, sizeof(ulong));
	os.write((const char*)&nEF, sizeof(ulong));
	os.write((const char*)&nEFHead, sizeof(uint));
	os.write((const char*)&nodRev, sizeof(ulong));
	os.write((const char*)&nLz, sizeof(ulong));
	os.write((const char*)&nRev, sizeof(ulong));
	os.write((const char*)&nTopk, sizeof(ulong));

	ulong sizeDSav = 9*sizeof(ulong) + 4*sizeof(uint) + sizeof(char);	// size for variables

	os.write((const char*)docLength, nDocs*sizeof(ulong));				// save docLength[]
	sizeDSav += nDocs*sizeof(ulong);
	if(showSize) cout << " .- docLength[] " << nDocs*sizeof(ulong) << " Bytes" << endl;

	os.write((const char*)LbRev, nodRev*sizeof(uchar));					// save LbRevF[]
	sizeDSav += nodRev*sizeof(uchar);
	if(showSize) cout << " .- LbRev[] " << nodRev*sizeof(uchar) << " Bytes" << endl;

	os.write((const char*)LbRevF, nEF*sizeof(uchar));					// save LbRevF[]
	sizeDSav += nEF*sizeof(uchar);
	if(showSize) cout << " .- LbRevF[] " << nEF*sizeof(uchar) << " Bytes" << endl;

	ulong size = nod*lgD/W64;
	if ((nod*lgD)%W64)
		size++;
	os.write((const char*)DocLZ, size*sizeof(ulong));				// save DocLZ[]
	sizeDSav += size*sizeof(ulong);
	if(showSize) cout << " .- DocLZ[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = nod*lgNod/W64;
	if ((nod*lgNod)%W64)
		size++;
	os.write((const char*)Node, size*sizeof(ulong));				// save Node[]
	sizeDSav += size*sizeof(ulong);
	if(showSize) cout << " .- Node[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = limDTk_b.bit_size()*lgD/W64;
	if ((limDTk_b.bit_size()*lgD)%W64)
		size++;
	os.write((const char*)docTopk, size*sizeof(ulong));				// save docTopk[]
	sizeDSav += size*sizeof(ulong);
	if(showSize) cout << " .- docTopk[] " << size*sizeof(ulong) << " Bytes" << endl;
	os.close();

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rrr.tk");
	store_to_file(fRev_rrr, fileName);
	sizeDSav += size_in_bytes(fRev_rrr);
	if(showSize) cout << " .- fRev_rrr " << size_in_bytes(fRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rank.tk");
	store_to_file(fRev_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rrr.tk");
	store_to_file(fictU_rrr, fileName);
	sizeDSav += size_in_bytes(fictU_rrr);
	if(showSize) cout << " .- fictU_rrr " << size_in_bytes(fictU_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rank.tk");
	store_to_file(fictU_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "listF_rrr.tk");
	store_to_file(listF_rrr, fileName);
	sizeDSav += size_in_bytes(listF_rrr);
	if(showSize) cout << " .- listF_rrr " << size_in_bytes(listF_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "lisF_sel.tk");
	store_to_file(lisF_sel, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rrr.tk");
	store_to_file(bTopk_rrr, fileName);
	sizeDSav += size_in_bytes(bTopk_rrr);
	if(showSize) cout << " .- bTopk_rrr " << size_in_bytes(bTopk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rank.tk");
	store_to_file(bTopk_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bDL_b.tk");
	store_to_file(bDL_b, fileName);
	sizeDSav += size_in_bytes(bDL_b);
	if(showSize) cout << " .- bDL_b " << size_in_bytes(bDL_b) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rrr.tk");
	store_to_file(limDTk_rrr, fileName);
	sizeDSav += size_in_bytes(limDTk_rrr);
	if(showSize) cout << " .- limDTk_rrr " << size_in_bytes(limDTk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rank.tk");
	store_to_file(limDTk_rank, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_sel.tk");
	store_to_file(limDTk_sel, fileName);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LzTrie.tk");
	cout << "***   Saving LzTrie... " << endl;
	sizeDSav += treeLz->saveDS(fileName, showSize);
	if(showSize) cout << " .- treeLz " << treeLz->sizeRMM << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "RevTrie.tk");
	cout << "***   Saving RevTrie... " << endl;

	sizeDSav += treeRev->saveDS(fileName, showSize);
	if(showSize) cout << " .- treeRev " << treeRev->sizeRMM << " Bytes" << endl;

	cout << "   Total bytes saved from data structure: " << sizeDSav << endl;
}

// load the Data Structure from the file 'fileName'
void TopkLZhq::loadDS(bool showSize){
	cout << " Load data structure from " << dirStore << endl;
	char *fileName = new char[300];
	sizeDS = 0;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rrr.tk");
	load_from_file(fRev_rrr, fileName);
	sizeDS += size_in_bytes(fRev_rrr);
	if(showSize) cout << " ** size of fRev_rrr " << size_in_bytes(fRev_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fRev_rank.tk");
	load_from_file(fRev_rank, fileName);
	util::init_support(fRev_rank, &fRev_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rrr.tk");
	load_from_file(fictU_rrr, fileName);
	sizeDS += size_in_bytes(fictU_rrr);
	if(showSize) cout << " ** size of fictU_rrr " << size_in_bytes(fictU_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "fictU_rank.tk");
	load_from_file(fictU_rank, fileName);
	util::init_support(fictU_rank, &fictU_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "listF_rrr.tk");
	load_from_file(listF_rrr, fileName);
	sizeDS += size_in_bytes(listF_rrr);
	if(showSize) cout << " ** size of listF_rrr " << size_in_bytes(listF_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "lisF_sel.tk");
	load_from_file(lisF_sel, fileName);
	util::init_support(lisF_sel, &listF_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rrr.tk");
	load_from_file(bTopk_rrr, fileName);
	sizeDS += size_in_bytes(bTopk_rrr);
	if(showSize) cout << " ** size of bTopk_rrr " << size_in_bytes(bTopk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bTopk_rank.tk");
	load_from_file(bTopk_rank, fileName);
	util::init_support(bTopk_rank, &bTopk_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "bDL_b.tk");
	load_from_file(bDL_b, fileName);
	sizeDS += size_in_bytes(bDL_b);
	if(showSize) cout << " ** size of bDL_b " << size_in_bytes(bDL_b) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rrr.tk");
	load_from_file(limDTk_rrr, fileName);
	sizeDS += size_in_bytes(limDTk_rrr);
	if(showSize) cout << " ** size of limDTk_rrr " << size_in_bytes(limDTk_rrr) << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_rank.tk");
	load_from_file(limDTk_rank, fileName);
	util::init_support(limDTk_rank, &limDTk_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "limDTk_sel.tk");
	load_from_file(limDTk_sel, fileName);
	util::init_support(limDTk_sel, &limDTk_rrr);

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "dataStructures.tk");
	ifstream is(fileName, ios::binary);

	is.read((char*)&n, sizeof(ulong));
	is.read((char*)&nDocs, sizeof(uint));
	is.read((char*)&lgD, sizeof(ulong));
	is.read((char*)&g, sizeof(uint));
	is.read((char*)&cutDoc, sizeof(char));
	is.read((char*)&nod, sizeof(ulong));
	is.read((char*)&lgNod, sizeof(uint));
	is.read((char*)&nF, sizeof(ulong));
	is.read((char*)&nEF, sizeof(ulong));
	is.read((char*)&nEFHead, sizeof(uint));
	is.read((char*)&nodRev, sizeof(ulong));
	is.read((char*)&nLz, sizeof(ulong));
	is.read((char*)&nRev, sizeof(ulong));
	is.read((char*)&nTopk, sizeof(ulong));

	sizeDS += 9*sizeof(ulong) + 4*sizeof(uint) + sizeof(char);	// size for variables

	docLength = new ulong[nDocs];
	is.read((char*)docLength, nDocs*sizeof(ulong));
	sizeDS += nDocs*sizeof(ulong);
	if(showSize) cout << " ** size of docLength[] " << nDocs*sizeof(ulong) << " Bytes" << endl;

	LbRev = new uchar[nodRev];
	is.read((char*)LbRev, nodRev*sizeof(uchar));
	sizeDS += nodRev*sizeof(uchar);
	if(showSize) cout << " ** size of LbRev[] " << nodRev*sizeof(uchar) << " Bytes" << endl;

	LbRevF = new uchar[nEF];
	is.read((char*)LbRevF, nEF*sizeof(uchar));
	sizeDS += nEF*sizeof(uchar);
	if(showSize) cout << " ** size of LbRevF[] " << nEF*sizeof(uchar) << " Bytes" << endl;

	ulong size = nod*lgD/W64;
	if ((nod*lgD)%W64)
		size++;
	DocLZ = new ulong[size];
	is.read((char*)DocLZ, size*sizeof(ulong));
	sizeDS += size*sizeof(ulong);
	if(showSize) cout << " ** size of DocLZ[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = nod*lgNod/W64;
	if ((nod*lgNod)%W64)
		size++;
	Node = new ulong[size];
	is.read((char*)Node, size*sizeof(ulong));
	sizeDS += size*sizeof(ulong);
	if(showSize) cout << " ** size of Node[] " << size*sizeof(ulong) << " Bytes" << endl;

	size = limDTk_rrr.size()*lgD/W64;
	if ((limDTk_rrr.size()*lgD)%W64)
		size++;
	docTopk = new ulong[size];
	is.read((char*)docTopk, size*sizeof(ulong));
	sizeDS += size*sizeof(ulong);
	if(showSize) cout << " ** size of docTopk[] " << size*sizeof(ulong) << " Bytes" << endl;
	is.close();

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "LzTrie.tk");
	treeLz = new RangeMMTree64(fileName, false);
	sizeDS += treeLz->sizeRMM;
	if(showSize) cout << " ** size of treeLz " << treeLz->sizeRMM << " Bytes" << endl;

	strcpy(fileName, "");
	strcpy(fileName, dirStore);
	strcat(fileName, "RevTrie.tk");
	treeRev = new RangeMMTree64(fileName, false);
	treeRev->labels = LbRev;
	sizeDS += treeRev->sizeRMM;
	if(showSize) cout << " ** size of treeRev " << treeRev->sizeRMM << " Bytes" << endl;

	cout << " Data Structure loaded !!" << endl;
}

TopkLZhq::~TopkLZhq() {
	sizeDS = n = nDocs = lgD = g = cutDoc = nod = lgNod = 0;
	nF = nEF = nEFHead = nodRev = nLz = nRev = nTopk = 0;

	delete []docLength;
	delete []EndDocs;
	delete []charTf;
	delete []PLz;
	delete []PRev;
	delete []LbRev;
	delete []LbRevF;
	delete []DocLZ;
	delete []Node;
	delete []docTopk;
	delete []dictD;
	delete []keyDoc;
	delete []keyFq;
	delete []heap;

	tries->~LZ78Tries64();
	cout << "destroyeing treeLZ..." << endl;
	treeLz->~RangeMMTree64();
	cout << "destroyeing treeRev..." << endl;
	treeRev->~RangeMMTree64();
}

