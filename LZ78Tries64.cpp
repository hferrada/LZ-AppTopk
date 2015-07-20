/*
 * LZ78Tries64.cpp
 *
 *  Created on: 01-07-2014
 *      Author: hector
 */

#include "LZ78Tries64.h"
bool LZ78Tries64::TRACE = false;
uint LZ78Tries64::SIGMA = 256;
uint LZ78Tries64::LENHASH = 16;

LZ78Tries64::LZ78Tries64(uint nDocs, uchar endSymbolDoc, ulong* symbolHTable) {
	lzTrie = new LZNode();
	revTrie = new RevNode();

	lzTrie->idNode = 0;
	lzTrie->nChildren = 0;
	lzTrie->preorder = 0;
	lzTrie->nextSibLZ = lzTrie->fChildLZ = NULL;

	revTrie->idNode = 0;
	revTrie->nChildren = 0;
	revTrie->preorder = 0;
	revTrie->fict = 0;
	revTrie->uSymbols = NULL;
	revTrie->nodeLZ = lzTrie;
	revTrie->nextSibRev = revTrie->fChildRev = NULL;

	this->cPosTb = symbolHTable;
	this->cutDoc = endSymbolDoc;
	this->nPhra = 0;
	this->nFictU = 0;
	this->nExtraFict = 0;
	this->nDocs = nDocs;
	this->bad = -1;
	this->lgD = ceilingLog64(nDocs, 2);
}

void LZ78Tries64::createTriesFromText(uchar* text, ulong n, ulong** nodeArr, ulong *endDocs, ulong *endDocsPos, bit_vector *separator_b){
	ulong posIni, textPos;
	uint nDoc;
	LZNode* nodLZ;

	nDoc = textPos = 0;
	while (textPos < n){
		nPhra++;
		(*separator_b)[textPos] = 1;

		if (text[textPos] == cutDoc){
			nPhra--;
			// we do not insert phrases with a unique symbol = cutDoc
			endDocs[nDoc] = nPhra;
			endDocsPos[nDoc] = textPos;
			(*separator_b)[textPos] = 0;
			if (nDoc % 10000 == 0) cout << "nDoc = " << nDoc << ", textPos = " << textPos << ", nPhra = " << nPhra << endl;
			nDoc++;
			textPos++;
			continue;
		}
		posIni = textPos;
		nodLZ = new LZNode();
		nodLZ->idNode = nPhra;
		nodLZ->nChildren = 0;
		nodLZ->nextSibLZ = nodLZ->fChildLZ = NULL;

		// create an insert the new child
		insertPhrase(text, &textPos, nodLZ);
		if (text[textPos] == cutDoc){
			endDocs[nDoc] = nPhra;
			endDocsPos[nDoc] = textPos;
			if (nDoc % 10000 == 0)
				cout << "nDoc = " << nDoc << ", textPos = " << textPos << ", nPhra = " << nPhra << endl;
			nDoc++;

			// we insert a special phrase for each terminal symbol
			RevNode *newNode = new RevNode();
			newNode->idNode = nPhra;
			newNode->symbol = cutDoc;
			newNode->nChildren = 0;
			newNode->fict = newNode->uPath = false;
			newNode->nodeLZ = nodLZ;
			newNode->nextSibRev = revTrie->fChildRev;
			revTrie->fChildRev = newNode;
			(revTrie->nChildren)++;

			if(textPos >= n)
				break;
		}else
			insertRevPhrase(text, posIni, textPos, nodLZ);
		textPos++;
	}
	cout << "Final nDoc = " << nDoc << ", textPos = " << textPos << ", nPhra = " << nPhra << endl;
	if(nDoc != nDocs){
		cout << "ERROR. Different number of documents, nDocs = " << nDocs << " != nDoc = "  << nDoc << endl;
		exit(1);
	}
	nPhra++;
	cout << "LZTrie OK ! (n'=" << nPhra << ")" << endl;
	lgPhr = ceilingLog64(nPhra, 2);

	sizeNode = nPhra*lgPhr/W64;
	if ((nPhra*lgPhr)%W64)
		sizeNode++;
	*nodeArr = new ulong[sizeNode];
	for(textPos=0; textPos<sizeNode; textPos++)
		(*nodeArr)[textPos] = 0;
	sizeNode *= sizeof(ulong);

	// sort nodes in each trie...
	sortLZNodes(this->lzTrie);
	cout << "LZNodes sorted !"<< endl;

	sortRevNodes(this->revTrie);
	cout << "RevNodes sorted !"<< endl;

	// set preorder values for each node in LZTrie...
	ulong pre = 0;
	cout << "To set preorder values for each node in LZTrie..."<< endl;
	setPreorderValues(this->lzTrie, &pre);
	if (nPhra != pre){
		cout << "total LZTrie preorder's = " << pre << " != nPhra = " << nPhra << endl;
		exit(1);
	}

	pre = 0;
	cout << "To generate NodeArray..."<< endl;
	genNodeArray(this->revTrie, &pre, *nodeArr);
	if (nPhra != pre){
		cout << "total RevTrie preorder's for true nodes = " << pre << " != nPhra = " << nPhra << endl;
		exit(1);
	}
}

void LZ78Tries64::insertPhrase(uchar* text, ulong* textPos, LZNode* newNode){
	uchar letter = text[*textPos];
	LZNode* father = lzTrie;
	bool doesntInsert = true;

	while(doesntInsert){
		if(father->nChildren == 0){						// 1.- father without any children
			father->hTabLZ = new LZNode*[LENHASH];
			for(uint i=0; i<LENHASH; i++)
				father->hTabLZ[i] = NULL;
			father->nChildren = 1;
			newNode->symbol = letter;
			father->hTabLZ[cPosTb[letter]] = newNode;
			doesntInsert = false;
		}else{
			suint pos = cPosTb[letter];
			if(father->hTabLZ[pos] == NULL){			// 2.- father without any child for this symbol 'letter'
				newNode->symbol = letter;
				father->hTabLZ[pos] = newNode;
				(father->nChildren)++;
				doesntInsert = false;
			}else{
				LZNode* child = father->hTabLZ[pos];

				if(child->symbol == letter){			// 3.- There is a child for this symbol 'letter'
					father = child;
					(*textPos)++;
					letter = text[*textPos];
					if (letter == cutDoc)
						break;
				}else{									// 4.- There are other children in this slot but with different symbol
					while(child->nextSibLZ && child->nextSibLZ->symbol != letter)
						child = child->nextSibLZ;

					if(child->nextSibLZ && child->nextSibLZ->symbol == letter){
						father = child->nextSibLZ;
						(*textPos)++;
						letter = text[*textPos];
						if (letter == cutDoc)
							break;
					}else{
						(father->nChildren)++;
						newNode->symbol = letter;
						child->nextSibLZ = newNode;
						doesntInsert = false;
					}
				}
			}
		}
	}

	if (letter == cutDoc){	// special case, end of document --> we will insert at the beginning
		newNode->symbol = cutDoc;
		(father->nChildren)++;
		if(father->fChildLZ == NULL)
			father->fChildLZ = newNode;
		else{
			newNode->nextSibLZ = father->fChildLZ;
			father->fChildLZ = newNode;
		}

		if (father->hTabLZ == NULL){
			father->hTabLZ = new LZNode*[LENHASH];
			for(uint i=0; i<LENHASH; i++)
				father->hTabLZ[i] = NULL;
		}
	}
}

void LZ78Tries64::insertRevPhrase(uchar* text, ulong posIni, ulong posFin, LZNode* nodeLZ){
	ulong cont, len = posFin - posIni + 1;
	uchar letter;
	RevNode *auxNode;
	RevNode *father = revTrie;
	RevNode *newNode = new RevNode();
	newNode->idNode = nPhra;
	newNode->nChildren = 0;
	newNode->fict = newNode->uPath = false;
	newNode->nodeLZ = nodeLZ;
	newNode->nextSibRev = newNode->fChildRev = NULL;

	for(cont=0; cont < len; cont++){
		letter = text[posFin-cont];

		if(father->nChildren == 0 || father->hTabRev[cPosTb[letter]] == NULL){
			if(father->nChildren == 0){
				father->hTabRev = new RevNode*[LENHASH];
				for(uint i=0; i<LENHASH; i++)
					father->hTabRev[i] = NULL;
			}
			for(;cont<len-1; cont++){
				(father->nChildren)++;
				auxNode = new RevNode();
				auxNode->idNode = bad--;
				auxNode->uPath = false;
				auxNode->nChildren = 0;
				auxNode->fict = true;
				auxNode->symbol = letter;
				auxNode->nextSibRev = auxNode->fChildRev = NULL;
				father->hTabRev[cPosTb[letter]] = auxNode;
				father = auxNode;
				father->hTabRev = new RevNode*[LENHASH];
				for(uint i=0; i<LENHASH; i++)
					father->hTabRev[i] = NULL;
				letter = text[posFin-cont-1];
			}
			(father->nChildren)++;
			newNode->symbol = letter;
			father->hTabRev[cPosTb[letter]] = newNode;
		}else{
			RevNode* child = father->hTabRev[cPosTb[letter]];

			if(child->symbol == letter)			// 3.- There is a child for this symbol 'letter'
				if (child->fict && cont == len-1){
					delete newNode;				// child is transformed in true node
					child->fict = false;
					child->idNode = nPhra;
					child->nodeLZ = nodeLZ;
				}else{
					father = child;
					if (cont == len-1){
						newNode->symbol = letter;
						father->nextSibRev = newNode;
					}
				}
			else{									// 4.- There are other children in this slot but with different symbol
				while(child->nextSibRev && child->nextSibRev->symbol != letter)
					child = child->nextSibRev;

				if(child->nextSibRev && child->nextSibRev->symbol == letter){
					if(child->nextSibRev->fict && cont == len-1){
						delete newNode;				// child->nextSibRev is transformed in true node
						child->nextSibRev->fict = false;
						child->nextSibRev->idNode = nPhra;
						child->nextSibRev->nodeLZ = nodeLZ;
					}else
						father = child->nextSibRev;
				}else{
					(father->nChildren)++;
					if(cont<len-1){
						auxNode = new RevNode();
						auxNode->idNode = bad--;
						auxNode->uPath = false;
						auxNode->nChildren = 0;
						auxNode->fict = true;
						auxNode->symbol = letter;
						auxNode->nextSibRev = auxNode->fChildRev = NULL;
						child->nextSibRev = auxNode;
						father = auxNode;
					}else{
						newNode->symbol = letter;
						child->nextSibRev = newNode;
					}
				}
			}
		}
	}
}

void LZ78Tries64::sortLZNodes(LZNode* nod){
	LZNode** lNodesLZ = new LZNode*[SIGMA];
	LZNode* aChild;
	LZNode* aux;
	uint i, cH0;
	for (i=0; i<SIGMA; i++)
		lNodesLZ[i] = NULL;

	cH0=0;
	if (nod->hTabLZ){
		for (i=0; i<LENHASH; i++){
			if(nod->hTabLZ[i]){
				aChild = nod->hTabLZ[i];
				while(aChild){
					if(aChild->symbol != cutDoc){
						cH0++;
						lNodesLZ[aChild->symbol] = aChild;
					}
					aChild = aChild->nextSibLZ;
				}
			}
		}
	}
	if (nod->fChildLZ){ // If nod->fChildLZ is NOT NULL then nod->fChildLZ->symbol == cutDoc, and also all its sibling
		aChild = nod->fChildLZ;
		cH0++;
		if(aChild->symbol != cutDoc){
			cout << "EEEEEEERRRRRRRRRRR 0, fChildLZ != cutDoc" << endl;
			exit(1);
		}
		while(aChild->nextSibLZ){
			if(aChild->symbol != cutDoc){
				cout << "EEEEEEERRRRRRRRRRR 0, fChildLZ != cutDoc" << endl;
				exit(1);
			}
			cH0++;
			aChild = aChild->nextSibLZ;
		}
		i=0;
	}else{
		for (i=0; i<SIGMA && lNodesLZ[i] == NULL; i++);
		nod->fChildLZ = lNodesLZ[i];
		aChild = nod->fChildLZ;
		i++;
	}

	if (cH0 != nod->nChildren){
		cout << "ERR. cH0 = " << cH0 << " != nod->nChildren = " << nod->nChildren << ", Id: " << nod->idNode << endl;
		exit(0);
	}

	aChild->nextSibLZ = NULL;
	for (; i<SIGMA; i++){
		if (lNodesLZ[i]){
			aux = lNodesLZ[i];
			aux->nextSibLZ = NULL;
			aChild->nextSibLZ = aux;
			aChild = aux;
		}
	}
	aChild->nextSibLZ = NULL;
	if (nod->hTabLZ)
		delete [] nod->hTabLZ;

	if (nod->nChildren){
		aChild = nod->fChildLZ;
		while(aChild && aChild->symbol == cutDoc)
			aChild = aChild->nextSibLZ;

		while(aChild){
			if(aChild->nChildren)
				sortLZNodes(aChild);
			aChild = aChild->nextSibLZ;
		}
	}
	delete [] lNodesLZ;
}

void LZ78Tries64::sortRevNodes(RevNode* nod){
	//cout << "Rev sort " << nod->idNode << " (" << nod->nChildren << "):" << endl;
	RevNode** lNodesRev = new RevNode*[SIGMA];
	RevNode* child;
	uint i;
	for (i=0; i<SIGMA; i++)
		lNodesRev[i] = NULL;

	uint cH0 = 0;
	if (nod->hTabRev){
		for (i=0; i<LENHASH; i++){
			if(nod->hTabRev[i]){
				child = nod->hTabRev[i];
				while(child){
					cH0++;
					lNodesRev[child->symbol] = child;
					child = child->nextSibRev;
				}
			}
		}
	}
	if (nod->fChildRev){ // If nod->fChildRev is NOT NULL then nod->fChildRev->symbol == cutDoc, and also all its sibling
		child = nod->fChildRev;
		if(child->symbol != cutDoc){
			cout << "EEEEEEERRRRRRRRRRR 1" << endl;
			exit(1);
		}
		cH0++;
		while(child->nextSibRev){
			cH0++;
			if(child->symbol != cutDoc){
				cout << "EEEEEEERRRRRRRRRRR 1" << endl;
				exit(1);
			}
			child = child->nextSibRev;
		}
		i=0;
	}else{
		for (i=0; i<SIGMA && lNodesRev[i] == NULL; i++);
		nod->fChildRev = lNodesRev[i];
		child=nod->fChildRev;
		i++;
	}

	if (cH0 != nod->nChildren){
		cout << "ERROR in sortRevNodes cH0 = " << cH0 << " != nod->nChildren = " << nod->nChildren << ", Id: " << nod->idNode << endl;
		exit(0);
	}

	for (; i<SIGMA; i++){
		if (lNodesRev[i]){
			child->nextSibRev = lNodesRev[i];
			child = lNodesRev[i];
		}
	}
	child->nextSibRev = NULL;
	if (nod->hTabRev)
		delete [] nod->hTabRev;

	if (nod->nChildren){
		child = nod->fChildRev;
		while(child->symbol == cutDoc)
			child = child->nextSibRev;

		while(child){
			if (child->nChildren)
				sortRevNodes(child);
			child = child->nextSibRev;
		}
	}
	delete [] lNodesRev;
}

// this assigns the preorder value for each LZTrie's node
void LZ78Tries64::setPreorderValues(LZNode* nod, ulong *pre){
	if (nod){
		//cout << "nod " << nod->idNode << " , pre: " << *pre << endl;
		LZNode* p = nod->fChildLZ;
		nod->preorder = *pre;
		(*pre)++;

		for(uint i=0; i<nod->nChildren; i++){
			setPreorderValues(p, pre);
			p = p->nextSibLZ;
		}
		if(p){
			cout << "ERROR in LzTrie: nod->nChildren = " << nod->nChildren << ". But there is another node p = " << p->idNode << endl;
			exit(1);
		}
	}
}

// create Node array
void LZ78Tries64::genNodeArray(RevNode* nod, ulong *preRev, ulong *Node){
	RevNode* p = nod->fChildRev;
	if (nod->idNode == 0){				// mejora: no preguntar por esto, hacerlo al inicio y luego llamar al metodo
		setNum64(Node, 0, lgPhr, 0);
		(*preRev)++;
	}else{
		if(nod->fict == false){
			setNum64(Node, (*preRev)*lgPhr, lgPhr, nod->nodeLZ->preorder);
			(*preRev)++;
		}
	}

	for(uint i=0; i<nod->nChildren; i++){
		genNodeArray(p, preRev, Node);
		p = p->nextSibRev;
	}
	if(p){
		cout << "ERROR in RevTrie: nod->nChildren = " << nod->nChildren << ". But there is another node p = " << p->idNode << endl;
		exit(1);
	}
}

//	count the fictitious nodes and extra fictitious nodes
void LZ78Tries64::countFictNodRevTrie(RevNode *node){
	RevNode *child;
	uint i;

	if (node->fict){
		// if I and my unique child are fictitious nodes --> count as extra fictitious node
		if (node->nChildren == 1 ){
			if (node->fChildRev->fict)
				nExtraFict++;
			else
				nFictU++;
		}else
			nFictU++;
	}

	if (node->nChildren){
		child = node->fChildRev;
		for(i=0; i<node->nChildren; i++){

			if(child == NULL){
				cout << "EEEEEEEEEEEEEEERRRRRRRRRRRRR child == NULL" << endl;
				cout << "i = " << i << ", node->nChildren = " << node->nChildren << ", id = " << node->idNode << endl;
				exit(1);
			}

			countFictNodRevTrie(child);
			child = child->nextSibRev;
		}
	}
}

//	create 	PRev, LbRev and LbRevF structures
void LZ78Tries64::genDFUDSSeqRevTrie(RevNode *node, ulong* PRev, ulong *pos, uchar *LbRev, ulong *pLbRev){
	RevNode *child;
	uint i;

	// set open parentheses
	for(i=0, child = node->fChildRev; i<node->nChildren; i++, (*pos)++, child=child->nextSibRev, (*pLbRev)++){
		LbRev[*pLbRev] = child->symbol;
		setBit64(PRev, *pos);
		//cout << "(";
	}

	// set close parenthesis
	cleanBit64(PRev, *pos);
	(*pos)++;
	//cout << ")";

	// recursive call for all children
	child = node->fChildRev;
	for(i=0; i<node->nChildren; i++){
		genDFUDSSeqRevTrie(child, PRev, pos, LbRev, pLbRev);
		child = child->nextSibRev;
	}
}

// return the preorder value of the most right leaf of this subtree in a LZTrie
ulong LZ78Tries64::getLastLZLeaf(LZNode *node){
	while(node->nChildren){
		node = node->fChildLZ;
		while(node->nextSibLZ)
			node = node->nextSibLZ;
	}

	return node->preorder;
}

// it returns the length of the subtree rooted at node
void LZ78Tries64::getLengthSubTree(RevNode *nod, ulong *len){
	(*len) += nod->nChildren;

	RevNode* p = nod->fChildRev;
	for(uint i=0; i<nod->nChildren; i++){
		getLengthSubTree(p, len);
		p = p->nextSibRev;
	}
}

//////////////////////////////////////////////////////////////////////////////7

void LZ78Tries64::listNodeLZ(LZNode* nod, ulong* preorder, uint deph){
	LZNode* p = nod->fChildLZ;
	uint i;
	uchar ss = nod->symbol;
	if(ss == cutDoc) ss = '$';
	if (*preorder){
		for (i=0; i<deph;i++) cout << "-";
		cout << "pre " << *preorder << ", '" << ss << "' Id: " << nod->idNode << ", children " << nod->nChildren << endl;
	}else
		cout << "pre 0, '#' Id: " << nod->idNode << ", children " << nod->nChildren << endl;
	(*preorder)++;
	for(i=0; i<nod->nChildren; i++){
		listNodeLZ(p, preorder, deph+1);
		p = p->nextSibLZ;
	}
}

void LZ78Tries64::listNodeRev(RevNode* nod, ulong* preorder, uint deph){
	RevNode* p = nod->fChildRev;
	uint i;
	uchar ss = nod->symbol;
	if(ss == cutDoc) ss = '$';
	if (*preorder){
		for (i=0; i<deph;i++) cout << "-";
		cout << "pre " << *preorder << ", '" << ss << "' Id: " << nod->idNode << ", children " << nod->nChildren << ", fict " << nod->fict << ", nOccT1: " << nod->nOccT1 << endl;
	}else
		cout << "pre 0, '#' Id: " << nod->idNode << ", children " << nod->nChildren << ", fict " << nod->fict << ", nOccT1: " << nod->nOccT1 << endl;
	(*preorder)++;
	for(i=0; i<nod->nChildren; i++){
		listNodeRev(p, preorder, deph+1);
		p = p->nextSibRev;
	}
}

// list LZTrie in preorder way
void LZ78Tries64::listLZTrie(){
	ulong preorder = 0;
	listNodeLZ(this->lzTrie, &preorder, 0);
}

// list RevTrie in preorder way
void LZ78Tries64::listRevTrie(){
	ulong preorder = 0;
	listNodeRev(this->revTrie, &preorder, 0);
}


LZ78Tries64::~LZ78Tries64() {
}

void LZ78Tries64::destroyNodeTrie(LZNode* nod){
	LZNode* p = nod->fChildLZ;
	LZNode* q = NULL;

	// delete all children of p
	for(uint i=0; i<nod->nChildren; i++){
		q = p->nextSibLZ;
		destroyNodeTrie(p);
		p = q;
	}

	// now destroy nod
	delete nod;
}

