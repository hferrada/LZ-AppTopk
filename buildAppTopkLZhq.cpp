//============================================================================
// Name        : buildAppTopkLZhq.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "TopkLZhq.h"
#include <ConfigFile.h>

bool TRACE = false;			// true: print all details for console
bool TEST = false;			// true: apply exhaustive test
uint N_REP = 10;
//ulong MAX_SIZE; 			// 1073741824 bytes = 1 GB

// Structure with all globals parameters program
typedef struct {
	string configFile;		// properties file

	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	char cutDoc;			// symbol to separate documents
	uint g;					// minimum length of document segments to store a topk list in a node. Maximum 8 bits
	TopkLZhq *index;

	ulong* EndDocs;			// this store the final phrase number for each document. This is no include in the final structure

	char inputFile[200];	// list of files
	bool lowercase;			// 1: transform to lowercase
	bool filesInList;			// 1: list of files, 0: Unique file
	char boundSymbol;		// original symbol delimiter of documents when we read all documents in 1 file.
	char dirStore[200];		// directory to save/load the data structure (files *.tk)

	// The following data structure are only for test !!
	ulong* patterns;
} ParamProgram;

void testSearchPatterns_CSA_LZ(ParamProgram *par);
bool searchPatternInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, ulong *occt2, ulong *occt3, ulong *fqCSA, ulong *fqCSAall);

int main(int argc, char *argv[]) {
	ParamProgram *par = new ParamProgram();
	char fileName[300];

	if(argc != 2){
		cout << "ERRORR !! " << endl;
		cout << "buildAppTopkLZhq's usage requires the Properties File as parameter !! " << endl;
		cout << "Example for the file 'config.txt' in the same directory: ./buildAppTopkLZhq config.txt" << endl;
		exit(1);
	}
	par->configFile = string(argv[1]);
	cout << "Congif file: " << par->configFile << endl;
	ConfigFile cf(par->configFile);

	TRACE = cf.Value("GLOBALS","TRACE");
	TEST = cf.Value("GLOBALS","TEST");
	N_REP = cf.Value("GLOBALS","N_REP");
	//MAX_SIZE = cf.Value("GLOBALS","MAX_SIZE"); 	// for trec

	par->g = cf.Value("TOPK","g");
	strcpy(par->inputFile, ((string)(cf.Value("TOPK","inputFile"))).c_str());
	par->filesInList = cf.Value("TOPK","filesInList");
	strcpy(par->dirStore, ((string)(cf.Value("TOPK","dirStore"))).c_str());
	par->lowercase = cf.Value("TOPK","lowercase");
	par->boundSymbol = cf.Value("TOPK","boundSymbol");
	par->cutDoc = cf.Value("TOPK","cutDoc");

	cout << "buildAppTopkLZ parameters..." << endl;
	cout << "g = " << (uint)par->g << endl;
	cout << "Input file: " << par->inputFile << endl;
	cout << "Files in list: " << par->filesInList << endl;		// 0:archivo unico para todocl por ejemplo
	cout << "dirStore: " << par->dirStore << endl;
	cout << "lowercase: " << par->lowercase << endl;
	cout << "boundSymbol code: " << (int)(par->boundSymbol) << endl;
	cout << "cutDoc code: " << (int)(par->cutDoc) << endl;
	//cout << "MAX_SIZE: " << MAX_SIZE << endl;		// for trec

	TopkLZhq::TRACE = TRACE;
	TopkLZhq::TEST = TEST;
	//par->index = new TopkLZhq(par->g, par->inputFile, par->fileType, par->cutDoc, par->lowercase, par->dirStore, par->boundSymbol, MAX_SIZE); // for trec
	par->index = new TopkLZhq(par->g, par->inputFile, par->filesInList, par->cutDoc, par->lowercase, par->dirStore, par->boundSymbol);
	par->n = par->index->n;

	par->EndDocs = par->index->EndDocs;
	cout << "____________________________________________________" << endl;
	cout << "***  Index size " << par->index->sizeDS << " bytes = " << (float)par->index->sizeDS*8.0/(float)par->n << " bpc" << endl;
	cout << "====================================================" << endl;

	par->index->saveDS(true);

	if (TEST){
		cout << " Data Structures only for test searchPattern (sep_rrr and FMI)..." << endl;

		par->index->sep_rrr = rrr_vector<127>(par->index->separator_b);
		par->index->sep_rank = rrr_vector<127>::rank_1_type(&(par->index->sep_rrr));
		std::cout << " sep_rrr uses " << size_in_bytes(par->index->sep_rrr) << " Bytes" << endl;

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rrr.test");
		store_to_file(par->index->sep_rrr, fileName);

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sep_rank.test");
		store_to_file(par->index->sep_rank, fileName);

		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "EndDocs.test");
		ofstream os2 (fileName, ios::binary);
		os2.write((const char*)par->EndDocs, par->index->nDocs*sizeof(ulong));					// save LbRev[]
		if(TRACE) cout << " .- EndDocs[] " << par->index->nDocs*sizeof(ulong) << " Bytes" << endl;
		os2.close();

		// load Sequence...
		strcpy(fileName, "");
		strcpy(fileName, par->dirStore);
		strcat(fileName, "sequence.test");
		ifstream is(fileName, ios::binary);
		par->seq = new uchar[par->n];
		is.read((char*)par->seq, par->n*sizeof(uchar));
		is.close();

		cout << "Running test searchPattern.." << endl;
		testSearchPatterns_CSA_LZ(par);
		cout << "Test searchPattern OK !!" << endl;
	}

	par->index->~TopkLZhq();

	cout << "$$$$$$$$$$$$$$$$$$$$$" << endl;
	return 0;
}

void createPatterns(ParamProgram *par, uint m, ulong repeat){
	ulong i, j, k;
	bool eq;

	//cout << "Creating patterns of length m=" << m << endl;
	for(k=0; k<repeat; k++){
		eq = true;
		while(eq){
			i = (rand() % (par->n-(m+1)))+1;
			for(j=i; j<i+(m+1); j++){
				if (par->seq[j] == par->cutDoc){
					i = (rand() % (par->n-(m+1)))+1;
					j = i-1;
				}
				if(j==0) break;
				else{
					if (j>i && par->seq[j-1] != par->seq[j])
						eq = false;
				}
			}
		}
		par->patterns[k] = i;
		//cout << "["<<k<<"] " << i << endl;
	}
	//cout << "Patterns created !!" << endl;
}

// YOU NEED A CSA FOR TEST SEARCHES !!
bool searchPatternInFMI(ParamProgram *par, uchar* pat, uint m, ulong *occt1, ulong *occt2, ulong *occt3, ulong *fqCSA, ulong *fqCSAall){
	ulong r, o1, o2, o3, doc;
	string query = string((char *)pat);
	size_t occs = sdsl::count(par->index->fmi, query.begin(), query.begin()+m);
	auto locations = locate(par->index->fmi, query.begin(), query.begin()+m);
	//cout << "Total occurrences found with FMI : " << occs << endl;

	o1 = o2 = o3 = 0;
	//cout << "locations..." << endl;
	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << endl;
		//cout << " : " << par->sep_rrr[locations[i]] << par->sep_rrr[locations[i]+1] << par->sep_rrr[locations[i]+2] << endl;
		//cout << " : " << par->index->seq[locations[i]] << par->index->seq[locations[i]+1] << par->index->seq[locations[i]+2] << endl;
		r = par->index->sep_rank.rank(locations[i]+m) - par->index->sep_rank.rank(locations[i]);
		ulong rank = par->index->sep_rank.rank(locations[i]+1);
		doc = par->index->searchDocument(rank);
		(fqCSAall[doc])++;
		if(par->index->sep_rrr[locations[i]]){
			if (r==1){
				o1++;
				(fqCSA[doc])++;
				//cout << locations[i] << ",(t1)doc:" << doc << endl;
			}else{
				if (r==2){
					o2++;
					//cout << locations[i] << ",(t2)doc:" << doc << endl;
				}else{
					o3++;
					//cout << locations[i] << ",(t3)doc:" << doc << endl;
				}
			}
		}else{
			if (r==0){
				o1++;
				(fqCSA[doc])++;
				//cout << locations[i] << ",(t1)doc:" << doc << endl;
			}else{
				if (r==1){
					o2++;
					//cout << locations[i] << ",(t2)doc:" << doc << endl;
				}else{
					o3++;
					//cout << locations[i] << ",(t3)doc:" << doc << endl;
				}
			}
		}
	}
	//cout << endl;
	*occt1 = o1;
	*occt2 = o2;
	*occt3 = o3;

	return true;
}

bool searchPatternRevTrie(ParamProgram *par, uchar* pat, uint m, long int *x){
	LZ78Tries64::RevNode* p = par->index->tries->revTrie;
	LZ78Tries64::RevNode *node = p->fChildRev;
	m--;
	while(node){
		if (pat[m] == node->symbol){
			if (m == 0){
				*x =node->idNode;
				return true;
			}else{
				m--;
				if(node->fict == true && node->uPath == true){
					// check unary path...
					uint i;
					for (i=0; i<node->lenUPath && m>0; i++, m--){
						if (pat[m] != node->uSymbols[i])
							return false;
					}
					if (m == 0){
						if (i<node->lenUPath){
							if (pat[0] == node->uSymbols[i]){
								*x =node->idNode;
								return true;
							}else
								return false;
						}
					}else{
						if(i < node->lenUPath)
							return false;
					}
				}
				node = node->fChildRev;
			}
		}else{
			if (pat[m] > node->symbol)
				node = node->nextSibRev;
			else
				return false;
		}

	}
	return false;
}

void createTopK_Test(ParamProgram *par, uint *docList, ulong *frqList, ulong *fqCSA, uint kst){
	ulong i, j, x, lenQ;
	//cout << "nDocs " << par->index->nDocs << ", kst=" << kst << endl;
	uint *maxQ = new uint[par->index->nDocs+1]; // queue for the D documents

	/*cout << "fqCSA:" << endl;
	for(i=0; i<par->index->nDocs; i++){
		cout << fqCSA[i] << " ";
	}
	cout << endl;*/

	for(i=0; fqCSA[i]==0; i++);
	maxQ[1] = i;	// first document with frequency > 0
	lenQ = 2;

	for(i++; i<par->index->nDocs; i++){
		// go down in maxQ
		if (fqCSA[i]){
			maxQ[lenQ] = i;
			j = lenQ;
			lenQ++;

			// j must to climb...
			while (j > 1 && fqCSA[maxQ[j/2]] < fqCSA[maxQ[j]]) {
				x = maxQ[j/2];
				maxQ[j/2] = maxQ[j];
				maxQ[j] = x;
				j /= 2;
			}
		}
	}

	/*cout << "maxQ:" << endl;
	for(i=1; i<lenQ; i++)
		cout << maxQ[i] << " ";
	cout << endl;*/

	// make answer...
	for(i=0; i<kst && lenQ>1; i++){
		docList[i] = maxQ[1];
		frqList[i] = fqCSA[maxQ[1]];
		lenQ--;
		j=2;
		x = maxQ[lenQ];

		// maxQ[j] must to go down...
		while (j < lenQ) {
			if(j+1<lenQ){
				if(fqCSA[maxQ[j+1]]>fqCSA[maxQ[j]])
					j++;
			}
			if(fqCSA[maxQ[j]] > fqCSA[x])
				maxQ[j/2] = maxQ[j];
			else{
				break;
			}
			j*=2;
		}
		maxQ[j/2] = x;
	}
	delete [] maxQ;
	//cout << "ok" << endl;

	/*cout << "**** docList:" << endl;
	for(i=0; i<par->kDocs && frqList[i]; i++)
		cout << docList[i] << " ";
	cout << endl;
	for(i=0; i<par->kDocs && frqList[i]; i++)
		cout << frqList[i] << " ";
	cout << endl;*/
}

void testSearchPatterns_CSA_LZ(ParamProgram *par){
	ulong j, k, t, pv, x;
	uint m, l, M=10;
	ulong occ1, occ2, occ3, diffK, cWrong;
	bool found, isInRev;
	uchar *pat = new uchar[M+1];
	long int nod1;
	uint *docList = new uint[par->index->nDocs+1];
	uint *docList_test = new uint[par->index->nDocs+1];
	ulong *frqList_test = new ulong[par->index->nDocs+1];
	ulong *fqCSAall = new ulong[par->index->nDocs+1];
	ulong *fqCSA = new ulong[par->index->nDocs+1];
	par->patterns = new ulong[N_REP];

	uint MAX_K = 10, PLUS_K = 1, INI_K = 1;
	if (par->index->nDocs >= 100){
		MAX_K = 100;
		PLUS_K = 20;
		INI_K = 20;
	}else{
		if (MAX_K > par->index->nDocs)
			MAX_K = par->index->nDocs;
	}

	for (uint kst=INI_K; kst<=MAX_K; kst+=PLUS_K){
		cout << " Test for k = " << kst << "..." << endl;
		for (m=2; m<=M; m+=2){
			createPatterns(par, m, N_REP);
			for (j=0; j<par->index->nDocs+1; j++)
				par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j] = par->index->dictDst[j] = 0;

			//cout << " ... test in patterns of length m = " << m << endl;

			for (t=cWrong=0; t<N_REP; t++){
				//if (kst == 3 && t==27 )
				//	TRACE = true;
				for (l=0; l<m; l++)
					pat[l] = par->seq[par->patterns[t]+l];
				pat[m] = '\0';
				for(k=0; k<=par->index->nDocs; k++)
					fqCSA[k] = docList_test[k] = frqList_test[k] = fqCSAall[k] = 0;

				searchPatternInFMI(par, pat, m, &occ1, &occ2, &occ3, fqCSA, fqCSAall);
				if (TRACE){
					cout << "pat=[" << pat << "], m=" << m << ", test=" << t << ", k*=" << kst << endl;
					cout << "occ1 = " << occ1 << " ..." << endl;
					for(k=0; k<par->index->nDocs; k++){
						if(fqCSA[k] > 0)
							cout << k << "(" << fqCSA[k] << ") ";
					}
					cout << endl;
					cout << "Total occ = " << occ1+occ2+occ3 << " ..." << endl;
					for(k=0; k<par->index->nDocs; k++){
						if(fqCSAall[k] > 0)
							cout << k << "(" << fqCSAall[k] << ") ";
					}
					cout << endl;
				}

				// TEST OCCURRENCES TYPE 1...
				nod1 = x = 1;
				//if(m==2 && t==1 && kst==2)
				//	cout << "";
				found = searchPatternRevTrie(par, pat, m, &nod1);
				isInRev = par->index->searchPattern_Rev(pat, m, &x, &pv);
				if (found != isInRev){
					if (found)
						cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", is in the node " << nod1 << " and NOT FOUND in RMM representation!!" << endl;
					else
						cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", NOT FOUND in RevTrie, but IS FOUND in RMM representation!!" << endl;
					exit(1);
				}
				if (found == false && occ1){
					cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", NOT FOUND in RevTrie, but IS FOUND in the FMI with " << occ1 << " occurrences !!" << endl;
					exit(1);
				}
				if (found && occ1==0){
					cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", is in the node " << nod1 << " but NOT FOUND int the FMI!!" << endl;
					exit(1);
				}
				if (found == false) continue;

				createTopK_Test(par, docList_test, frqList_test, fqCSAall, kst);

				if (TRACE){
					cout << "Docs and Frequencies..." << endl;
					for(k=0; k<kst && frqList_test[k]; k++)
						cout << k+1 << ", (" << docList_test[k] << "," << frqList_test[k] << ") " << endl;
				}

				if (isInRev){
					// topk approximate...
					diffK = par->index->topKDocument(pat, m, docList, kst);
					if (TRACE){
						cout << "diffK= " << diffK << endl;
						cout << "List App..." << endl;
						for(k=0; k<diffK; k++)
							cout << k+1 << ", (" << docList[k] << ") " << endl;
					}
					for(k=0; k<diffK; k++){
						if (fqCSAall[docList[k]] == 0){
							cout << "ERROR. In pat = [" << pat << "], t = " << t << ", kst = " << kst << ", App report the document: " << docList[k] << " in position " << k << ", but it has not frequency!!" << endl;
							exit(1);
						}
					}

				}else
					cWrong++;

			}
			// All patterns is considered and have the same importance, these are independent of the their tf values !!
			//cout << "Patterns without occurrences type 1: " << cWrong << endl;
		}
	}
	delete [] docList;
	delete [] docList_test;
	delete [] frqList_test;
	delete [] fqCSA;
	delete [] fqCSAall;
	delete [] (par->patterns);
}
