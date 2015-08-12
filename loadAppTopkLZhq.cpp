/*
 * loadAppTopkLZhq.cpp
 *
 *  Created on: 15-09-2014
 *      Author: hector
 */

#include "TopkLZhq.h"
#include <ConfigFile.h>

bool PRINT = false;			// true: print all details for console
bool TEST_IND = true;			// true: apply exhaustive test
bool RUN_EXP = false;
bool CREATE_TABLE = false;
bool RECALL = true;
bool TODOCL = false;
bool GOODPATT = false;
uint REPEAT = 5;				// number of repetitions for each experiment
uint MAX_M = 10;			// maximum length pattern value to compute quality

// Structure with all globals parameters program
typedef struct {
	string configFile;		// properties file

	uchar *seq;				// original sequence (1 byte for symbol)
	ulong n;				// Length of generalize Text = T1$T2$...TD$
	char cutDoc;			// symbol to separate documents
	bool lowerCPatt;		// 1: transform the patterns to lowercase
	uint g;

	TopkLZhq *index;

	ulong* EndDocs;			// this store the final phrase number for each document. This is no include in the final structure

	bool pattFromFile;			// 0: random patterns, 1: from file (as in todoCL)
	char dirStore[300];		// directory to save/load the data structure (files *.tk)
	char dirResult[300];	// directory to save tables

	char fileTodoCL1W[300];
	char fileTodoCL2W[300];

	ulong* patterns;
	uchar **patt;
	suint* lenPat;
} ParamProgram;

void testTopk(ParamProgram *par);
void runExperimentsTwo(ParamProgram *par, uint k);
void runExperiments(ParamProgram *par);
void createPatterns(ParamProgram *par, uint m);
void loadPatterns(ParamProgram *par, char *fileQueries);
void loadPatternsAll(ParamProgram *par, char *fileQueries);
void createTopK_Test_load(ParamProgram *par, uint *docList, ulong *frqList, ulong *fqCSA, uint kst);
void createTopK_Test_load2(ParamProgram *par, ulong *frqList, ulong *fqCSA, uint kst);
void createTableQTodoCl(ParamProgram *par, uint kDocs, char *fileRes, bool goodPatterns);
void createTableQ(ParamProgram *par, uint kDocs, char *fileRes, bool goodPatterns);
void createTableQRecall(ParamProgram *par, uint kDocs);

int main(int argc, char *argv[]) {
	ParamProgram *par = new ParamProgram();
	char fileName[300];

	if(argc != 2){
		cout << "ERRORR !! " << endl;
		cout << "loadAppTopkLZhq's usage requires the Properties File as parameter !! " << endl;
		cout << "Example for the file 'config.txt' in the same directory: ./loadAppTopkLZhq config.txt" << endl;
		exit(1);
	}
	par->configFile = string(argv[1]);
	cout << "Congif file: " << par->configFile << endl;
	ConfigFile cf(par->configFile);

	PRINT = cf.Value("GLOBALS","TRACE");
	TEST_IND = cf.Value("GLOBALS","TEST");
	REPEAT = cf.Value("GLOBALS","N_REP");
	RUN_EXP = cf.Value("GLOBALS","RUN_EXP");
	MAX_M = cf.Value("GLOBALS","MAX_M");

	par->g = cf.Value("TOPK","g");
	TODOCL = cf.Value("TOPK","TODOCL");
	CREATE_TABLE = cf.Value("TOPK","CREATE_TABLE");
	GOODPATT = cf.Value("TOPK","GOODPATT");
	RECALL = cf.Value("TOPK","RECALL");
	strcpy(par->dirStore, ((string)(cf.Value("TOPK","dirStore"))).c_str());
	strcpy(par->dirResult, ((string)(cf.Value("TOPK","dirResult"))).c_str());
	par->pattFromFile = cf.Value("TOPK","pattFromFile");	// 0:random
	par->lowerCPatt = cf.Value("TOPK","lowercase");
	if(par->pattFromFile){
		strcpy(par->fileTodoCL1W, ((string)(cf.Value("TOPK","fileTodoCL1W"))).c_str());
		strcpy(par->fileTodoCL2W, ((string)(cf.Value("TOPK","fileTodoCL2W"))).c_str());
	}

	cout << "loadAppTopkLZ parameters..." << endl;
	cout << "dirStore: " << par->dirStore << endl;
	cout << "dirResult: " << par->dirResult << endl;
	cout << "patterns from file: " << par->pattFromFile << endl;
	cout << "lowercase patterns: " << par->lowerCPatt << endl;
	if(par->pattFromFile){
		cout << "fileTodoCL1W: " << par->fileTodoCL1W << endl;
		cout << "fileTodoCL2W: " << par->fileTodoCL2W << endl;
	}

	par->index = new TopkLZhq(par->dirStore, true);
	par->n = par->index->n;
	par->cutDoc = par->index->cutDoc;

	cout << "____________________________________________________" << endl;
	cout << "***  Index size " << par->index->sizeDS << " bytes = " << (float)par->index->sizeDS*8.0/(float)par->n << " bpc" << endl;
	cout << "====================================================" << endl;

	strcpy(fileName, "");
	strcpy(fileName, par->dirStore);
	strcat(fileName, "sequence.test");
	ifstream is(fileName, ios::binary);
	par->seq = new uchar[par->n];
	is.read((char*)par->seq, par->n*sizeof(uchar));
	is.close();

	strcpy(fileName, "");
	strcpy(fileName, par->dirStore);
	strcat(fileName, "sep_rrr.test");
	load_from_file(par->index->sep_rrr, fileName);

	strcpy(fileName, "");
	strcpy(fileName, par->dirStore);
	strcat(fileName, "sep_rank.test");
	load_from_file(par->index->sep_rank, fileName);
	util::init_support(par->index->sep_rank, &par->index->sep_rrr);

	strcpy(fileName, "");
	strcpy(fileName, par->dirStore);
	strcat(fileName, "EndDocs.test");
	ifstream is2(fileName, ios::binary);
	par->EndDocs = new ulong[par->index->nDocs];
	is2.read((char*)par->EndDocs, par->index->nDocs*sizeof(ulong));
	is2.close();
	par->index->EndDocs = par->EndDocs;

	strcpy(fileName, "");
	strcpy(fileName, par->dirStore);
	strcat(fileName, "fmi.test");
	load_from_file(par->index->fmi, fileName);
    cout << "Index construction complete, index requires " << size_in_mega_bytes(par->index->fmi) << " MiB." << endl;

	if (TEST_IND){
		cout << "Test Index..." << endl;
		testTopk(par);
		cout << "Test Index OK !!" << endl;
	}

	if (RECALL){
		cout << endl << "* Create Recall/Quality Table from random patterns..." << endl;
		// normal case...
		createTableQRecall(par, 10);
		createTableQRecall(par, 100);
		if (RUN_EXP)
			runExperiments(par);
	}else{
		if(CREATE_TABLE){
			char suffTable[30];
			if (TODOCL){
				cout << endl << "* ALL PATTERNS ... Create Quality Table for TodoCl (1 word) from [" << par->fileTodoCL1W << "]" << endl;
				loadPatternsAll(par, par->fileTodoCL1W);
				strcpy(suffTable, "");
				strcpy(suffTable, "QT1W_allP.tk");
				createTableQTodoCl(par, 10, suffTable, GOODPATT);
				createTableQTodoCl(par, 100, suffTable, GOODPATT);
				if (RUN_EXP){
					// This initialization is done only one time.
					for (uint i=0; i<par->index->nDocs; i++)
						par->index->dictD[i] = par->index->keyFq[i] = 0;
					runExperimentsTwo(par, 10);
					runExperimentsTwo(par, 100);
				}
				delete [] (par->patt);
				delete [] (par->lenPat);

				cout << endl << "* ALL PATTERNS ... Create Quality Table for TodoCl (2 words) from [" << par->fileTodoCL2W << "]" << endl;
				loadPatternsAll(par, par->fileTodoCL2W);
				strcpy(suffTable, "");
				strcpy(suffTable, "QT2W_allP.tk");
				createTableQTodoCl(par, 10, suffTable, GOODPATT);
				createTableQTodoCl(par, 100, suffTable, GOODPATT);
				if (RUN_EXP){
					// This initialization is done only one time.
					for (uint i=0; i<par->index->nDocs; i++)
						par->index->dictD[i] = par->index->keyFq[i] = 0;
					runExperimentsTwo(par, 10);
					runExperimentsTwo(par, 100);
				}
				delete [] par->patt;
				delete [] par->lenPat;
			}else{
				char str[100];
				strcpy(suffTable, "");

				if (GOODPATT){
					cout << endl << "* Create Quality Table for GOOD random patterns..." << endl;
					sprintf(str, "QTable_goodP_g%d_k10.txt", par->g);
					strcpy(suffTable, str);
				}else{
					cout << endl << "* Create Quality Table for NORMAL random patterns..." << endl;
					strcpy(suffTable, "");
					sprintf(str, "QTable_normalP_g%d_k10.txt", par->g);
					strcpy(suffTable, str);
				}
				createTableQ(par, 10, suffTable, GOODPATT);

				if (par->index->nDocs > 100){
					strcpy(suffTable, "");
					if (GOODPATT){
						sprintf(str, "QTable_goodP_g%d_k100.txt", par->g);
						strcpy(suffTable, str);
					}else{
						sprintf(str, "QTable_normalP_g%d_k100.txt", par->g);
						strcpy(suffTable, str);
					}
					createTableQ(par, 100, suffTable, GOODPATT);
				}

				if (RUN_EXP)
					runExperiments(par);
			}
		}else{
			if (RUN_EXP)
				runExperiments(par); // this store the patterns in par->patterns[]
		}
	}

	cout << "$$$$$$$$$$$$$$$$$$$$$" << endl;
	return 0;
}

// YOU NEED A CSA FOR TEST SEARCHES !!
uint searchPatternInFMI_all(ParamProgram *par, uchar* pat, uint m, ulong *occ, ulong *fqCSA){
	uint doc, diffReal=0;
	bool *DOCS = new bool[par->index->nDocs];
	string query = string((char *)pat);
	size_t occs = sdsl::count(par->index->fmi, query.begin(), query.begin()+m);
	auto locations = locate(par->index->fmi, query.begin(), query.begin()+m);
	//cout << "Total occurrences found with FMI : " << occs << endl;

	for(doc=0; doc<par->index->nDocs; doc++)
		DOCS[doc] = false;

	for(ulong i=0; i<occs; i++){
		//cout << locations[i] << " ";
		ulong rank = par->index->sep_rank.rank(locations[i]+1);
		doc = par->index->searchDocument(rank);
		if (DOCS[doc] == 0){
			diffReal++;
			DOCS[doc] = 1;
		}
		(fqCSA[doc])++;
	}
	//cout << endl;
	*occ = occs;

	delete [] DOCS;
	return diffReal;
}

void createTableQTodoCl(ParamProgram *par, uint kDocs, char *fileRes, bool goodPatterns){
	uint m, j, t, k, diffDocsReal, badPatt;;
	ulong occ, diffK, cWrong, pv, x;
	bool isInRev;
	uchar *patron;
	ulong *fqCSA = new ulong[par->index->nDocs+1];
	uint *docList = new uint[par->index->nDocs+1];
	ulong *frqList = new ulong[par->index->nDocs+1];
	uint *docList_csa = new uint[par->index->nDocs+1];
	ulong *frqList_csa = new ulong[par->index->nDocs+1];
	float quality, tf, tf_prima, avgM = 0.0;

	cout << "____________________________________________________" << endl;

	quality = 0.0;
	for (t=cWrong=badPatt=0; t<REPEAT; t++){
		for(k=0; k<kDocs; k++)
			docList[k] = frqList[k] = docList_csa[k] = frqList_csa[k] = 0;

		for(k=0; k<par->index->nDocs; k++)
			fqCSA[k] = 0;

		patron = par->patt[t];
		m = par->lenPat[t];

		//cout << "search in CSA/FMI...patron = [" << patron << "], m=" << m << ", test=" << t << endl;
		diffDocsReal = searchPatternInFMI_all(par, patron, m, &occ, fqCSA); // always occ > 0 because there is this pattern in the original text
		if(goodPatterns && diffDocsReal < kDocs){
			badPatt++;
			continue;
		}
		avgM += m;

		//cout << "patron=[" << patron << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
		if (PRINT){
			cout << "> pat=[" << patron << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
			if (occ){
				for(k=0; k<par->index->nDocs; k++)
					if (fqCSA[k])
						cout << k << "(" << fqCSA[k] << ") ";
				cout << endl;
			}
		}

		//cout << "pat=[" << pat << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
		x = 1;
		isInRev = par->index->searchPattern_Rev(patron, m, &x, &pv);
		for (k=0; k<kDocs; k++)
			frqList_csa[k] = frqList[k] = 0;
		// createTopK_Test for all occ in CSA, that is the general and correct topK !!
		createTopK_Test_load(par, docList_csa, frqList_csa, fqCSA, kDocs);
		if (PRINT){
			cout << "Real k: Id-tf..." << endl;
			for(k=0; k<kDocs && frqList_csa[k]; k++)
				cout << k+1 << ":" << docList_csa[k] << "-"<< frqList_csa[k] << " ";
			cout << endl;
		}
		tf = tf_prima = 0.0;
		for(k=0; k<kDocs && frqList_csa[k]; k++)
			tf += frqList_csa[k];

		if (isInRev){
			for (j=0; j<par->index->nDocs; j++)
				par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j] = 0;

			//cout << "search topKDocument for pat=[" << pat << endl;
			diffK = par->index->topKDocument(patron, m, docList, kDocs);
			if (PRINT){
				cout << "App: k, Id-tf..." << endl;
				for(k=0; k<diffK; k++)
					cout << k+1 << ", " << docList[k] << endl;
			}
			for(k=0; k<diffK; k++)
				tf_prima += fqCSA[docList[k]];

			if(tf < tf_prima){
				cout << "ERRRRR. tf = " << tf << " < tf_prima = " << tf_prima << endl;
				cout << "Real k: Id-tf..." << endl;
				for(k=0; k<kDocs && frqList_csa[k]; k++)
					cout << k+1 << ":" << docList_csa[k] << "-"<< frqList_csa[k] << " ";
				cout << endl;
				cout << "App: k, Id-tf..." << endl;
				for(k=0; k<diffK; k++)
					cout << k+1 << ", " << docList[k] << endl;

				cout << "t=" << t << endl;
				exit(0);
			}

			quality += tf_prima/tf;
		}else
			cWrong++;
	}
	REPEAT -= badPatt;
	quality /= (float)REPEAT;
	avgM /= (float)REPEAT;

	char fileTable[300];
	strcpy(fileTable, "");
	strcpy(fileTable, par->dirResult);
	strcat(fileTable, fileRes);
	FILE *fTable = fopen(fileTable, "a+");
	cout << "File " << fileTable << " opened" << endl;
	// All patterns are considered with the same importance !! these are independent of their tf values.
	cout << "There are " << cWrong << " (from " << REPEAT << ") patterns without occurrences type one" << endl;
	cout << "[k] [g] [size] [avgM] [ratio quality] = [" << kDocs << "] [" << par->index->g << "] [" <<  par->index->sizeDS*8.0/par->n << "] [" << avgM << "] [" << quality << "]" << endl;
	fprintf(fTable, "%d %d %f %f %f\n", kDocs, par->index->g, (float)par->index->sizeDS*8.0/(float)par->n, avgM, quality);
	fclose(fTable);

	delete [] docList;
	delete [] frqList;
	delete [] docList_csa;
	delete [] frqList_csa;
	delete [] fqCSA;
}

void createTableQ(ParamProgram *par, uint kDocs, char *fileRes, bool goodPatterns){
	uint m, i, j, t, k, diffDocsReal;
	ulong occ, diffK, cWrong, pv, x;
	bool isInRev, eq;
	uchar *pat = new uchar[MAX_M+1];
	ulong *fqCSA = new ulong[par->index->nDocs+1];
	uint *docList = new uint[par->index->nDocs+1];
	ulong *frqList = new ulong[par->index->nDocs+1];
	uint *docList_csa = new uint[par->index->nDocs+1];
	ulong *frqList_csa = new ulong[par->index->nDocs+1];
	float quality, tf, tf_prima;
	char fileTable[300];
	strcpy(fileTable, "");
	strcpy(fileTable, par->dirResult);
	strcat(fileTable, fileRes);
	FILE *fTable = fopen(fileTable, "a+");
	cout << "____________________________________________________" << endl;
	cout << "File " << fileTable << " opened" << endl;

	par->patterns = new ulong[REPEAT];

	for (m=2; m<=MAX_M; m++){
		//cout << "createPatterns m = " << m << endl;
		cout << " ... :) row for patterns of length m = " << m << endl;
		quality = 0.0;
		for (t=cWrong=0; t<REPEAT; t++){
			if (goodPatterns){
				diffDocsReal=0;
				while(diffDocsReal < kDocs){
					for(k=0; k<par->index->nDocs; k++)
						fqCSA[k] = 0;

					eq = true;
					while(eq){
						i = (rand() % (par->n-(m+1)))+1;
						for(j=i; j<i+(m+1); j++){
							if (par->seq[j] <= par->cutDoc){
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
					for (j=0; j<m; j++)
						pat[j] = par->seq[i+j];
					pat[m] = '\0';

					diffDocsReal = searchPatternInFMI_all(par, pat, m, &occ, fqCSA); // always occ > 0 because there is this pattern in the original text
				}
			}else{
				for(k=0; k<par->index->nDocs; k++)
					fqCSA[k] = 0;

				eq = true;
				while(eq){
					i = (rand() % (par->n-(m+1)))+1;
					for(j=i; j<i+(m+1); j++){
						if (par->seq[j] <= par->cutDoc){
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
				for (j=0; j<m; j++)
					pat[j] = par->seq[i+j];
				pat[m] = '\0';

				diffDocsReal = searchPatternInFMI_all(par, pat, m, &occ, fqCSA); // always occ > 0 because there is this pattern in the original text
			}

			// search in CSA...
			for(k=0; k<kDocs; k++)
				docList[k] = frqList[k] = docList_csa[k] = frqList_csa[k] = 0;

			//cout << "patron=[" << patron << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
			if (PRINT){
				cout << "pat=[" << pat << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
				if (fqCSA[0]){
					for(k=0; k<par->index->nDocs && fqCSA[k]; k++)
						cout << k << "(" << fqCSA[k] << ") ";
					cout << endl;
				}
			}

			//cout << "pat=[" << pat << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
			x = 1;
			isInRev = par->index->searchPattern_Rev(pat, m, &x, &pv);
			for (k=0; k<kDocs; k++)
				frqList_csa[k] = frqList[k] = 0;
			// createTopK_Test for all occ in CSA, that is the general and correct topK !!
			createTopK_Test_load(par, docList_csa, frqList_csa, fqCSA, kDocs);
			if (PRINT){
				cout << "Real k: Id-tf..." << endl;
				for(k=0; k<kDocs && frqList_csa[k]; k++)
					cout << k+1 << ":" << docList_csa[k] << "-"<< frqList_csa[k] << " ";
				cout << endl;
			}
			tf = tf_prima = 0.0;
			for(k=0; k<kDocs && frqList_csa[k]; k++)
				tf += frqList_csa[k];

			if (isInRev){
				for (j=0; j<par->index->nDocs; j++)
					par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j] = 0;

				//cout << "search topKDocument for pat=[" << pat << endl;
				diffK = par->index->topKDocument(pat, m, docList, kDocs);
				if (PRINT){
					cout << "App: k, Id-tf..." << endl;
					for(k=0; k<diffK; k++)
						cout << k+1 << ", " << docList[k] << endl;
				}
				for(k=0; k<diffK; k++)
					tf_prima += fqCSA[docList[k]];

				if(tf < tf_prima){
					cout << "ERRRRR. tf = " << tf << " < tf_prima = " << tf_prima << endl;
					cout << "Real k: Id-tf..." << endl;
					for(k=0; k<kDocs && frqList_csa[k]; k++)
						cout << k+1 << ":" << docList_csa[k] << "-"<< frqList_csa[k] << " ";
					cout << endl;
					cout << "App: k, Id-tf..." << endl;
					for(k=0; k<diffK; k++)
						cout << k+1 << ", " << docList[k] << endl;

					cout << "t=" << t << endl;
					exit(0);
				}

				quality += tf_prima/tf;
			}else
				cWrong++;

		}
		quality /= (float)(REPEAT);

		// All patterns are considered with the same importance !! these are independent of their tf values.
		//cout << "For m="<<m<<" there are " << cWrong << " patterns without occurrences type one" << endl;
		cout << "[k] [m] [ratio quality] = [" << kDocs << "] [" << m << "] [" << quality << "]" << endl;
		fprintf(fTable, "%d %d %f\n", kDocs, m, quality);
	}

	cout << "There are " << cWrong << " (from " << REPEAT << ") patterns without occurrences type one" << endl;

	fclose(fTable);
	delete [] docList;
	delete [] docList_csa;
	delete [] frqList_csa;
	delete [] frqList;
	delete [] fqCSA;
}

// ordena en 'sortedApp' la lista App por frecuencia real
void sortAppReal(uint *docList, ulong *fqCSA, uint diffK, ulong *sortedApp){
	uint i, k;
	ulong frec;

	sortedApp[0] = fqCSA[docList[0]];
	// insertion sort
	for (k=1; k<diffK; k++){
		frec = fqCSA[docList[k]];
		for (i=k; i>0 && sortedApp[i-1] < frec; i--)
			sortedApp[i] = sortedApp[i-1];
		sortedApp[i] = frec;
	}
}

void createTableQRecall(ParamProgram *par, uint kDocs){
	uint m, i, j, t, k, cValApp, diffDocsReal;
	ulong occ, diffK, cWrong, pv, x;
	bool eq;
	uchar *pat = new uchar[MAX_M+1];
	ulong *fqCSA = new ulong[par->index->nDocs+1];
	uint *docList = new uint[par->index->nDocs+1];

	ulong *frqList_app = new ulong[par->index->nDocs+1];
	ulong *frqList_real = new ulong[par->index->nDocs+1];

	float *tf = new float[par->index->nDocs];
	float *tf_app = new float[par->index->nDocs];
	float *recallTK = new float[par->index->nDocs];
	float *quality = new float[par->index->nDocs];
	ulong *sortedApp = new ulong[par->index->nDocs];

	for (m=6; m<=MAX_M; m+=4){
		for (k=0; k<par->index->nDocs; k++)
			recallTK[k] = quality[k] = 0.0;

		for (t=cWrong=0; t<REPEAT; t++){
			diffDocsReal=0;
			while(diffDocsReal < kDocs){
				for(k=0; k<par->index->nDocs; k++)
					fqCSA[k] = 0;
				eq = true;		// this while ensure to found a pattern that appear at least in kDocs documents !!
				while(eq){
					i = (rand() % (par->n-(m+1)))+1;
					for(j=i; j<i+(m+1); j++){
						if (par->seq[j] <= par->cutDoc){
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
				for (j=0; j<m; j++)
					pat[j] = par->seq[i+j];
				pat[m] = '\0';
				//cout << "search in CSA/FMI...patron = [" << patron << "], m=" << m << ", test=" << t << endl;
				diffDocsReal = searchPatternInFMI_all(par, pat, m, &occ, fqCSA); // always occ > 0 because there is this pattern in the original text
			}

			for(k=0; k<par->index->nDocs; k++){
				sortedApp[k] = docList[k] = frqList_app[k] = frqList_real[k] = 0;
				tf[k] = 0.0;
			}

			//cout << "patron=[" << patron << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
			if (PRINT){
				cout << endl << "pat=[" << pat << "], m=" << m << ", test=" << t << ", occ = " << occ << endl;
				for(k=0; k<par->index->nDocs; k++)
					if (fqCSA[k])
						cout << k << "(" << fqCSA[k] << ") ";
				cout << endl;

			}

			// createTopK_Test for all occ in CSA, that is the general and correct topK !!
			createTopK_Test_load2(par, frqList_real, fqCSA, kDocs);
			if (PRINT){
				cout << "Real k: Id-tf..." << endl;
				for(k=0; k<kDocs && frqList_real[k]; k++)
					cout << k+1 << ":" << "("<< frqList_real[k] << ") ";
				cout << endl;
			}

			for(k=0; k<kDocs; k++){
				if (frqList_real[k] > 0){
					if (k==0)
						tf[0] = frqList_real[0];
					else
						tf[k] = tf[k-1] + frqList_real[k];
				}else
					tf[k] = tf[k-1];
				//cout << "tf["<<k<<"] = " << tf[k] << endl;
			}

			x = 1;
			if (par->index->searchPattern_Rev(pat, m, &x, &pv)){
				for (j=0; j<par->index->nDocs; j++)
					par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j] = 0;

				diffK = par->index->topKDocument(pat, m, docList, kDocs);
				if (PRINT){
					cout << "* App: k, docId, freq..." << endl;
					for(k=0; k<diffK; k++)
						cout << k+1 << ", " << docList[k] << ", " << fqCSA[docList[k]] << endl;
				}
				// sort the app-list 'docList' by real tf.
				sortAppReal(docList, fqCSA, diffK, sortedApp);
				if (PRINT){
					cout << " ** Sort App: freq ..." << endl;
					for(k=0; k<diffK; k++){
						cout << k+1 << ", " << sortedApp[k] << endl;
					}
				}

				for(k=0; k<kDocs; k++){
					if(k < diffK){
						if (k == 0)
							tf_app[0] = sortedApp[0];
						else
							tf_app[k] = tf_app[k-1] + sortedApp[k];
					}else{
						tf_app[k] = tf_app[k-1];
						//cout << "Final tf_app = " << tf_app[diffK-1] << ", tf = " << tf[k] << " % = " << tf_app[diffK-1] / tf[k] << endl ;
					}
					(quality[k]) += (tf_app[k] / tf[k]);
				}

				cValApp=0;
				for(k=0; k<diffK; k++){
					if(sortedApp[cValApp] >= frqList_real[k])
						cValApp++;
					(recallTK[k]) += (float)cValApp/(float)(k+1);
				}

				for(; k<kDocs; k++){
					(recallTK[k]) += (float)cValApp/(float)(k+1);
				}

			}else
				cWrong++;
		}

		cout << "There are " << cWrong << " (from " << REPEAT << ") patterns without occurrences type one" << endl;

		char fileTable[200];
		char str[20];
		cout << "____________________________________________________" << endl;
		strcpy(fileTable, "");
		strcpy(fileTable, par->dirResult);
		sprintf(str, "Recall_g%d_k%d_m%d", par->g, kDocs, m);
		strcat(fileTable, str);
		FILE *fTable = fopen(fileTable, "w");
		cout << "File " << fileTable << " opened, columns:" << endl;
		cout << "[ranking] [recall] [tfApp/tf]" << endl;
		for(k=0; k<kDocs; k++){
			quality[k] /= REPEAT;
			recallTK[k] /= REPEAT;

			cout << k+1 << ", " << recallTK[k] << ", " << quality[k] << endl;
			fprintf(fTable, "%d %f %f\n", k+1, recallTK[k], quality[k]);
		}
		fclose(fTable);
	}

	delete [] docList;
	delete [] fqCSA;
}

bool searchPatternInFMI_load(ParamProgram *par, uchar* pat, uint m, ulong *occt1, ulong *occt2, ulong *occt3, ulong *fqCSA, ulong *fqCSAall){
	ulong r, o1, o2, o3, doc, rank;
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
		r = par->index->sep_rank.rank(locations[i]+m) - par->index->sep_rank.rank(locations[i]+1);
		rank = par->index->sep_rank.rank(locations[i]+1);
		doc = par->index->searchDocument(rank);
		(fqCSAall[doc])++;
		if (r==0){
			o1++;
			(fqCSA[doc])++;
			//cout << locations[i] << ",doc:" << doc << endl;
		}else{
			if (r==1)
				o2++;
			else{
				if (r==2 && par->index->sep_rrr[locations[i]])
					o2++;
				else
					o3++;
			}
		}
	}
	//cout << endl;
	*occt1 = o1;
	*occt2 = o2;
	*occt3 = o3;

	return true;
}

void createTopK_Test_load2(ParamProgram *par, ulong *frqList, ulong *fqCSA, uint kst){
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

void createTopK_Test_load(ParamProgram *par, uint *docList, ulong *frqList, ulong *fqCSA, uint kst){
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

void testTopk(ParamProgram *par){
	ulong j, k, t, pv, x;
	uint m, l, M=10;
	ulong occ1, occ2, occ3, diffK, cWrong;
	bool isInRev;
	uchar *pat = new uchar[M+1];
	uint *docList = new uint[par->index->nDocs+1];
	uint *docList_test = new uint[par->index->nDocs+1];
	ulong *frqList_test = new ulong[par->index->nDocs+1];
	ulong *fqCSAall = new ulong[par->index->nDocs+1];
	ulong *fqCSA = new ulong[par->index->nDocs+1];

	par->patterns = new ulong[REPEAT];

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
		cout << " Test for k = " << kst << endl;
		for (m=2; m<=M; m+=2){
			createPatterns(par, m);
			for (j=0; j<par->index->nDocs+1; j++){
				par->index->dictD[j] = par->index->keyFq[j] = par->index->keyDoc[j];
				par->index->dictDst[j] = 0;
			}

			cout << " ... test in patterns of length m=" << m << endl;
			if(false){
				pat[0] = 'a';
				pat[1] = 'p';
				pat[2] = '\0';

				diffK = par->index->topKDocument(pat, m, docList, 4);
				cout << "diffK= " << diffK << endl;
				cout << "List App..." << endl;
				for(k=0; k<diffK; k++)
					cout << k+1 << ", (" << docList[k] << ") " << endl;
				exit(0);
			}

			for (t=cWrong=0; t<REPEAT; t++){
				for (l=0; l<m; l++)
					pat[l] = par->seq[par->patterns[t]+l];
				pat[m] = '\0';
				for(k=0; k<=par->index->nDocs; k++){
					fqCSA[k] = fqCSAall[k] = frqList_test[k] = 0;
					docList_test[k] = 0;
				}

				searchPatternInFMI_load(par, pat, m, &occ1, &occ2, &occ3, fqCSA, fqCSAall);
				if (PRINT){
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
				x = 1;
				isInRev = par->index->searchPattern_Rev(pat, m, &x, &pv);
				if (isInRev == false && occ1){
					cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", NOT FOUND in RevTrie, but IS FOUND in the FMI with " << occ1 << " occurrences !!" << endl;
					exit(1);
				}
				if (isInRev && occ1==0){
					cout << "ERROR, [" << pat << "], m=" << m << ", k* =" << kst << ", test=" << t << ", is in the node with preorder " << x << " but NOT FOUND int the FMI!!" << endl;
					exit(1);
				}
				if (isInRev == false) continue;

				createTopK_Test_load(par, docList_test, frqList_test, fqCSAall, kst);

				if (PRINT){
					cout << "Docs and Frequencies..." << endl;
					for(k=0; k<kst && frqList_test[k]; k++)
						cout << k+1 << ", (" << docList_test[k] << "," << frqList_test[k] << ") " << endl;
				}

				if (isInRev){
					// topk approximate...
					diffK = par->index->topKDocument(pat, m, docList, kst);
					if (PRINT){
						cout << "diffK= " << diffK << endl;
						cout << "List App..." << endl;
						for(k=0; k<diffK; k++)
							cout << k+1 << ", (" << docList[k] << ") " << endl;
					}

					for(k=0; k<diffK; k++){
						if (fqCSAall[docList[k]] == 0){
							cout << "ERROR. In pat = [" << pat << "], t = " << t << ", App report the document " << docList[k] << " but it has not frequency!!" << endl;
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

void createPatterns(ParamProgram *par, uint m){
	ulong i, j, k;
	bool eq;

	cout << "Creating " << REPEAT << " patterns of length m=" << m << endl;
	for(k=0; k<REPEAT; k++){
		eq = true;
		while(eq){
			i = (rand() % (par->n-(m+1)))+1;
			for(j=i; j<i+(m+1); j++){
				if (par->seq[j] <= par->cutDoc){
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

void loadPatternsAll(ParamProgram *par, char *fileQueries){
	uint i, j, k, len, totLines;
	string line;
	ifstream input (fileQueries);
	assert(input.good()); 				// check open
	if (input.is_open()){
		totLines = len = 0;
		while (input.good()){
			totLines++;
			getline(input,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			len++;
		}
		input.close();
	} else{
		cout << "Unable to open file " << fileQueries << endl;
		exit(0);
	}

	cout << "Total lines: " << totLines << ", patterns selected: " << len << endl;
	REPEAT = len;
	par->patt = new uchar*[len];
	par->lenPat = new suint[len];
	ifstream input2(fileQueries);// open file
	if (input2.is_open()){
		i = 0;
		while (input2.good() && i < len){
			getline(input2,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			par->patt[i] = new uchar[k+1];
			par->lenPat[i] = k;
			for(j=0; j<k; j++){
				if (par->lowerCPatt)
					par->patt[i][j] = (uchar)tolower(line.at(j));
				else
					par->patt[i][j] = (uchar)(line.at(j));
			}
			par->patt[i][k] = '\0';
			i++;
			//cout << par->patt[i] << endl;
		}
		input2.close();
	}
	cout << i << " patterns loaded !!" << endl;
}

void loadPatterns(ParamProgram *par, char *fileQueries){
	uint i, j, k, len, numspaces, totLines;
	string line;
	ifstream input (fileQueries);
	assert(input.good()); 				// check open
	if (input.is_open()){
		totLines = len =0;
		while (input.good()){
			totLines++;
			getline(input,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			for(j=numspaces=0; j<k; j++){
				if (isspace(line[j]))
					numspaces++;
			}
			if(numspaces <= 1)
				len++;
		}
		input.close();
	} else{
		cout << "Unable to open file " << fileQueries << endl;
		exit(0);
	}

	cout << "Total lines: " << totLines << ", patterns selected: " << len << endl;
	REPEAT = len;
	par->patt = new uchar*[len];
	par->lenPat = new suint[len];
	ifstream input2(fileQueries);// open file
	if (input2.is_open()){
		i = 0;
		while (input2.good() && i < len){
			getline(input2,line);
			k = line.size();
			if (k > 99 || k < 3) continue;
			for(j=numspaces=0; j<k; j++){
				if (isspace(line[j]))
					numspaces++;
			}
			if(numspaces <= 1){
				par->patt[i] = new uchar[k+1];
				par->lenPat[i] = k;
				for(j=0; j<k; j++){
					if (par->lowerCPatt)
						par->patt[i][j] = (uchar)tolower(line.at(j));
					else
						par->patt[i][j] = (uchar)(line.at(j));
				}
				par->patt[i][k] = '\0';
				i++;
				//cout << par->patt[i] << endl;
			}
		}
		input2.close();
	}
	cout << i << " patterns loaded !!" << endl;
}

void runExperimentsOne(ParamProgram *par, uint m, uint k){
	double t, avgTime=0.0;
	char aFile[300];
	uint *docList = new uint[par->index->nDocs];

	cout << "____________________________________________________" << endl;
	cout << "Start Top-k for " << REPEAT << " patterns of length m = " << m << " and k = " << k << endl;

	for (uint j=0; j<REPEAT; j++){
		t = getTime_ms();
		par->index->topKDocument(par->seq+par->patterns[j], m, docList, k);
		avgTime += getTime_ms() - t;
	}

	avgTime /= (double)REPEAT;
	cout << "Average CPU time for execution, Search type 1 : " << avgTime*1000.0 << " Microseconds" << endl;
	cout << "Size : " << par->index->sizeDS*8.0/(float)par->n << endl;
	cout << "____________________________________________________" << endl;

	strcpy(aFile, "");
	strcpy(aFile, par->dirStore);
	strcat(aFile, "resume.topk");
	FILE *fp = fopen(aFile, "a+" );
	// [m] [k] [g] [size] [time]
	fprintf(fp, "%d %d %d %f %G\n", m, k, par->index->g, par->index->sizeDS*8.0/(float)par->n, avgTime*1000.0);
	fclose(fp);
}

void runExperimentsTwo(ParamProgram *par, uint k){
	double t, avgTime=0.0;
	char aFile[300];
	uint m;
	uchar *patron;
	uint *docList = new uint[par->index->nDocs];
	float avgM = 0.0;

	cout << "____________________________________________________" << endl;
	cout << "Start Top-k for k = " << k << endl;

	for (uint j=0; j<REPEAT; j++){
		patron = par->patt[j];
		m = par->lenPat[j];
		t = getTime_ms();
		par->index->topKDocument(patron, m, docList, k);
		avgTime += getTime_ms() - t;
		avgM += m;
	}

	avgM /= (float)REPEAT;
	avgTime /= (double)REPEAT;
	cout << "Average CPU time for execution, Search type 1 : " << avgTime*1000.0 << " Microseconds" << endl;
	cout << "Average pattern length : " << avgM << endl;
	cout << "Size : " << par->index->sizeDS*8.0/(float)par->n << endl;
	cout << "____________________________________________________" << endl;

	strcpy(aFile, "");
	strcpy(aFile, par->dirStore);
	strcat(aFile, "resume.topk");
	FILE *fp = fopen(aFile, "a+" );
	// [avgM] [k] [g] [size] [time]
	fprintf(fp, "%f %d %d %f %G\n", avgM, k, par->index->g, par->index->sizeDS*8.0/(float)par->n, avgTime*1000.0);
	fclose(fp);
}

void runExperiments(ParamProgram *par){
	cout << " RUN TOP-K DOCUMENT LISTING..." << endl;
	cout << "====================================================" << endl;
	uint j, m;

	// This initialization is done only one time.
	for (j=0; j<par->index->nDocs; j++)
		par->index->dictD[j] = par->index->keyFq[j] = 0;

	if (par->pattFromFile){
		cout << " Todocl for patterns with 1 word..." << endl;
		loadPatternsAll(par, par->fileTodoCL1W);
		for (uint kst=10; kst<=100; kst+=90)
			runExperimentsTwo(par, kst);
		delete [] par->patt;
		delete [] par->lenPat;

		cout << " Todocl for patterns with 2 words..." << endl;
		loadPatternsAll(par, par->fileTodoCL2W);
		for (uint kst=10; kst<=100; kst+=90)
			runExperimentsTwo(par, kst);
		delete [] par->patt;
		delete [] par->lenPat;

	}else{
		par->patterns = new ulong[REPEAT];

		m=6;
		uint kst=10;
		createPatterns(par, m);
		runExperimentsOne(par, m, kst);

		m=10;
		createPatterns(par, m);
		runExperimentsOne(par, m, kst);

		if (par->index->nDocs > 100){
			kst=100;
			m=6;
			createPatterns(par, m);
			runExperimentsOne(par, m, kst);
			m=10;
			createPatterns(par, m);
			runExperimentsOne(par, m, kst);
		}

		delete [] (par->patterns);
	}
}

