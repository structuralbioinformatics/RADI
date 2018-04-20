#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include"./svdcmp.h"

/*S'ha d'optimitzar el el calcul de pseudocount sense fer servir lambda_psc, que després faig que sigui igual al Meff! Vigilar en el codi final!!*/

void MaxEnt (int, int, double**, double****, double**);
void MI_Compute (int, int, double**, double****, double**);
void Pseudoinverse (int, double**);

/*
q = number of possible positions for each position in the alignment
L = length of the MSA
numseqs = number of homolog sequences in the MSA
*/

/*ARGV[]:
0.- Program name.
1.- MSA input file
2.- DI output file
3.- -f or -c for clustered or full DCA algorithm (obviable and by default clustered?)
(4.- -ss (for loading a secondary structure file type))
*/


int main(int argc, char *argv[]){
	int i, j, k, A, B, id, aux, aux2, L, l, linel, lold, q, numseqs, cluster;
	int *pos;
	double Max_id, Max_lid, Max_gaps, Meff, lambda_psc, lambda_Meff, daux;
	double *gaps, *seqs, **DI, **MI, **f, ****ff;
	char aux_aa;
	char *aa, *mini_aa, *line, *ss, *MSA_filename, *DI_filename, *SS_filename;
	char **MSA;
	clock_t t;
	FILE *MSA_File, *DI_File, *SS_File;

	Max_id = 0.8;
	lambda_psc = 0.5;
	Max_gaps = 0.15;
	linel = 131072;

//By default, the algorithm runs with its clustered version	
	q = 9;
	cluster = 1;
	SS_filename = "-";

	//We read the execution parameters
	if (argc == 3){
		MSA_filename = argv[1];
		DI_filename = argv[2];
	}
	else if (argc == 4){
		MSA_filename = argv[1];
		DI_filename = argv[2];
		if (strcmp(argv[3],"-f") == 0){
			q = 21;
			cluster = 0;
		}
		else if (strcmp(argv[3],"-c") == 0){
			q = 9;
			cluster = 1;
		}
		else{
			printf("Given arguments are invalid! Better contact with Stephen Hawking!\n");
			exit(0);
		}
	}
	else if (argc == 5){
		MSA_filename = argv[1];
		DI_filename = argv[2];
		if (strcmp(argv[3],"-ss") == 0){
			SS_filename = argv[4];
		}
		else{
			printf("Given arguments are invalid! Better contact with Stephen Hawking!\n");
			exit(0);
		}
	}
	else if (argc == 6){
		MSA_filename = argv[1];
		DI_filename = argv[2];
		if (strcmp(argv[3],"-f") == 0){
			q = 21;
			cluster = 0;
			if (strcmp(argv[4],"-ss") == 0){
				SS_filename = argv[5];
			}
			else{
				printf("Given arguments are invalid! Better contact with Stephen Hawking!\n");
				exit(0);
			}
		}
		else if (strcmp(argv[3],"-c") == 0){
			q = 9;
			cluster = 1;
			if (strcmp(argv[4],"-ss") == 0){
				SS_filename = argv[5];
			}
			else{
				printf("Given arguments are invalid! Better contact with Stephen Hawking!\n");
				exit(0);
			}
		}
		else if (strcmp(argv[3],"-ss") == 0){
			SS_filename = argv[4];
			if (strcmp(argv[5],"-f") == 0){
				q = 21;
				cluster = 0;
			}
			else if (strcmp(argv[5],"-c") == 0){
				q = 9;
				cluster = 1;
			}
			else{
				printf("Given arguments are invalid! Better contact with Stephen Hawking!\n");
				exit(0);
			}
		}
		else{
			printf("Given arguments are invalid! Better contact with Stephen Hawking!\n");
			exit(0);
		}

	}
	else{
		printf("Given arguments are invalid! Better contact with Stephen Hawking!\n");
		exit(0);
	}

	//We define the aa dictionary
	mini_aa = ".rhkdestnqavlimfwycgp";
	if (cluster == 0){
		aa = "-RHKDESTNQAVLIMFWYCGP";
	}
	else if (cluster == 1){
		aa = "012345678";
	}

	//Here we infer L and numseqs from the MSA
	MSA_File = fopen(MSA_filename, "r");
	if (MSA_File == NULL){
		printf("Problems opening the MSA file\n");
		exit(2);
	}
	numseqs = 0;
	L = 0;
	line = (char*) malloc(linel*sizeof(char));
	while (fgets(line,10000,MSA_File) != NULL){
		if (line[0] == '>'){
			numseqs++;
			continue;
		}
		if (numseqs == 1){
			i = 0;
			while (line[i] != '\n'){
				i++;
			}
			L = L+i;
		}
	}
	fclose(MSA_File);

	if (L==0 || numseqs == 0){
		printf("MSA file is empty: nothing to do here!\n");
		exit(4);
	}

	//Ara reservem memoria i llegim el MSA:
	MSA = (char**) malloc(numseqs*sizeof(char*));
	if (MSA == NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (i=0; i<numseqs; i++){
		MSA[i] = (char*) malloc(L*sizeof(char));
		if (MSA[i] == NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}
	MSA_File = fopen(MSA_filename, "r");
	if (MSA_File == NULL){
		printf("Problems opening the MSA file\n");
		exit(2);
	}

	k = -1;
	aux = 0;
	while (fgets(line,linel,MSA_File) != NULL){
		if (line[0] == '>'){
			aux = 1;
			i = 0;
			k++;
			continue;
		}
		j = 0;
		if (aux > 0){
			while (line[j] != '\n'){
				if (cluster == 0){ //Control for non-std aa
					aux = 2;
					for (A=0; A<q; A++){
						if (line[j] == aa[A]){
							aux = 1;
							break;
						}
						if (line[j] == mini_aa[A]){
							aux = 1;
							break;
						}
					}
					if (aux == 2){
						line[j] = '-'; //Convertim els non-std aa en gaps
//						printf("Non-standard amino acids were found in the MSA: %c\n", line[j]);
//						exit (3);
					}
				}
				else if (cluster == 1){
					aux = 2;
					for (A=0; A<21; A++){
						if (line[j] == mini_aa[A]){
							aux = 1;
							break;
						}
					}
					//Clustering process
					if (aux == 2){
						if (line[j] == '-'){ //Gap
							line[j] = '0';
							aux = 1;
						}
						else if (line[j] == 'R' || line[j] == 'H' || line[j] == 'K'){ //Positive charge
							line[j] = '1';
							aux = 1;
						}
						else if (line[j] == 'D' || line[j] == 'E'){ //Negative charge
							line[j] = '2';
							aux = 1;
						}
						else if (line[j] == 'S' || line[j] == 'T' || line[j] == 'N' || line[j] == 'Q'){ //Polar
							line[j] = '3';
							aux = 1;
						}
						else if (line[j] == 'A' || line[j] == 'V' || line[j] ==	'L' || line[j] == 'I' || line[j] == 'M'){ //Non-polar (aliphatic)
							line[j] = '4';
							aux = 1;
						}
						else if (line[j] == 'F' || line[j] == 'W' || line[j] == 'Y'){ //Aromatic rings
							line[j] = '5';
							aux = 1;
						}
						else if (line[j] == 'C'){ //Disulphide bonds
							line[j] = '6';
							aux = 1;
						}
						else if (line[j] == 'G'){ //Little
							line[j] = '7';
							aux = 1;
						}
						else if (line[j] == 'P'){ //Turn after beta-sheet
							line[j] = '8';
							aux = 1;
						}
					}
					if (aux == 2){ //Control for non-std aa
						line[j] = '0'; //Convertim els non-std aa en gaps
//						printf("Non-standard amino acids were found in the MSA: %c\n", line[j]);
//						exit (3);
					}
				}
				MSA[k][i] = (char) line[j];
				i++;
				j++;
			}
		}
	}
	free(line);
	fclose(MSA_File);

//We check if the alignment has been well read
//for (k=0; k<numseqs; k++){ for (i=0; i<L; i++){ printf("%c",MSA[k][i]); } printf("\n"); } printf("\n");

	gaps = (double*) calloc(L,sizeof(double));
	seqs = (double*) calloc(numseqs,sizeof(double));
	pos = (int*) malloc(L*sizeof(int));
	for (i=0; i<L; i++){
		pos[i] = i;
	}

	/*Those columns with lowercase are not well aligned, such as those with gap in the seed sequence, we erase them (and keep record on pos)*/
	j = 0;
	for (i=0; i<L; i++){
		aux = 0;
		if (MSA[0][i] == aa[0]){
			continue;
		}
		for (A=0; A<q; A++){
			if (MSA[0][i] == mini_aa[A]){
				aux = 1;
				break;
			}
		}
		if (aux == 0){
			pos[j] = i;
			j++;
		}
	}
	l = j;
/*for(k=0; k<numseqs; k++){
for(i=0; i<l; i++){
printf("%c",MSA[k][pos[i]]);
}
printf("\n");
}*/
	/*Here we only keep these columns that are well aligned and with no more than a Max_gaps frequency of gaps*/
	lold = 0;
	while (l != lold){
		Max_lid = Max_id*l;
		Meff = 0;
		for (i=0; i<numseqs; i++){
			seqs[i] = 0;
			for (j=0; j<numseqs; j++){
				id = 0;
				for (k=0; k<l; k++){
					if (MSA[i][pos[k]] == MSA[j][pos[k]]){
						id++;
					}
				}
				if (id >= Max_lid){
					seqs[i]++;
				}
			}
			seqs[i] = 1./seqs[i];
			Meff = Meff + seqs[i];
		}
	lambda_psc = Meff; //Best value for lambda_psc
		lambda_Meff = lambda_psc + Meff;
		daux = lambda_psc/q;
		for (i=0; i<l; i++){
			gaps[i] = daux;
			for (k=0; k<numseqs; k++){
				if (MSA[k][pos[i]] == aa[0]){
					gaps[i] = gaps[i] + seqs[k];
				}
			}
			gaps[i] = gaps[i]/lambda_Meff;
		}

		j = 0;
		for (i=0; i<l; i++){
			aux = 1;
//printf("%d->%g\n", pos[i], gaps[i]);
			if (gaps[i] > Max_gaps){
//printf("Erased (gaps)\n");
				aux = 0;
				continue;
			}
			if (aux == 1){
				aux = pos[i];
				pos[j] = aux;
				j++;
			}
		}
		lold = l;
		l = j;
	}
	free(gaps);

	// Reserve memory for f, ff, and count pseudofrequencies
	f = (double**) malloc(l*sizeof(double*));
	ff = (double****) malloc(l*sizeof(double***)); 
	if (f == NULL || ff == NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (i=0; i<l; i++){
		f[i] = (double*) calloc(q,sizeof(double));
		ff[i] = (double***) malloc(l*sizeof(double**));
		if (f[i] == NULL || ff[i] == NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
		for (j=0; j<l; j++){
			ff[i][j] = (double**) malloc(q*sizeof(double*));
			if (ff[i][j] == NULL){
				printf("It was not possible to allocate memory\n");
				exit(1);
			}
			for (A=0; A<q; A++){
				ff[i][j][A] = (double*) calloc(q,sizeof(double));
				if (ff[i][j][A] == NULL){
					printf("It was not possible to allocate memory\n");
					exit(1);
				}
			}				
		}
	}

	Meff = 0;
	Max_lid = Max_id*l;
	for (i=0; i<numseqs; i++){
		seqs[i] = 0;
		for (j=0; j<numseqs; j++){
			id = 0;
			for (k=0; k<l; k++){
				if (MSA[i][pos[k]] == MSA[j][pos[k]]){
					id++;
				}
			}
			if (id >= Max_lid){
				seqs[i]++;
			}
		}
		seqs[i] = 1./seqs[i];
		Meff = Meff + seqs[i];
	}
printf("L=%d (l=%d) numseqs=%d Meff=%g q=%d\n", L, l, numseqs, Meff, q);

	lambda_psc = Meff; //Best value for lambda_psc
	lambda_Meff = lambda_psc+Meff;

	lambda_psc = lambda_psc/q;
	for (i=0; i<l; i++){
		for (A=0; A<q; A++){
			f[i][A] = lambda_psc;
			ff[i][i][A][A] = lambda_psc; //això es per comptar diferent els elements de la diagonal... pq? nobody knows, però ho fan així
			for (k=0; k<numseqs; k++){
				if (MSA[k][pos[i]] == aa[A]){
					f[i][A] = f[i][A] + seqs[k];
				}
			}
			f[i][A] = f[i][A]/lambda_Meff;
		}
	}

	lambda_psc = lambda_psc/q;
	for (i=0; i<l; i++){
		for (j=0; j<l; j++){
			for (A=0; A<q; A++){
				for (B=0; B<q; B++){
					if (i != j){ //això es per comptar diferent els elements de la diagonal... pq? nobody knows, però ho fan així
						ff[i][j][A][B] = lambda_psc;
					}
					for (k=0; k<numseqs; k++){
						if((MSA[k][pos[i]] == aa[A]) && (MSA[k][pos[j]] == aa[B])){
							ff[i][j][A][B] = ff[i][j][A][B] + seqs[k];
						}
					}
					ff[i][j][A][B] = ff[i][j][A][B]/lambda_Meff;
				}
			}
		}
	}
	free(seqs);

	for (i=0; i<numseqs; i++){
		free(MSA[i]);
	}
	free(MSA);

	/* Reserve memory for MI and DI */
	DI = (double**) malloc(l*sizeof(double*));
	MI = (double**) malloc(l*sizeof(double*));
	if (DI == NULL || MI == NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (i=0; i<l; i++){
		DI[i] = (double*) calloc(l,sizeof(double));
		MI[i] = (double*) calloc(l,sizeof(double));
		if (DI[i] == NULL || MI[i] == NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}

	MI_Compute(l, q, f, ff, MI);

	t = clock();

printf("Pre-MaxEnt at time t=%d\n", ((float)t)/CLOCKS_PER_SEC);

	MaxEnt (l, q, f, ff, DI);

t = clock()-t;
t = ((float)t)/CLOCKS_PER_SEC;

printf("Post-MaxEnt after time t=%d\n", t);

	line = (char*) malloc(linel*sizeof(char));
	ss = (char*) malloc(linel*sizeof(char));

/*SS_filename = "-";
	if (SS_filename != "-"){
		SS_File = fopen(SS_filename, "r");
		if (SS_File == NULL){
			printf("Problems opening the secondary structure file\n");
			exit(2);
		}
		aux = 0;
		while (fgets(line,linel,MSA_File) != NULL){
			if (line[0] == '>'){
				aux = 1;
				continue;
			}
			//LLEGIR SS I GUARDAR LES POSICIONS ADIENTS AL *ss
		}
		//AQUI, A MESURA QUE AVANCEM LLEGINT EL *ss, MIRO CANVIS EN ESTRUCTURES SECUNDARIES. PER TOTS ELS PARELLS DINS LA MATEIXA ESTRUCTURA FORCEM DI=0
				


	}*/
	free(line);
	free(ss);

//RETOCAR LA FORMA DE GUARDAR LA DI I INCORPORAR LES POSICIONS CORRECTES A L'ALINEAMENT! FER-HO EN PLAN COM EL MORCOS.
	/* Copy into a file the matrix of Direct Information */
	DI_File = fopen(DI_filename, "w");
	if (DI_File == NULL){
		printf("Problems opening the writing DI file\n");
		exit(2);
	}
	fprintf(DI_File, ">>L=%d, l=%d, numseqs=%d, Meff=%.2f\n", L, l, numseqs, Meff);
	fprintf(DI_File, ">>Time: %d\n", t);
	fprintf(DI_File, ">>Positions: ");
	for (i=0; i<l; i++){
		fprintf(DI_File, "%d,", pos[i]);
	}
	fprintf(DI_File, "\n\n");
	fprintf(DI_File, "#Pos\tMI\t\tDI\n");
	for (i=0; i<l; i++){
		for (j=i+1; j<l; j++){
			fprintf(DI_File, "%d %d\t%.6f\t%.6f\n", pos[i], pos[j], MI[i][j], DI[i][j]);
		}
	}
	free(pos);
	fclose(DI_File);

	/* Free memory for DI, f and ff */
	for (i=0; i<l; i++){
		for (j=0; j<l; j++){
			for (A=0; A<q; A++){
				free(ff[i][j][A]);
			}
			free(ff[i][j]);
		}
		free(DI[i]);
		free(MI[i]);
		free(f[i]);
		free(ff[i]);
	}
	free(DI);
	free(MI);
	free(f);
	free(ff);

	return 0;
}





void MaxEnt (int L, int q, double **f, double ****ff, double **DI){

	int i, j, k, l, q1, Lq, A, B, A1, B1, count, Maxcount, auxexit;
	double c, daux, mu0, sum, sum1, sum2, dist, eps, Z;
	double *mu1, *mu2, *scra1, *scra2, *new1, *new2;
	double **C, **h, **h_new, **Pdir;
	double ****e;

	Maxcount = 1e8;
	eps = 1e-8;

	q1 = q-1;
	Lq = L*q1;

	C = (double**) malloc(Lq*sizeof(double*));
	e = (double****) malloc(L*sizeof(double***));
	if (C == NULL || e == NULL){
		printf("It was not possible to allocate memory \n");
		exit(1);
	}
	for (i=0; i<L; i++){
		e[i] = (double***) malloc(L*sizeof(double**));
		if (e[i] == NULL){
			printf("It was not possible to allocate memory (e)\n");
			exit(1);
		}
		for (j=0; j<L; j++){
			e[i][j] = (double**) malloc(q*sizeof(double*));
			if (e[i][j] == NULL){
				printf("It was not possible to allocate memory (e)\n");
				exit(1);
			}
			for (A=0; A<q; A++){
				e[i][j][A] = (double*) calloc(q,sizeof(double));
				if (e[i][j][A] == NULL){
					printf("It was not possible to allocate memory (e)\n");
					exit(1);
				}
			}
		}
	}
	for(i=0; i<Lq; i++){
		C[i] = (double*) calloc(Lq,sizeof(double));
		if (C[i] == NULL){
			printf("It was not possible to allocate memory (C)\n");
			exit(1);
		}
	}

	/* Gauge invariance method for obtaining interaction lagrange multipliers e */
	for (i=0; i<L; i++){
		for (j=0; j<L; j++){
			for (A=1; A<q; A++){
				for (B=1; B<q; B++){
					A1 = A-1;
					B1 = B-1;
					C[i*q1+A1][j*q1+B1] = ff[i][j][A][B]-f[i][A]*f[j][B];
				}
			}
		}
	}

	Pseudoinverse(Lq, C);

	for (i=0; i<Lq; i++){
		for (j=0; j<Lq; j++){
			k = i/q1;
			A = i%q1+1;
			l = j/q1;
			B = j%q1+1;
			e[k][l][A][B] = exp(-C[i][j]);
		}
	}
	for (i=0; i<L; i++){
		for (j=0; j<L; j++){
			e[i][j][0][0] = 1.;
			for (A=0; A<q; A++){
				e[i][j][A][0] = 1.;
				e[i][j][0][A] = 1.;
			}
		}
	}

	for (i=0; i<Lq; i++){
		free(C[i]);
	}
	free(C);

	h = (double**) malloc(L*sizeof(double*));
	h_new = (double**) malloc(L*sizeof(double*));
	if (h == NULL || h_new ==NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (i=0; i<L; i++){
		h[i] = (double*) malloc(q*sizeof(double));
		h_new[i] = (double*) malloc(q*sizeof(double));
		if (h[i] == NULL || h_new[i] ==  NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}
	Pdir = (double**) malloc(q*sizeof(double*));
	if (Pdir == NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (A=0; A<q; A++){
		Pdir[A] = (double*) calloc(q,sizeof(double));
		if (Pdir == NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}
	mu1 = (double*) malloc(q*sizeof(double));
	mu2 = (double*) malloc(q*sizeof(double));
	scra1 = (double*) malloc(q*sizeof(double));
	scra2 = (double*) malloc(q*sizeof(double));
	new1 = (double*) malloc(q*sizeof(double));
	new2 = (double*) malloc(q*sizeof(double));
	mu0 = exp(1./q);
//mu0 = 1./q; //Triga les mateixes iteracions comenci on comenci!!!!???? L'anterior makes more sense...
	for (i=0; i<L-1; i++){
		DI[i][i] = 1.;
		for (j=i+1; j<L; j++){
			for (A=0; A<q; A++){
				mu1[A] = mu0;
				mu2[A] = mu0;
			}
			for (count=0; count<Maxcount; count++){
				dist = 0;
				for (A=0; A<q; A++){
					scra1[A] = 0;
					scra2[A] = 0;
					for (B=0; B<q; B++){
						scra1[A] = scra1[A] + e[i][j][A][B]*mu2[B];
						scra2[A] = scra2[A] + e[i][j][B][A]*mu1[B];
					}
				}
				sum1 = 0;
				sum2 = 0;
				for (A=0; A<q; A++){
					new1[A] = f[i][A]/scra1[A];
					new2[A] = f[j][A]/scra2[A];
					sum1 = sum1 + new1[A];
					sum2 = sum2 + new2[A];
				}
				for (A=0; A<q; A++){
					new1[A] = new1[A]/sum1;
					daux = fabs(new1[A]-mu1[A]);
					if (daux > dist){
						dist = daux;
					}
					new2[A] = new2[A]/sum2;
					daux = fabs(new2[A]-mu2[A]);
					if (daux > dist){
						dist = daux;
					}
					mu1[A] = new1[A];
					mu2[A] = new2[A];
				}
				if (dist < eps){
					break;
				}
			}
//printf("count[%d][%d] = %d\t", i,j, count);
			if (count == Maxcount){
				printf("Lagrange single-site marginals could not be obtained in %d iterations\n", Maxcount);
				exit(9);
			}
			Z = 0.;		
			for (A=0; A<q; A++){
				for (B=0; B<q; B++){
					Z = Z + e[i][j][A][B]*mu1[A]*mu2[B];
				}
			}
			for (A=0; A<q; A++){
				for (B=0; B<q; B++){
					Pdir[A][B] = (e[i][j][A][B]*mu1[A]*mu2[B])/Z;
				}
			}
			for(A=0; A<q; A++){
				for (B=0; B<q; B++){
					daux = Pdir[A][B]/(f[i][A]*f[j][B]);
					if (daux < 1e-100){
						continue;
					}
					DI[i][j] = DI[i][j] + Pdir[A][B]*log(daux);
				}
			}
			DI[j][i] = DI[i][j];
//printf("DI: %.2g\n",DI[i][j]);
		}
	}
	free(mu1);
	free(mu2);
	free(scra1);
	free(scra2);
	free(new1);
	free(new2);

	for (i=0; i<L; i++){
		for (j=0; j<L; j++){
			for (A=0; A<q; A++){
				free(e[i][j][A]);
			}
			free(e[i][j]);
		}
		free(e[i]);
		free(h[i]);
		free(h_new[i]);
	}
	free(e);
	free(h);
	free(h_new);

	/* Free memory for Pdir */
	for (A=0; A<q; A++){
		free(Pdir[A]);
	}
	free(Pdir);

	return;	
}

void MI_Compute (int L, int q, double **f, double ****ff, double **MI){
	int i, j, A, B;
	double daux;
	
	for (i=0; i<L; i++){
		MI[i][i] = 0;
		for (j=i+1; j<L; j++){
			for (A=0; A<q; A++){
				for (B=0; B<q; B++){
					daux = ff[i][j][A][B]/(f[i][A]*f[j][B]);
					if (daux < 1e-10){
						continue;
					}
					MI[i][j] = MI[i][j] + ff[i][j][A][B]*log(daux);
				}
			}
			MI[j][i] = MI[i][j];
		}
	}	

	return;
}

/* Routine to compute Moore-Penrose pseudoinverse of A using SVD procedure, and store it into A */
void Pseudoinverse (int n, double **A){
	int i, j, k;
	double E, aux;
	double *D, **V, **AUX;

	E = 1e-3;

	D = (double*) calloc(n,sizeof(double)); /* Vector for the diagonal matrix of SVD descomposition */
	AUX = (double**) malloc(n*sizeof(double*));
	V = (double**) malloc(n*sizeof(double*));
	if (D == NULL || AUX == NULL || V == NULL){
		printf("It was not possible to allocate memory\n");
		exit(1);
	}
	for (i=0; i<n; i++){
		AUX[i] = (double*) calloc(n,sizeof(double));
		V[i] = (double*) calloc(n,sizeof(double));
		if (V[i] == NULL || AUX == NULL){
			printf("It was not possible to allocate memory\n");
			exit(1);
		}
	}
	/* Let's perform a SVD decomposition A=UDV^t, where U will be stored into A */
	/* Ojo! V is not transposed in the result! */
	svdcmp(A,n,n,D,V);
	/* Compute the pseudoinverse of the diagonal matrix D */
	for (i=0; i<n; i++){
		if (fabs(D[i])>E){
			D[i] = 1./D[i];
		}
		else{
			D[i] = 0.;
		}
	}

	/* Undo the basis change for SVD descomposition, with the pseudoinverse of D, in order to obtain the pseudoinverse of A */
	/* Remember that we have to transpose V in order to undo the basis change! And that D is a diagonal matrix! */
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			A[i][j] = A[i][j]*D[j];
		}
	}
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			AUX[i][j] = 0.;
			for (k=0; k<n; k++){
				AUX[i][j] = AUX[i][j] + A[i][k]*V[j][k];
			}
		}
	}
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			A[i][j] = AUX[i][j];
		}
	}

	/* Free memory for V and D */
	for (i=0; i<n; i++){
		free(AUX[i]);
		free(V[i]);
	}
	free(AUX);
	free(V);
	free(D);

	return;
}
