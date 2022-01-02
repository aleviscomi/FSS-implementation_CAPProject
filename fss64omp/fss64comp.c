/**************************************************************************************
*
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
*
* Progetto dell'algoritmo Fish School Search 221 231 a
* in linguaggio assembly x86-64 + SSE
*
* Fabrizio Angiulli, aprile 2019
*
**************************************************************************************/

/*
*
* Software necessario per l'esecuzione:
*
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
*
* entrambi sono disponibili come pacchetti software
* installabili mediante il packaging tool del sistema
* operativo; per esempio, su Ubuntu, mediante i comandi:
*
*    sudo apt-get install nasm
*    sudo apt-get install gcc
*
* potrebbe essere necessario installare le seguenti librerie:
*
*    sudo apt-get install lib64gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
*
* Per generare il file eseguibile:
*
* nasm -f elf64 fss64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o fss64.o fss64c.c -o fss64c -lm && ./fss64c $pars
*
* oppure
*
* ./runfss64
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <stdbool.h>
#include <omp.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
	MATRIX x; //posizione dei pesci
	VECTOR xh; //punto associato al minimo di f, soluzione del problema
	VECTOR c; //coefficienti della funzione
	VECTOR r; //numeri casuali
	int np; //numero di pesci, quadrato del parametro np
	int d; //numero di dimensioni del data set
	int iter; //numero di iterazioni
	type stepind; //parametro stepind
	type stepvol; //parametro stepvol
	type wscale; //parametro wscale
	int display;
	int silent;
} params;

/*
*
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate
* 	mediante un array (double*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere
* 	memorizzate mediante array di array (double**).
*
* 	In entrambi i casi il candidato dovr? inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente ? che le matrici siano in row-major order.
*
*/

void* get_block(int size, int elements) {
	return _mm_malloc(elements*size,32);
}

void free_block(void* p) {
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
	free_block(mat);
}

/*
*
* 	load_data
* 	=========
*
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
*
* 	Codifica del file:
* 	primi 4 byte: numero di colonne (N) --> numero intero
* 	successivi 4 byte: numero di righe (M) --> numero intero
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri doubleing-point a precisione singola
*
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice.
*****************************************************************************
*
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;

	fp = fopen(filename, "rb");

	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}

	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);

	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);

	*n = rows;
	*k = cols;

	return data;
}

/*
* 	save_data
* 	=========
*
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
*
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 64 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 64 bit
* 	successivi N*M*4 byte: matrix data in row-major order --> numeri interi o doubleing-point a precisione singola
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

// PROCEDURE ASSEMBLY

extern void prodottoScalare(VECTOR x, VECTOR y, int size, type* ps);
extern void distanzaEuclidea(type* v1, type* v2, int dim, type* dist);
extern void mediaPesata(params* input, MATRIX matrix, VECTOR vect, VECTOR mediaP, type sumVect);
extern void sommaVettoreMatrice(params* input, MATRIX matrix, VECTOR vect);
extern void generaPossibileMovimento(params* input, int* randIndex, int i, VECTOR y);
extern void muoviPesce(MATRIX x, VECTOR y, MATRIX dx, int i, int d);
extern void mantieniPosizionePesce(MATRIX dx, int i, int d);
extern void faiMovimentoVolitivo(params* input, VECTOR b, type dist, type randNum, bool weightGain, int i);
extern void sommaElementiVettore(params* input, VECTOR df, type* sumdf);


// FUNZIONI DI UTILITA'

/*
 * Funzione che calcola un numero getRandom tra 'from' e 'to' a partire
 * dal numero 'input->r[randIndex]
 */
type getRandom(params* input, int* randIndex, int from, int to) {
    return from + (to-from)*input->r[(*randIndex)++];
    //return from + (to-from)*((float) rand() / (float) RAND_MAX); //calcolo dei getRandom alternativo utile per fare test
} //getRandom

/*
 * Funzione che calcola e restituisce, dati i vettori c e x, e la dimensione size di
 * questi due, il valore di funzione:
 * f(x) = e^(x^2) + x^2 - c o x
 */
type function(VECTOR c, VECTOR x, int size) {
    int i;
    type xQuad;
    type cx;
	
    //calcolo x^2 e c o x (in C)
    xQuad = 0;
	cx = 0;
	for(i=0; i < size; i++) {
        xQuad += x[i] * x[i];
		cx += c[i] * x[i];
	}

    return exp(xQuad) + xQuad - cx;
} //function

/*
 * Funzione che calcola (x) = e^(x^2) + x^2 - c o x, con le ottimizzazioni Assembly
 */
type functionAssembly(VECTOR c, VECTOR x, int size) {
    int i;
    type xQuad;
    type cx;
    
    //calcolo x^2 e c o x (in assembly)
    xQuad = 0;
	cx = 0;
    prodottoScalare(x, x, size, &xQuad);	
    prodottoScalare(c, x, size, &cx);

    return exp(xQuad) + xQuad - cx;
} //functionAssembly


/*
 * Funzione che calcola e restituisce, dato un vettore v e la sua dimensione size,
 * il minimo di v
 */
type min(VECTOR v, int size){
    int i;
    type min;

    min = v[0];
    for(i=0; i<size; i++)
        if(v[i] < min)
            min = v[i];
    return min;
} //min




// FUNZIONI DELL'ALGORITMO FSS

/*
 * Funzione che calcola il movimento individuale del pesce i-esimo.
 * Restituisce in df e in dx, rispettivamente, la variazione di funzione
 * e di posizione del pesce i-esimo
 */
void individualMovement(params* input, int* randIndex, int i, VECTOR df, MATRIX dx) {
    int j;
    VECTOR y; // y[i], vettore che contiene la possibile nuova posizione del pesce i-esimo
    type fx; //f(x[i])
    type fy; //f(y[i])

    y = (VECTOR) alloc_matrix(1, input->d);

    //calcolo la possibile variazione di posizione del pesce i-esimo, come: y[j] = x[i][j] + rand(-1, 1) * stepind
    for(j=0; j<input->d; j++)
		y[j] = input->x[input->d * i + j] + getRandom(input, randIndex, -1, 1) * input->stepind;
	
    //calcolo i valori di funzione con x e y
    fx = function(input->c, &input->x[input->d * i], input->d);
    fy = function(input->c, y, input->d);
	
	//se la nuova posizione è migliore allora mi sposto, altrimenti mantengo la precedente
    if(fy < fx) {
        df[i] = fy - fx;
        for(j=0; j<input->d; j++) {
            dx[input->d * i + j] = y[j] - input->x[input->d * i + j];
            input->x[input->d * i + j] = y[j];
        }
    }
    else {
        df[i] = 0;
        for (j = 0; j < input->d; j++)
            dx[input->d * i + j] = 0;
    }

    dealloc_matrix(y);
} //individualMovement

/*
 * Funzione che calcola il movimento individuale con l'aggiunta delle ottimizzazioni Assembly
 */
void individualMovementAssembly(params* input, int* randIndex, int i, VECTOR df, MATRIX dx) {
    int j;
    VECTOR y; // y[i], vettore che contiene la possibile nuova posizione del pesce i-esimo
    type fx; //f(x[i])
    type fy; //f(y[i])

    y = (VECTOR) alloc_matrix(1, input->d);
	
    //calcolo la possibile variazione di posizione del pesce i-esimo, come: y[j] = x[i][j] + rand(-1, 1) * stepind
	generaPossibileMovimento(input, randIndex, i, y);
	
    //calcolo i valori di funzione con x e y
    fx = functionAssembly(input->c, &input->x[input->d * i], input->d);
    fy = functionAssembly(input->c, y, input->d);

	
    //se la nuova posizione è migliore allora mi sposto, altrimenti mantengo la precedente
    if(fy < fx) {
        df[i] = fy - fx;
		muoviPesce(input->x, y, dx, i, input->d);
    }
    else {
        df[i] = 0;
		mantieniPosizionePesce(dx, i, input->d);
    }

    dealloc_matrix(y);
} //individualMovementAssembly


/*
 * Funzione che aggiorna il peso di tutti i pesci nel caso in cui
 * questi si siano spostati in una posizione migliore (df<0).
 *
 * Se almeno un pesce aggiornerà  il proprio peso, verrà  settata
 * a true la variabile booleana 'weightGain', per indicare che
 * il peso del banco è aumentato.
 * Ciò servirà  per il movimento volitivo, per indicare se il banco
 * debba contrarsi verso il baricentro, o meno
 */
void feedOperator(params* input, VECTOR w, VECTOR df, bool* weightGain){
    int i;
    type dfmin; //min{df}

    dfmin = min(df, input->np);

    for(i=0; i < input->np; i++)
        if(df[i] < 0) { //se il pesce i-esimo si avvicina al minimo, si è spostato e può mangiare
            w[i] += df[i] / dfmin;
            *weightGain = true;
        }

} //feedOperator


/*
 * Funzione che permette a tutti i pesci di essere attratti verso
 * i pesci che hanno avuto un miglioramento maggiore.
 * Se nessun pesce si è mosso (sumdf=0) i pesci rimangono tutti
 * fermi.
 */
void instinctiveMovement(params* input, MATRIX dx, VECTOR df) {
    int i,j;
    VECTOR vectI; //vettore movimento istintivo (I)
    type sumdf; //sommatoria df[i]

    //inizializzazione variabili
    vectI = (VECTOR) alloc_matrix(1, input->d);
    for(i=0; i<input->d; i++)
        vectI[i] = 0;
	
	//SUM(df)
    sumdf=0;
    for(i=0; i < input->np; i++)
        sumdf+=df[i];

    //calcolo I
 	if(sumdf != 0)
		for(i=0; i < input->np; i++)
			for(j=0; j < input->d; j++)
				vectI[j] += dx[input->d * i + j] * df[i] / sumdf;
	
    //sommo I ai vettori x[i]
    for(i=0; i < input->np; i++)
        for(j=0; j < input->d; j++)
            input->x[input->d * i + j] += vectI[j];
	
	
    dealloc_matrix(vectI);
} //instinctiveMovement

/*
 * Funzione che esegue il movimento istintivo con le ottimizzazioni Assembly
 */
void instinctiveMovementAssembly(params* input, MATRIX dx, VECTOR df) {
    int i,j;
    VECTOR vectI; //vettore movimento istintivo (I)
    type sumdf; //sommatoria df[i]

    //inizializzazione variabili
    vectI = (VECTOR) alloc_matrix(1, input->d);
    for(i=0; i<input->d; i++)
        vectI[i] = 0;
	
	//SUM(df)
    sommaElementiVettore(input, df, &sumdf);

	//calcolo I
	if(sumdf != 0)
		mediaPesata(input, dx, df, vectI, sumdf);
	
    //sommo I ai vettori x[i]
    sommaVettoreMatrice(input, input->x, vectI);
	
	
    dealloc_matrix(vectI);
} //instinctiveMovementAssembly


/*
 * Funzione che calcola il baricentro dei pesci
 */
void calculateB(params* input, VECTOR w, VECTOR b) {
    int i,j;
    type sumW; //sommatoria W[i]

    //inizializzazione variabili
    for(i=0; i<input->d; i++)
        b[i] = 0;
	
	//SUM(w)
    sumW=0;
    for(i=0; i < input->np; i++)
        sumW+=w[i];

    //calcolo B
    for(i=0; i < input->np; i++) 
        for(j=0; j < input->d; j++)
            b[j] += input->x[input->d * i + j] * w[i] / sumW;
    
} //calculateB

/*
 * Funzione che calcola il baricentro dei pesci con le ottimizzazioni Assembly
 */
void calculateBAssembly(params* input, VECTOR w, VECTOR b) {
    int i,j;
    type sumW; //sommatoria W[i]

    //inizializzazione variabili
    for(i=0; i<input->d; i++)
        b[i] = 0;
	
	//SUM(w)
    sommaElementiVettore(input, w, &sumW);
	
	//calcolo B
    mediaPesata(input, input->x, w, b, sumW);
	
} //calculateBAssembly


/*
 * Funzione che calcola il movimento volitivo, facendo
 * contrarre il banco al baricentro se 'weightGain'=true
 * o allontanarlo se pari a false
 */
void volitiveMovement(params* input, int* randIndex, VECTOR b, bool weightGain){
    //dichiarazione
    int i, j;
    VECTOR numerator;
    type dist;
    type randNum;

    //inizializzazione
    numerator = (VECTOR) alloc_matrix(1,input->d);
	
    //calcolo, per ogni pesce, la nuova posizione dovuta al movimento volitivo
    for(i=0; i<input->np; i++){
        //calcolo di un numero getRandom per ogni pesce
        randNum = getRandom(input, randIndex, 0, 1);

        dist = 0;
        //calcolo della distanza euclidea: sqrt(sum(((B[j]) - (x[i][j]))^2))
        for(j=0; j<input->d; j++)
            dist += (b[j] - input->x[i * input->d + j]) * (b[j] - input->x[i * input->d + j]);
        dist = sqrt(dist);

		//nel caso in cui il vettore Xi e B sono uguali (dunque distEuclidea=0) non sommo nulla ad Xi, per cui posso passare al movimento volitivo del prossimo pesce
		if(dist==0) continue;
		
        //calcolo del movimento per ogni coordinata del pesce i-esimo
        for(j=0; j<input->d; j++) {
            //(stepvol * rand(0,1) * (x[i][j] - B[j])) / (distanzaEuclidea)
            numerator[j] = (input->stepvol) * randNum * (input->x[i * input->d + j] - b[j]);
            numerator[j] /= dist;

            //se il peso del banco è aumentato allora i pesci si contraggono verso il baricentro (segno -), altrimenti viceversa (segno +)
            if (weightGain)
                input->x[i * input->d + j] -= numerator[j];
            else
                input->x[i * input->d + j] += numerator[j];
        }
    }

    dealloc_matrix(numerator);
} //volitiveMovement

/*
 * Funzione che calcola il movimento volitivo con le ottimizzazioni Assembly
 */
void volitiveMovementAssembly(params* input, int* randIndex, VECTOR b, bool weightGain){
    //dichiarazione
    int i, j;
    VECTOR numerator;
    type dist;
    type randNum;

    //inizializzazione
    numerator = (VECTOR) alloc_matrix(1,input->d);
	
    //calcolo, per ogni pesce, la nuova posizione dovuta al movimento volitivo
    for(i=0; i<input->np; i++){
        //calcolo di un numero getRandom per ogni pesce
        randNum = getRandom(input, randIndex, 0, 1);

        dist = 0;
		//calcolo della distanza euclidea: sqrt(sum(((B[j]) - (x[i][j]))^2))
        distanzaEuclidea(&input->x[i * input->d], b, input->d, &dist);

		//nel caso in cui il vettore Xi e B sono uguali (dunque distEuclidea=0) non sommo nulla ad Xi, per cui posso passare al movimento volitivo del prossimo pesce
		if(dist==0) continue;
	    
        //calcolo del movimento per ogni coordinata del pesce i-esimo
		faiMovimentoVolitivo(input, b, dist, randNum, weightGain, i);
    }

    dealloc_matrix(numerator);
} //volitiveMovementAssembly



// ALGORITMO FSS

void fss(params* input){
    /** dichiarazione */
    int it, i;
    int randIndex; //indice per il file contenente i getRandom
    VECTOR w; //vettore dei pesi dei pesci
    VECTOR df; //vettore delle variazioni di funzione dei pesci
    MATRIX dx; //matrice delle variazioni di posizione dei pesci
    VECTOR b;
    type stepindInit; //stepind iniziale
    type stepvolInit; //stepvol iniziale
    bool weightGain; //booleano che vale true se il peso totale del banco ï¿½ aumentato (cioï¿½ se almeno un pesce mangia)
    type min; //memorizza il min degli f(x) finali
    int minIndex; //memorizza l'indice del min degli f(x) finali

    double t;
    double t1=0;
    double t2=0;
    double t3=0;
    double t4=0;
    double t5=0;

    /** inizializzazione */
    // le posizioni dei pesci sono giï¿½ inizializzate leggendo la matrice x64_x_x.ds2
    w = (VECTOR) alloc_matrix(1, input->np);
    df = (VECTOR) alloc_matrix(1, input->np);
    dx = alloc_matrix(input->np, input->d);
    b = (VECTOR) alloc_matrix(1, input->d);
    stepindInit = input->stepind;
    stepvolInit = input->stepvol;
    randIndex = 0;
    //srand(time(NULL)); //utile in fase di testing per calcolare ad ogni esecuzione numeri getRandom differenti
    weightGain = false;

    //inizializzazione pesi pesci a wscale/2
    for (i = 0; i < input->np; i++)
        w[i] = input->wscale / 2;

    /** inizio fss */
    it = 0;
    while (it < input->iter) {
        /** movimento individuale */
        t = omp_get_wtime();
		#pragma omp parallel for
        for(i=0; i<input->np; i++)
            individualMovementAssembly(input, &randIndex, i, df, dx);
        t = omp_get_wtime() - t;
        t1 += t;

        /** operatore di alimentazione */
        t = omp_get_wtime();
        feedOperator(input, w, df, &weightGain);
        t = omp_get_wtime() - t;
        t2 += t;

        /** movimento istintivo */
        t = omp_get_wtime();
        instinctiveMovementAssembly(input, dx, df);
        t = omp_get_wtime() - t;
        t3 += t;

        /** calcolo baricentro */
        t = omp_get_wtime();
        calculateBAssembly(input, w, b);
        t = omp_get_wtime() - t;
        t4 += t;

        /** movimento volitivo */
        t = omp_get_wtime();
        volitiveMovementAssembly(input, &randIndex, b, weightGain);
        t = omp_get_wtime() - t;
        t5 += t;

        /** aggiorno parametri */
        input->stepind = input->stepind - stepindInit / input->iter;
        input->stepvol = input->stepvol - stepvolInit / input->iter;
        weightGain = false;
        it++;
    }
    printf("\nIndividual Movement: %f\n\n",t1);
    printf("\nFeed Operator: %f\n\n",t2);
    printf("\nInstinctive Movement: %f\n\n",t3);
    printf("\nBaricentre: %f\n\n",t4);
    printf("\nVolitive Movement: %f\n\n",t5);

    /** assegno ad xh, argmin(f) */
    min = function(input->c, &input->x[0], input->d);
    minIndex = 0;
    for (i = 1; i < input->np; i++) {
        type f = function(input->c, &input->x[input->d * i], input->d);
        if (f < min) {
            min = f;
            minIndex = i;
        }
    }
    input->xh = (VECTOR) alloc_matrix(1, input->d);
    for (i = 0; i < input->d; i++)
        input->xh[i] = input->x[input->d * minIndex + i];

    //COMMENTARE QUESTA RIGA!!!
    printf("f(xh) = %lf\n", min); //stampa di debug per visualizzare il valore ottimo trovato

    //deallocazione
    dealloc_matrix(w);
    dealloc_matrix(df);
    dealloc_matrix(dx);
    dealloc_matrix(b);
} //fss

int main(int argc, char** argv) {

	char fname[256];
	char* coefffilename = NULL;
	char* randfilename = NULL;
	char* xfilename = NULL;
	int i, j, k;
	double t;
	float time;

	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->x = NULL;
	input->xh = NULL;
	input->c = NULL;
	input->r = NULL;
	input->np = 25;
	input->d = 2;
	input->iter = 350;
	input->stepind = 1;
	input->stepvol = 0.1;
	input->wscale = 10;

	input->silent = 0;
	input->display = 0;

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 1){
		printf("%s -c <c> -r <r> -x <x> -np <np> -si <stepind> -sv <stepvol> -w <wscale> -it <itmax> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tc: il nome del file ds2 contenente i coefficienti\n");
		printf("\tr: il nome del file ds2 contenente i numeri casuali\n");
		printf("\tx: il nome del file ds2 contenente le posizioni iniziali dei pesci\n");
		printf("\tnp: il numero di pesci, default 25\n");
		printf("\tstepind: valore iniziale del parametro per il movimento individuale, default 1\n");
		printf("\tstepvol: valore iniziale del parametro per il movimento volitivo, default 0.1\n");
		printf("\twscale: valore iniziale del peso, default 10\n");
		printf("\titmax: numero di iterazioni, default 350\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-c") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing coefficient file name!\n");
				exit(1);
			}
			coefffilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-r") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing random numbers file name!\n");
				exit(1);
			}
			randfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-x") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing initial fish position file name!\n");
				exit(1);
			}
			xfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-np") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing np value!\n");
				exit(1);
			}
			input->np = atoi(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-si") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing stepind value!\n");
				exit(1);
			}
			input->stepind = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-sv") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing stepvol value!\n");
				exit(1);
			}
			input->stepvol = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-w") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing wscale value!\n");
				exit(1);
			}
			input->wscale = atof(argv[par]);
			par++;
		} else if (strcmp(argv[par],"-it") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing iter value!\n");
				exit(1);
			}
			input->iter = atoi(argv[par]);
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

	if(coefffilename == NULL || strlen(coefffilename) == 0){
		printf("Missing coefficient file name!\n");
		exit(1);
	}

	if(randfilename == NULL || strlen(randfilename) == 0){
		printf("Missing random numbers file name!\n");
		exit(1);
	}

	if(xfilename == NULL || strlen(xfilename) == 0){
		printf("Missing initial fish position file name!\n");
		exit(1);
	}

	int x,y;
	input->c = load_data(coefffilename, &input->d, &y);
	input->r = load_data(randfilename, &x, &y);
	input->x = load_data(xfilename, &x, &y);

	if(input->np < 0){
		printf("Invalid value of np parameter!\n");
		exit(1);
	}

	if(input->stepind < 0){
		printf("Invalid value of si parameter!\n");
		exit(1);
	}

	if(input->stepvol < 0){
		printf("Invalid value of sv parameter!\n");
		exit(1);
	}

	if(input->wscale < 0){
		printf("Invalid value of w parameter!\n");
		exit(1);
	}

	if(input->iter < 0){
		printf("Invalid value of it parameter!\n");
		exit(1);
	}

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Coefficient file name: '%s'\n", coefffilename);
		printf("Random numbers file name: '%s'\n", randfilename);
		printf("Initial fish position file name: '%s'\n", xfilename);
		printf("Dimensions: %d\n", input->d);
		printf("Number of fishes [np]: %d\n", input->np);
		printf("Individual step [si]: %f\n", input->stepind);
		printf("Volitive step [sv]: %f\n", input->stepvol);
		printf("Weight scale [w]: %f\n", input->wscale);
		printf("Number of iterations [it]: %d\n", input->iter);
	}

    // COMMENTARE QUESTA RIGA!
    // prova(input);
    //

	//
	// Fish School Search
	//

    t = omp_get_wtime();
    fss(input);
    t = omp_get_wtime() - t;
    time = t;

	if(!input->silent)
		printf("FSS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato di xh
	//
	sprintf(fname, "xh64_%d_%d_%d.ds2", input->d, input->np, input->iter);
	save_data(fname, input->xh, 1, input->d);
	if(input->display){
		if(input->xh == NULL)
			printf("xh: NULL\n");
		else{
			printf("xh: [");
			for(i=0; i<input->d-1; i++)
				printf("%f,", input->xh[i]);
			printf("%f]\n", input->xh[i]);
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

    return 0;
} //main
