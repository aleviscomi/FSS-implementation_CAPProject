; ---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf32 fss32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati

	dim	equ	4		; dimensione in byte di un singolo dato (4 se float, 8 se double)
	p	equ	4		; grado di parallelismo SIMD (4 se float, 2 se double)

section .bss			; Sezione contenente dati non inizializzati

section .text			; Sezione contenente il codice macchina


; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------


;
; procedura di calcolo del prodotto scalare
;

global prodottoScalare

	v1	equ	8		; puntatore al vettore dei coefficienti
	v2	equ	12		; puntatore al vettore x
	n			equ	16		; dimensione vettori
	ris		equ	20		; puntatore alla variabile contenente il risultato

prodottoScalare:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP				; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	EBX				; salvo i registri da preservare
	PUSH	ESI
	PUSH	EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV	EAX, [EBP+v1]		; v1
	MOV	EBX, [EBP+v2]		; v2
	MOV	ECX, [EBP+n]		; n
	MOV	EDX, [EBP+ris]		; ris

	;
	; corpo della funzione
	;

	XORPS	XMM0, XMM0		; ps = 0
	XOR 	ESI, ESI			; i = 0

.i:
	MOV 	EDI, ESI			; indTemp = i
	ADD 	EDI, p				; indTemp+=p
	CMP 	EDI, ECX			; (indTemp > n) ?
	JG		.i_scalar		; se vero passa a lavorare con scalari, anziché vettori

	MOVAPS	XMM1, [EAX+ESI*dim]	; v1[i, ..., i+p-1]
	MULPS 	XMM1, [EBX+ESI*dim]	; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	ADDPS	XMM0, XMM1			; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	ADD		ESI, p			; i+=p
	JMP		.i

.i_scalar:
	CMP	ESI, ECX			; (i >= n) ?
	JGE		.end

	MOVSS	XMM1, [EAX+ESI*dim]	; v1[i]
	MULSS	XMM1, [EBX+ESI*dim]	; temp[i] = v1[i] * v2[i]
	ADDSS	XMM0, XMM1			; ps[i, ..., i+p-1] += temp[i]

	INC		ESI				; i++
	JMP		.i_scalar

.end:
	HADDPS	XMM0, XMM0		; effettuo le due somme orizzontali rimanenti
	HADDPS	XMM0, XMM0

	MOVSS	[EDX], XMM0		; *ris = ps

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	MOV	ESP, EBP			; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET						; ritorna alla funzione chiamante



;
; procedura di calcolo della distanza euclidea
;

global distanzaEuclidea

	v1	equ	8
	v2 	equ	12
	n	equ	16
	ris	equ	20

distanzaEuclidea:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP
	MOV	EBP, ESP
	PUSH 	EBX
	PUSH 	ESI
	PUSH 	EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV	EAX, [EBP + v1]	; puntatore a v1
	MOV	EBX, [EBP + v2]	; puntatore a v2
	MOV	ECX, [EBP + n]	; n° elementi
	MOV	EDX, [EBP + ris]	; puntatore alla variabile contenente il risultato

	; INIZIALIZZAZIONI

	XOR 	ESI, ESI				; i = 0
	XORPS	XMM1, XMM1			; conterrà  le somme parziali (v1[i] - v2[i])^2


	;
	; corpo della funzione
	;

.i:
	MOV 	EDI, ESI				; temp = i
	ADD 	EDI, p				; temp += p		per controllare che siano presenti almeno p elementi nella prossima iterazione
	CMP	EDI, ECX				; temp < n
	JGE		.i_scalar
	MOVAPS	XMM0, [EAX + ESI*dim]	; XMM0 = v1[i ... i+p-1]
	SUBPS	XMM0, [EBX + ESI*dim]	; XMM0 -= v2[i ... i+p-1]
	MULPS	XMM0, XMM0			; XMM0 = XMM0^2
	ADDPS	XMM1, XMM0			; XMM1 += XMM0
	ADD		ESI, p				; i += p
	JMP		.i

.i_scalar:
	CMP	ESI, ECX				; i < n		gli elementi restanti in nÂ° < p
	JGE		.end
	MOVSS	XMM0, [EAX + ESI*dim]	; XMM0 = v1[i]
	SUBSS	XMM0, [EBX + ESI*dim]	; XMM0 -= v2[i]
	MULSS	XMM0, XMM0			; XMM0 = XMM0^2
	ADDSS	XMM1, XMM0			; XMM1 += XMM0
	INC		ESI					; i++
	JMP		.i_scalar

.end:
	HADDPS	XMM1, XMM1
	HADDPS	XMM1, XMM1
	SQRTSS	XMM1, XMM1

	MOVSS	[EDX], XMM1			; *ris = XMM1

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI
	POP		ESI
	POP		EBX
	MOV	ESP, EBP
	POP		EBP
	RET



;
; procedura di calcolo della media pesata
;

global mediaPesata

	input		equ	8			; puntatore al vettore dei parametri
	matrix		equ	12		; puntatore alla matrice
	vect			equ	16		; puntatore al vettore vect
	mediaP		equ	20		; puntatore al vettore contenente il risultato
	sumVect	equ	24		; contiene la somma degli elementi di vect

mediaPesata:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP				; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	EBX				; salvo i registri da preservare
	PUSH	ESI
	PUSH	EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV 	EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
			; [EAX]	input->x
			; [EAX + 4] input->xh
			; [EAX + 8] input->c
			; [EAX + 12] input->r
			; [EAX + 16] input->nx
			; [EAX + 20] input->d
			; [EAX + 24] input->iter
			; [EAX + 28] input->stepind
			; [EAX + 32] input->stepvol
			; [EAX + 36] input->wscale

	MOV	EBX, [EBP+matrix]	; indirizzo del primo elemento della matrice

	;
	; corpo della funzione

	MOV	ESI, 0		; i = 0
.i:
	CMP	ESI, [EAX+16]	; (i<nx) ?
	JGE		.end

	MOV	EDI, 0			; j = 0
.j:
	ADD		EDI, p
	CMP	EDI, [EAX+20]		; (j<d) ?
	JG		.end_j
	SUB		EDI, p

	MOV	EDX, ESI					; i
	IMUL	EDX, [EAX+20]		; i*d
	ADD		EDX, EDI					; i*d + j
	IMUL	EDX, dim					; i*d*dim + j*dim

	MOV	ECX, [EBP+vect]		; &v

	MOVAPS	XMM0, [EBX + EDX]		; matrix[i][j, ..., j+p-1]
	MOVSS		XMM1, [ECX+ESI*dim]	; v[i]
	SHUFPS 	XMM1, XMM1, 0				; v[i, ..., i+p-1] = v[i]
	MULPS		XMM0, XMM1					; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1]

	MOVSS		XMM1, [EBP+sumVect]	; sumV
	SHUFPS 	XMM1, XMM1, 0				; sumV[i, ..., i+p-1] = sumV
	DIVPS		XMM0, XMM1					; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]

	MOV		ECX, [EBP+mediaP]			; &mediaP
	ADDPS		XMM0, [ECX+EDI*dim] 	; mediaP[j, ..., j+p-1] += matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]
	MOVAPS	[ECX+EDI*dim], XMM0

	ADD		EDI, p					; j+=4
	JMP		.j

.end_j:
	SUB		EDI, p

.j_scalar:
	CMP	EDI, [EAX+20]
	JGE		.update_i

	MOV	EDX, ESI					; i
	IMUL	EDX, [EAX+20]		; i*d
	ADD		EDX, EDI					; i*d + j
	IMUL	EDX, dim					; i*d*dim + j*dim

	MOV	ECX, [EBP+vect]		; &v

	MOVSS	XMM0, [EBX + EDX]		; matrix[i][j]
	MOVSS	XMM1, [ECX+ESI*dim]	; v[i]
	MULSS	XMM0, XMM1					; matrix[i][j] * v[i]

	DIVSS	XMM0, [EBP+sumVect]	; matrix[i][j] * v[i] / sumV

	MOV	ECX, [EBP+mediaP]			; &mediaP
	ADDSS	XMM0, [ECX+EDI*dim] 	; mediaP[j] += matrix[i][j] * v[i] / sumV[i]
	MOVSS	[ECX+EDI*dim], XMM0

	INC		EDI							; j++
	JMP		.j_scalar

.update_i:
	INC		ESI							; i++
	JMP		.i

.end:
	; a questo punto in mediaP è stato già messo il risultato e posso quindi terminare


	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	MOV	ESP, EBP		; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET							; ritorna alla funzione chiamante




;
; procedura che somma un vettore ad ogni riga di una matrice (aggiornando la matrice)
;

global sommaVettoreMatrice

	input	equ	8		; puntatore al vettore dei parametri
	matrix	equ	12		; puntatore alla matrice
	vect		equ	16		; puntatore al vettore vect

sommaVettoreMatrice:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP				; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	EBX				; salvo i registri da preservare
	PUSH	ESI
	PUSH	EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV 	EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
			; [EAX]	input->x
			; [EAX + 4] input->xh
			; [EAX + 8] input->c
			; [EAX + 12] input->r
			; [EAX + 16] input->nx
			; [EAX + 20] input->d
			; [EAX + 24] input->iter
			; [EAX + 28] input->stepind
			; [EAX + 32] input->stepvol
			; [EAX + 36] input->wscale

	MOV	EBX, [EBP+matrix]	; indirizzo del primo elemento della matrice

	MOV	ECX, [EBP+vect]		; indirizzo di vect

	;
	; corpo della funzione
	;

	MOV	ESI, 0		; i = 0
.i:
	CMP	ESI, [EAX+16]	; (i<nx) ?
	JGE		.end

	MOV	EDI, 0			; j = 0
.j:
	ADD		EDI, p
	CMP	EDI, [EAX+20]		; (j<d) ?
	JG		.end_j
	SUB		EDI, p

	MOV	EDX, ESI					; i
	IMUL	EDX, [EAX+20]		; i*d
	ADD		EDX, EDI					; i*d + j
	IMUL	EDX, dim					; i*d*dim + j*dim

	MOVAPS	XMM0, [EBX + EDX]		; matrix[i][j, ..., j+p-1]
	MOVAPS	XMM1, [ECX+EDI*dim]	; v[j, ..., j+p-1]
	ADDPS		XMM0, XMM1					; matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]
	MOVAPS	[EBX + EDX], XMM0		; matrix[i][j, ..., j+p-1]  = matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]

	ADD		EDI, p						; j+=4
	JMP		.j

.end_j:
	SUB		EDI, p

.j_scalar:
	CMP	EDI, [EAX+20]
	JGE		.update_i

	MOV	EDX, ESI					; i
	IMUL	EDX, [EAX+20]		; i*d
	ADD		EDX, EDI					; i*d + j
	IMUL	EDX, dim					; i*d*dim + j*dim

	MOVSS	XMM0, [EBX + EDX]	; matrix[i][j]
	MOVSS	XMM1, [ECX+EDI*dim]	; v[j]
	ADDSS	XMM0, XMM1			; matrix[i][j] + v[j]
	MOVSS	[EBX + EDX], XMM0	; matrix[i][j]  = matrix[i][j] + v[j]

	INC		EDI							; j++
	JMP		.j_scalar

.update_i:
	INC		ESI							; i++
	JMP		.i

.end:
	; a questo punto la matrice è stata già modificata e posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	MOV	ESP, EBP			; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET						; ritorna alla funzione chiamante



;
; procedura di generazione della possibile nuova posizione dei pesci nel movimento individuale
;

global generaPossibileMovimento

	align 16
	uno				dd	1.0, 1.0, 1.0, 1.0		; 1
	align 16
	due				dd	2.0, 2.0, 2.0, 2.0		; 2

	input			equ	8		; puntatore al vettore dei parametri
	randIndex	equ	12	; puntatore all'indice che scorre il file dei numeri random
	i					equ	16	; indice del pesce i-esimo
	y					equ	20	; puntatore al vettore delle possibili nuove posizioni

generaPossibileMovimento:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP					; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	EBX					; salvo i registri da preservare
	PUSH	ESI
	PUSH	EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV 	EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
			; [EAX]	input->x
			; [EAX + 4] input->xh
			; [EAX + 8] input->c
			; [EAX + 12] input->r
			; [EAX + 16] input->nx
			; [EAX + 20] input->d
			; [EAX + 24] input->iter
			; [EAX + 28] input->stepind
			; [EAX + 32] input->stepvol
			; [EAX + 36] input->wscale

	MOV	EBX, [EAX]				; indirizzo del primo elemento della matrice

	MOV	EDX, [EBP+randIndex]		; indirizzo dell'indice dei numeri random
	MOV	ECX, [EDX]				; indice dei numeri random

	;
	; corpo della funzione
	;

	XOR		EDI, EDI					; j = 0
.j:
	ADD		EDI, p
	CMP	EDI, [EAX + 20]		; (j<d) ?
	JG		.end_j
	SUB		EDI, p

	MOV		EDX, [EAX+12]			; &input->r
	MOVAPS	XMM0, [EDX+ECX*dim]		; input->r[ri, ..., ri+p-1]
	MULPS		XMM0, [due]				; input->r[ri, ..., ri+p-1] * 2
	SUBPS		XMM0, [uno]				; input->r[ri, ..., ri+p-1] * 2 - 1

	ADD			ECX, p

	MOVSS		XMM1, [EAX+28]		; input->stepind
	SHUFPS	XMM1, XMM1, 0			; input->stepind[j, ..., j+p-1]
	MULPS		XMM0, XMM1				; (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]

	MOV		ESI, [EBP+i]					; i
	IMUL		ESI, [EAX+20]				; i*d
	ADD			ESI, EDI						; i*d + j
	IMUL		ESI, ESI, dim				; i*d*dim + j*dim

	MOVAPS	XMM1, [EBX+ESI]		; input->x[i*d*dim + j*dim]
	ADDPS		XMM0, XMM1				; input->x[i*d*dim + j*dim] + (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]

	MOV		EDX, [EBP+y]
	MOVAPS	[EDX+EDI*dim], XMM0

	ADD		EDI, p
	JMP		.j

.end_j:
	SUB		EDI, p

.j_scalar:
	CMP	EDI, [EAX+20]
	JGE		.end

	MOV		EDX, [EAX+12]			; &input->r
	MOVSS		XMM0, [EDX+ECX*dim]		; input->r[j]
	MULSS		XMM0, [due]				; input->r[ri] * 2
	SUBSS		XMM0, [uno]				; input->r[ri] * 2 - 1

	INC			ECX

	MOVSS		XMM1, [EAX+28]		; input->stepind
	MULSS		XMM0, XMM1				; (input->r[ri] * 2 - 1) * input->stepind

	MOV		ESI, [EBP+i]					; i
	IMUL		ESI, [EAX+20]				; i*d
	ADD			ESI, EDI						; i*d + j
	IMUL		ESI, ESI, dim				; i*d*dim + j*dim

	MOVSS		XMM1, [EBX+ESI]		; input->x[i*d*dim + j*dim]
	ADDSS		XMM0, XMM1				; input->x[i*d*dim + j*dim] + (input->r[jri * 2 - 1) * input->stepind[j]

	MOV		EDX, [EBP+y]
	MOVSS		[EDX+EDI*dim], XMM0

	INC			EDI
	JMP			.j_scalar

.end:
	MOV	EDX, [EBP+randIndex]		; aggiorno randIndex
	MOV	[EDX], ECX
	; a questo punto y è stato già modificato e posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	MOV	ESP, EBP			; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET						; ritorna alla funzione chiamante




;
; procedura che muove il pesce i-esimo nella nuova posizione y (aggiornando dx)
;

global muoviPesce

	x					equ	8		; puntatore al vettore dei parametri
	y2				equ	12	; puntatore al vettore delle possibili nuove posizioni
	deltaX			equ	16	; puntatore alla matrice delle variazioni di posizione
	i2					equ	20	; indice del pesce i-esimo
	d					equ	24	; input->d

muoviPesce:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP					; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	EBX					; salvo i registri da preservare
	PUSH	ESI
	PUSH	EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV	EAX, [EBP+x]				; matrice x
	MOV	EBX, [EBP+y2]				; vettore y
	MOV	ECX, [EBP+deltaX]		; matrice dx

	;
	; corpo della funzione
	;

	XOR		EDI, EDI					; j = 0
.j:
	ADD		EDI, p
	CMP	EDI, [EBP+d]		; (j<d) ?
	JG		.end_j
	SUB		EDI, p

	MOVAPS	XMM0, [EBX+EDI*dim]	; y[j, ..., j+p-1]
	MOVAPS	XMM2, XMM0					; copio y[j, ..., j+p-1] per salvarlo in x e non rientrare in memoria
	MOV		ESI, [EBP+i2]					; i
	IMUL		ESI, [EBP+d]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim
	MOVAPS	XMM1, [EAX+ESI]			; input->x[i*d*dim + j*dim]
	SUBPS		XMM0, XMM1					; y[j, ..., j+p-1] - input->x[i*d*dim + j*dim]

	MOVAPS	[ECX+ESI], XMM0			; dx[i*d*dim+j*dim] = y[j, ..., j+p-1] - input->x[i*d*dim + j*dim]
	MOVAPS	[EAX+ESI], XMM2			; input->x[i*d*dim + j*dim] = y[j, ..., j+p-1]

	ADD		EDI, p
	JMP		.j

.end_j:
	SUB		EDI, p

.j_scalar:
	CMP	EDI, [EBP+d]
	JGE		.end

	MOVSS		XMM0, [EBX+EDI*dim]	; y[j]
	MOVSS		XMM2, XMM0					; copio y[j] per salvarlo in x e non rientrare in memoria
	MOV		ESI, [EBP+i2]					; i
	IMUL		ESI, [EBP+d]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim
	MOVSS		XMM1, [EAX+ESI]			; input->x[i*d*dim + j*dim]
	SUBSS		XMM0, XMM1					; y[j] - input->x[i*d*dim + j*dim]

	MOVSS		[ECX+ESI], XMM0			; dx[i*d*dim+j*dim] = y[j] - input->x[i*d*dim + j*dim]
	MOVSS		[EAX+ESI], XMM2			; input->x[i*d*dim + j*dim] = y[j]

	INC			EDI
	JMP			.j_scalar

.end:
	; a questo punto dx e x sono stati già modificati e posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	MOV	ESP, EBP			; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET						; ritorna alla funzione chiamante


;
; procedura che lascia il pesce i-esimo nella sua posizione (dx[i][j]=0)
;

global mantieniPosizionePesce

	align 16
	zero				dd	0.0, 0.0, 0.0, 0.0		; 0

	deltaX2		equ	8		; puntatore alla matrice delle variazioni di posizione
	i3					equ	12	; indice del pesce i-esimo
	d2				equ	16	; input->d


mantieniPosizionePesce:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP					; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	EBX					; salvo i registri da preservare
	PUSH	ESI
	PUSH	EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV	EAX, [EBP+deltaX2]		; matrice dx

	;
	; corpo della funzione
	;

	XOR		EDI, EDI					; j = 0
.j:
	ADD		EDI, p
	CMP	EDI, [EBP+d2]		; (j<d) ?
	JG		.end_j
	SUB		EDI, p

	MOV		ESI, [EBP+i3]					; i
	IMUL		ESI, [EBP+d2]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim
	MOVAPS	XMM0, [zero]					; 0
	MOVAPS	[EAX+ESI], XMM0			; dx[i*d*dim + j*dim] = 0

	ADD			EDI, p
	JMP			.j

.end_j:
	SUB			EDI, p

.j_scalar:
	CMP		EDI, [EBP+d2]
	JGE			.end

	MOV		ESI, [EBP+i3]					; i
	IMUL		ESI, [EBP+d2]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim
	MOVSS		XMM0, [zero]					; 0
	MOVSS		[EAX+ESI], XMM0			; dx[i*d*dim + j*dim] = 0

	INC			EDI
	JMP			.j_scalar

.end:
	; a questo punto dx e x sono stati già modificati e posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	MOV	ESP, EBP			; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET						; ritorna alla funzione chiamante




;
; procedura che esegue il movimento volitivo per tutti i pesci
;

global faiMovimentoVolitivo

	input				equ	8
	b						equ	12
	distEuclidea		equ	16
	randNum			equ	20
	weightGain		equ	24
	i4						equ	28


faiMovimentoVolitivo:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP					; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	EBX					; salvo i registri da preservare
	PUSH	ESI
	PUSH	EDI

	;
	; lettura dei parametri dal record di attivazione
	;


	MOV 	EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
			; [EAX]	input->x
			; [EAX + 4] input->xh
			; [EAX + 8] input->c
			; [EAX + 12] input->r
			; [EAX + 16] input->nx
			; [EAX + 20] input->d
			; [EAX + 24] input->iter
			; [EAX + 28] input->stepind
			; [EAX + 32] input->stepvol
			; [EAX + 36] input->wscale

	MOV	EBX, [EAX]				; indirizzo della matrice x

	MOV	ECX, [EBP+b]			; indirizzo del vettore B

	;
	; corpo della funzione
	;

	XOR			EDI, EDI				; j = 0
.j:
	ADD			EDI, p
	CMP		EDI, [EAX+20]		; (j<d) ?
	JG			.end_j
	SUB			EDI, p

	MOV		ESI, [EBP+i4]		; i
	IMUL		ESI, [EAX+20]		; i*d
	ADD			ESI, EDI				; i*d + j
	IMUL		ESI, ESI, dim		; i*d*dim + j*dim

	MOVAPS	XMM0, [EBX+ESI]			; x[i*d+j]
	MOVAPS	XMM1, [ECX+EDI*dim]	; B[j]
	SUBPS		XMM0, XMM1					; x[i*d+j] - B[j]

	MOVSS		XMM1, [EBP+randNum]	; randNum
	SHUFPS	XMM1, XMM1, 0				; randNum[j, ..., j+p-1]
	MULPS		XMM0, XMM1					; (x[i*d+j]-B[j])*randNum

	MOVSS 	XMM1, [EAX+32]			; stepvol
	SHUFPS	XMM1, XMM1, 0				; stepvol[j, ..., j+p-1]
	MULPS		XMM0, XMM1					; (x[i*d+j]-B[j])*randNum*stepvol

	MOVSS		XMM1, [EBP+distEuclidea]			; dist
	SHUFPS	XMM1, XMM1, 0				; dist[j, ..., j+p-1]
	DIVPS		XMM0, XMM1					; (x[i*d+j]-B[j])*randNum*stepvol / dist (= numerator)


	MOV		EDX, [EBP+weightGain]
	CMP		EDX, 0								; if (weightGain)
	JE				.false
	MOVAPS	XMM1, [EBX+ESI]			; x[i*d+j]
	SUBPS		XMM1, XMM0					; x[i*d+j] - numerator
	MOVAPS	[EBX+ESI], XMM1			; x[i*d+j] = x[i*d+j] - numerator
	JMP			.end_false
.false:
	MOVAPS	XMM1, [EBX+ESI]			; x[i*d+j]
	ADDPS		XMM1, XMM0					; x[i*d+j] + numerator
	MOVAPS	[EBX+ESI], XMM1			; x[i*d+j] = x[i*d+j] - numerator
.end_false:

	ADD			EDI, p
	JMP			.j

.end_j:
	SUB			EDI, p

.j_scalar:
	CMP		EDI, [EAX+20]
	JGE			.end

	MOV		ESI, [EBP+i4]					; i
	IMUL		ESI, [EAX+20]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim

	MOVSS		XMM0, [EBX+ESI]			; x[i*d+j]
	MOVSS		XMM1, [ECX+EDI*dim]	; B[j]
	SUBSS		XMM0, XMM1					; x[i*d+j] - B[j]

	MULSS		XMM0, [EBP+randNum]	; (x[i*d+j]-B[j])*randNum

	MULSS 	XMM0, [EAX+32]			; (x[i*d+j]-B[j])*randNum*stepvol

	DIVSS		XMM0, [EBP+distEuclidea]			; (x[i*d+j]-B[j])*randNum*stepvol / dist (= numerator)

	MOV		EDX, [EBP+weightGain]
	CMP		EDX, 0								; if (weightGain)
	JE				.false_scalar
	MOVSS		XMM1, [EBX+ESI]			; x[i*d+j]
	SUBSS		XMM1, XMM0					; x[i*d+j] - numerator
	MOVSS		[EBX+ESI], XMM1			; x[i*d+j] = x[i*d+j] - numerator
	JMP			.end_false_scalar
.false_scalar:
	MOVSS		XMM1, [EBX+ESI]			; x[i*d+j]
	ADDSS		XMM1, XMM0					; x[i*d+j] + numerator
	MOVSS		[EBX+ESI], XMM1			; x[i*d+j] = x[i*d+j] - numerator
.end_false_scalar:

	INC			EDI
	JMP			.j_scalar

.end:
	; a questo punto posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	MOV	ESP, EBP			; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET						; ritorna alla funzione chiamante

