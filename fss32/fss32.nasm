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

	dim			equ	4		; dimensione in byte di un singolo dato (4 se float, 8 se double)
	p				equ	4		; grado di parallelismo SIMD (4 se float, 2 se double)
	unroll		equ	4		; fattore di unroll
	blocksize	equ	32	; fattore di blocking

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

	v1		equ	8			; puntatore al vettore dei coefficienti
	v2		equ	12		; puntatore al vettore x
	n			equ	16		; dimensione vettori
	ris		equ	20		; puntatore alla variabile contenente il risultato

prodottoScalare:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH		EBP				; salvo il Base Pointer
	MOV		EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH		EBX				; salvo i registri da preservare
	PUSH		ESI
	PUSH		EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV		EAX, [EBP+v1]		; v1
	MOV		EBX, [EBP+v2]		; v2
	MOV		ECX, [EBP+n]		; n
	MOV		EDX, [EBP+ris]		; ris

	;
	; corpo della funzione
	;

	XORPS		XMM0, XMM0		; ps = 0
	XORPS		XMM1, XMM1		; ps = 0
	XORPS		XMM2, XMM2		; ps = 0
	XORPS		XMM3, XMM3		; ps = 0
	XOR 		ESI, ESI				; i = 0

.i:
	MOV 		EDI, ESI				; indTemp = i
	ADD 		EDI, p*unroll		; indTemp+=p
	CMP 		EDI, ECX				; (indTemp > n) ?
	JG			.i_no_unroll			; se vero passa a lavorare senza loop unrolling

	MOVUPS	XMM4, [EAX+ESI*dim+p*0*dim]	; v1[i, ..., i+p-1]
	MOVUPS 	XMM5, [EBX+ESI*dim+p*0*dim]	; v2[i, ..., i+p-1]
	MULPS 		XMM4, XMM5									; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	ADDPS		XMM0, XMM4									; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	MOVUPS	XMM4, [EAX+ESI*dim+p*1*dim]	; v1[i, ..., i+p-1]
	MOVUPS 	XMM5, [EBX+ESI*dim+p*1*dim]	; v2[i, ..., i+p-1]
	MULPS 		XMM4, XMM5									; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	ADDPS		XMM1, XMM4									; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	MOVUPS	XMM4, [EAX+ESI*dim+p*2*dim]	; v1[i, ..., i+p-1]
	MOVUPS 	XMM5, [EBX+ESI*dim+p*2*dim]	; v2[i, ..., i+p-1]
	MULPS 		XMM4, XMM5									; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	ADDPS		XMM2, XMM4									; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	MOVUPS	XMM4, [EAX+ESI*dim+p*3*dim]	; v1[i, ..., i+p-1]
	MOVUPS 	XMM5, [EBX+ESI*dim+p*3*dim]	; v2[i, ..., i+p-1]
	MULPS 		XMM4, XMM5									; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	ADDPS		XMM3, XMM4									; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	ADD			ESI, p*unroll		; i+=p*unroll
	JMP			.i

.i_no_unroll:
	MOV 		EDI, ESI			; indTemp = i
	ADD 		EDI, p				; indTemp+=p
	CMP 		EDI, ECX			; (indTemp > n) ?
	JG			.i_scalar			; se vero passa a lavorare con scalari, anziché vettori

	MOVUPS	XMM4, [EAX+ESI*dim]	; v1[i, ..., i+p-1]
	MOVUPS 	XMM5, [EBX+ESI*dim]	; v2[i, ..., i+p-1]
	MULPS 		XMM4, XMM5					; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	ADDPS		XMM0, XMM4					; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	ADD			ESI, p				; i+=p
	JMP			.i_no_unroll


.i_scalar:
	CMP		ESI, ECX			; (i >= n) ?
	JGE			.end

	MOVSS		XMM4, [EAX+ESI*dim]	; v1[i]
	MULSS		XMM4, [EBX+ESI*dim]	; temp[i] = v1[i] * v2[i]
	ADDSS		XMM0, XMM4					; ps[i, ..., i+p-1] += temp[i]

	INC			ESI				; i++
	JMP			.i_scalar

.end:
	ADDPS		XMM0, XMM1
	ADDPS		XMM0, XMM2
	ADDPS		XMM0, XMM3

	HADDPS	XMM0, XMM0
	HADDPS	XMM0, XMM0

	MOVSS		[EDX], XMM0		; *ris = ps

	;
	;	sequenza di uscita dalla funzione
	;

	POP			EDI				; ripristina i registri da preservare
	POP			ESI
	POP			EBX
	MOV		ESP, EBP		; ripristina lo Stack Pointer
	POP			EBP				; ripristina il Base Pointer
	RET								; ritorna alla funzione chiamante



;
; procedura di calcolo della distanza euclidea
;

global distanzaEuclidea

	v1	equ	8
	v2 	equ	12
	n		equ	16
	ris	equ	20

distanzaEuclidea:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH		EBP
	MOV		EBP, ESP
	PUSH 		EBX
	PUSH 		ESI
	PUSH 		EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV		EAX, [EBP + v1]	; puntatore a v1
	MOV		EBX, [EBP + v2]	; puntatore a v2
	MOV		ECX, [EBP + n]	; n° elementi
	MOV		EDX, [EBP + ris]	; puntatore alla variabile contenente il risultato


	;
	; corpo della funzione
	;


	XORPS		XMM0, XMM0		; conterrà  le somme parziali (v1[i] - v2[i])^2: sono le 4 somme parziali che porto avanti per l'unrolling
	XORPS		XMM1, XMM1
	XORPS		XMM2, XMM2
	XORPS		XMM3, XMM3

	XOR 		ESI, ESI				; i = 0
.i:
	MOV 		EDI, ESI				; temp = i
	ADD 		EDI, p*unroll		; temp += p*unroll		per controllare che siano presenti almeno p*unroll elementi nella prossima iterazione
	CMP		EDI, ECX				; temp < n
	JG			.i_no_unroll

	MOVUPS	XMM4, [EAX + ESI*dim+p*0*dim]	; XMM4 = v1[i ... i+p-1]
	MOVUPS	XMM5, [EAX + ESI*dim+p*1*dim]
	MOVUPS	XMM6, [EAX + ESI*dim+p*2*dim]
	MOVUPS	XMM7, [EAX + ESI*dim+p*3*dim]

	SUBPS		XMM4, [EBX + ESI*dim+p*0*dim] 	; XMM4 -= v2[i ... i+p-1]
	MULPS		XMM4, XMM4 								; XMM4 = XMM4^2
	ADDPS		XMM0, XMM4									; XMM0 += XMM4

	SUBPS		XMM5, [EBX + ESI*dim+p*1*dim]
	MULPS		XMM5, XMM5
	ADDPS		XMM1, XMM5

	SUBPS		XMM6, [EBX + ESI*dim+p*2*dim]
	MULPS		XMM6, XMM6
	ADDPS		XMM2, XMM6

	SUBPS		XMM7, [EBX + ESI*dim+p*3*dim]
	MULPS		XMM7, XMM7
	ADDPS		XMM3, XMM7

	ADD			ESI, p*unroll									; i += p*unroll
	JMP			.i

.i_no_unroll:
	MOV 		EDI, ESI				; temp = i
	ADD 		EDI, p					; temp += p		per controllare che siano presenti almeno p elementi nella prossima iterazione
	CMP		EDI, ECX				; temp < n
	JG			.i_scalar

	MOVUPS	XMM4, [EAX + ESI*dim]	; XMM4 = v1[i ... i+p-1]
	SUBPS		XMM4, [EBX + ESI*dim]	; XMM4 -= v2[i ... i+p-1]
	MULPS		XMM4, XMM4 				; XMM4 = XMM4^2
	ADDPS		XMM0, XMM4					; XMM0 += XMM4

	ADD			ESI, p								; i += p
	JMP			.i_no_unroll

.i_scalar:
	CMP		ESI, ECX							; i < n		gli elementi restanti in n° < p
	JGE			.end

	MOVSS		XMM4, [EAX + ESI*dim]	; XMM4 = v1[i ... i+p-1]
	SUBSS		XMM4, [EBX + ESI*dim]	; XMM4 -= v2[i ... i+p-1]
	MULSS		XMM4, XMM4 				; XMM4 = XMM4^2
	ADDSS		XMM0, XMM4					; XMM0 += XMM4

	INC			ESI									; i++
	JMP			.i_scalar

.end:
	ADDPS		XMM0, XMM1		; somma delle 4 somme parziali dovute all'unrolling
	ADDPS		XMM0, XMM2
	ADDPS		XMM0, XMM3

	HADDPS	XMM0, XMM0
	HADDPS	XMM0, XMM0

	SQRTSS	XMM0, XMM0

	MOVSS		[EDX], XMM0					; *ris = XMM0

	;
	;	sequenza di uscita dalla funzione
	;

	POP			EDI
	POP			ESI
	POP			EBX
	MOV		ESP, EBP
	POP			EBP
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

	PUSH		EBP				; salvo il Base Pointer
	MOV		EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH		EBX				; salvo i registri da preservare
	PUSH		ESI
	PUSH		EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV 		EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
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

	MOV		EBX, [EBP+matrix]	; indirizzo del primo elemento della matrice

	;
	; corpo della funzione


	MOVSS		XMM1, [EBP+sumVect]	; XMM1 = sumV
	SHUFPS 	XMM1, XMM1, 0				; XMM1 = sumV[i, ..., i+p-1] = sumV

	MOV		ESI, 0		; i = 0
.i:
	CMP		ESI, [EAX+16]	; (i<nx) ?
	JGE			.end

	MOV		ECX, [EBP+vect]		; &v

	MOVSS		XMM2, [ECX+ESI*dim]	; XMM2 = v[i]
	SHUFPS 	XMM2, XMM2, 0				; XMM2 = v[i, ..., i+p-1] = v[i]

	MOV		ECX, [EBP+mediaP]			; &mediaP

	MOV		EDI, 0			; j = 0
.j:
	ADD			EDI, p
	CMP		EDI, [EAX+20]		; (j<d) ?
	JG			.end_j
	SUB			EDI, p

	MOV		EDX, ESI					; i
	IMUL		EDX, [EAX+20]		; i*d
	ADD			EDX, EDI					; i*d + j
	IMUL		EDX, dim					; i*d*dim + j*dim

	MOVUPS	XMM0, [EBX + EDX]		; matrix[i][j, ..., j+p-1]
	MULPS		XMM0, XMM2 				; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1]
	DIVPS		XMM0, XMM1					; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]

	ADDPS		XMM0, [ECX+EDI*dim] 	; mediaP[j, ..., j+p-1] += matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]
	MOVAPS	[ECX+EDI*dim], XMM0

	ADD			EDI, p					; j+=4
	JMP			.j

.end_j:
	SUB			EDI, p

.j_scalar:
	CMP		EDI, [EAX+20]
	JGE			.update_i

	MOV		EDX, ESI					; i
	IMUL		EDX, [EAX+20]		; i*d
	ADD			EDX, EDI					; i*d + j
	IMUL		EDX, dim					; i*d*dim + j*dim

	MOVSS		XMM0, [EBX + EDX]		; matrix[i][j]
	MULSS		XMM0, XMM2					; matrix[i][j] * v[i]
	DIVSS		XMM0, XMM1					; matrix[i][j] * v[i] / sumV

	ADDSS		XMM0, [ECX+EDI*dim] 	; mediaP[j] += matrix[i][j] * v[i] / sumV[i]
	MOVSS		[ECX+EDI*dim], XMM0

	INC			EDI							; j++
	JMP			.j_scalar

.update_i:
	INC			ESI							; i++
	JMP			.i

.end:
	; a questo punto in mediaP è stato già messo il risultato e posso quindi terminare


	;
	;	sequenza di uscita dalla funzione
	;

	POP			EDI				; ripristina i registri da preservare
	POP			ESI
	POP			EBX
	MOV		ESP, EBP		; ripristina lo Stack Pointer
	POP			EBP				; ripristina il Base Pointer
	RET							; ritorna alla funzione chiamante




;
; procedura che somma un vettore ad ogni riga di una matrice (aggiornando la matrice)
;

global sommaVettoreMatrice

	input	equ	8			; puntatore al vettore dei parametri
	matrix	equ	12		; puntatore alla matrice
	vect		equ	16		; puntatore al vettore vect

sommaVettoreMatrice:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH		EBP				; salvo il Base Pointer
	MOV		EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH		EBX				; salvo i registri da preservare
	PUSH		ESI
	PUSH		EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV 		EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
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

	MOV		EBX, [EBP+matrix]	; indirizzo del primo elemento della matrice

	MOV		ECX, [EBP+vect]		; indirizzo di vect

	;
	; corpo della funzione
	;

	MOV		ESI, 0						; i = 0
.i:
	CMP		ESI, [EAX+16]			; (i<nx) ?
	JGE			.end

	MOV		EDI, 0						; j = 0
.j:
	ADD			EDI, p*unroll
	CMP		EDI, [EAX+20]			; (j<d) ?
	JG			.end_j
	SUB			EDI, p*unroll

	MOV		EDX, ESI					; i
	IMUL		EDX, [EAX+20]		; i*d
	ADD			EDX, EDI					; i*d + j
	IMUL		EDX, dim					; i*d*dim + j*dim

	MOVUPS	XMM0, [EBX + EDX+p*0*dim]		; matrix[i][j, ..., j+p-1]
	MOVAPS	XMM4, [ECX+EDI*dim+p*0*dim]	; v[j, ..., j+p-1]
	ADDPS		XMM0, XMM4									; matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]

	MOVUPS	XMM1, [EBX + EDX+p*1*dim]
	MOVAPS	XMM4, [ECX+EDI*dim+p*1*dim]
	ADDPS		XMM1, XMM4

	MOVUPS	XMM2, [EBX + EDX+p*2*dim]
	MOVAPS	XMM4, [ECX+EDI*dim+p*2*dim]
	ADDPS		XMM2, XMM4

	MOVUPS	XMM3, [EBX + EDX+p*3*dim]
	MOVAPS	XMM4, [ECX+EDI*dim+p*3*dim]
	ADDPS		XMM3, XMM4

	MOVUPS	[EBX + EDX+p*0*dim], XMM0		; matrix[i][j, ..., j+p-1]  = matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]
	MOVUPS	[EBX + EDX+p*1*dim], XMM1
	MOVUPS	[EBX + EDX+p*2*dim], XMM2
	MOVUPS	[EBX + EDX+p*3*dim], XMM3


	ADD			EDI, p*unroll						; j+=4
	JMP			.j

.end_j:
	SUB			EDI, p*unroll

.j_no_unroll:
	ADD			EDI, p
	CMP		EDI, [EAX+20]			; (j<d) ?
	JG			.end_j_no_unroll
	SUB			EDI, p

	MOV		EDX, ESI					; i
	IMUL		EDX, [EAX+20]		; i*d
	ADD			EDX, EDI					; i*d + j
	IMUL		EDX, dim					; i*d*dim + j*dim

	MOVUPS	XMM0, [EBX + EDX]		; matrix[i][j, ..., j+p-1]
	MOVAPS	XMM1, [ECX+EDI*dim]	; v[j, ..., j+p-1]
	ADDPS		XMM0, XMM1					; matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]
	MOVUPS	[EBX + EDX], XMM0		; matrix[i][j, ..., j+p-1]  = matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]

	ADD			EDI, p						; j+=4
	JMP			.j_no_unroll

.end_j_no_unroll:
	SUB			EDI, p

.j_scalar:
	CMP		EDI, [EAX+20]
	JGE			.update_i

	MOV		EDX, ESI					; i
	IMUL		EDX, [EAX+20]		; i*d
	ADD			EDX, EDI					; i*d + j
	IMUL		EDX, dim					; i*d*dim + j*dim

	MOVSS		XMM0, [EBX + EDX]		; matrix[i][j]
	MOVSS		XMM1, [ECX+EDI*dim]	; v[j]
	ADDSS		XMM0, XMM1				; matrix[i][j] + v[j]
	MOVSS		[EBX + EDX], XMM0		; matrix[i][j]  = matrix[i][j] + v[j]

	INC			EDI							; j++
	JMP			.j_scalar

.update_i:
	INC			ESI							; i++
	JMP			.i

.end:
	; a questo punto la matrice è stata già modificata e posso terminare

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
; procedura di generazione della possibile nuova posizione dei pesci nel movimento individuale
;

global generaPossibileMovimento

	align 16
	uno				dd	1.0, 1.0, 1.0, 1.0		; 1

	input			equ	8		; puntatore al vettore dei parametri
	randIndex	equ	12	; puntatore all'indice che scorre il file dei numeri random
	i					equ	16	; indice del pesce i-esimo
	y					equ	20	; puntatore al vettore delle possibili nuove posizioni

generaPossibileMovimento:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH		EBP					; salvo il Base Pointer
	MOV		EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH		EBX					; salvo i registri da preservare
	PUSH		ESI
	PUSH		EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV 		EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
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

	MOV		EBX, [EAX]				; indirizzo del primo elemento della matrice

	MOV		EDX, [EBP+randIndex]		; indirizzo dell'indice dei numeri random
	MOV		ECX, [EDX]				; indice dei numeri random

	;
	; corpo della funzione
	;


	MOVSS		XMM4, [EAX+28]		; input->stepind
	SHUFPS	XMM4, XMM4, 0			; input->stepind[j, ..., j+p-1]


	XOR			EDI, EDI					; j = 0
.j:
	ADD			EDI, p*unroll
	CMP		EDI, [EAX + 20]		; (j<d) ?
	JG			.end_j
	SUB			EDI, p*unroll

	MOV		EDX, [EAX+12]									; &input->r
	MOVUPS	XMM0, [EDX+ECX*dim+p*0*dim]		; input->r[ri, ..., ri+p-1]
	ADDPS		XMM0, XMM0										; input->r[ri, ..., ri+p-1] * 2
	SUBPS		XMM0, [uno]										; input->r[ri, ..., ri+p-1] * 2 - 1

	MOVUPS	XMM1, [EDX+ECX*dim+p*1*dim]
	ADDPS		XMM1, XMM1
	SUBPS		XMM1, [uno]

	MOVUPS	XMM2, [EDX+ECX*dim+p*2*dim]
	ADDPS		XMM2, XMM2
	SUBPS		XMM2, [uno]

	MOVUPS	XMM3, [EDX+ECX*dim+p*3*dim]
	ADDPS		XMM3, XMM3
	SUBPS		XMM3, [uno]

	ADD			ECX, p*unroll

	MULPS		XMM0, XMM4				; (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]
	MULPS		XMM1, XMM4
	MULPS		XMM2, XMM4
	MULPS		XMM3, XMM4

	MOV		ESI, [EBP+i]					; i
	IMUL		ESI, [EAX+20]				; i*d
	ADD			ESI, EDI						; i*d + j
	IMUL		ESI, ESI, dim				; i*d*dim + j*dim

	MOVUPS	XMM5, [EBX+ESI+p*0*dim]	; input->x[i*d*dim + j*dim]
	ADDPS		XMM0, XMM5							; input->x[i*d*dim + j*dim]+(input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]

	MOVUPS	XMM5, [EBX+ESI+p*1*dim]
	ADDPS		XMM1, XMM5

	MOVUPS	XMM5, [EBX+ESI+p*2*dim]
	ADDPS		XMM2, XMM5

	MOVUPS	XMM5, [EBX+ESI+p*3*dim]
	ADDPS		XMM3, XMM5

	MOV		EDX, [EBP+y]
	MOVAPS	[EDX+EDI*dim+p*0*dim], XMM0
	MOVAPS	[EDX+EDI*dim+p*1*dim], XMM1
	MOVAPS	[EDX+EDI*dim+p*2*dim], XMM2
	MOVAPS	[EDX+EDI*dim+p*3*dim], XMM3

	ADD			EDI, p*unroll
	JMP			.j

.end_j:
	SUB			EDI, p*unroll

.j_no_unroll:
	ADD			EDI, p
	CMP		EDI, [EAX + 20]		; (j<d) ?
	JG			.end_j_no_unroll
	SUB			EDI, p

	MOV		EDX, [EAX+12]					; &input->r
	MOVUPS	XMM0, [EDX+ECX*dim]		; input->r[ri, ..., ri+p-1]
	ADDPS		XMM0, XMM0						; input->r[ri, ..., ri+p-1] * 2
	SUBPS		XMM0, [uno]						; input->r[ri, ..., ri+p-1] * 2 - 1

	ADD			ECX, p

	MULPS		XMM0, XMM4				; (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]

	MOV		ESI, [EBP+i]					; i
	IMUL		ESI, [EAX+20]				; i*d
	ADD			ESI, EDI						; i*d + j
	IMUL		ESI, ESI, dim				; i*d*dim + j*dim

	MOVUPS	XMM5, [EBX+ESI]		; input->x[i*d*dim + j*dim]
	ADDPS		XMM0, XMM5				; input->x[i*d*dim + j*dim] + (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]

	MOV		EDX, [EBP+y]
	MOVAPS	[EDX+EDI*dim], XMM0

	ADD			EDI, p
	JMP			.j_no_unroll

.end_j_no_unroll:
	SUB			EDI, p

.j_scalar:
	CMP		EDI, [EAX+20]
	JGE			.end

	MOV		EDX, [EAX+12]					; &input->r
	MOVSS		XMM0, [EDX+ECX*dim]		; input->r[j]
	ADDSS		XMM0, XMM0						; input->r[ri] * 2
	SUBSS		XMM0, [uno]						; input->r[ri] * 2 - 1

	INC			ECX

	MULSS		XMM0, XMM4				; (input->r[ri] * 2 - 1) * input->stepind

	MOV		ESI, [EBP+i]					; i
	IMUL		ESI, [EAX+20]				; i*d
	ADD			ESI, EDI						; i*d + j
	IMUL		ESI, ESI, dim				; i*d*dim + j*dim

	MOVSS		XMM5, [EBX+ESI]		; input->x[i*d*dim + j*dim]
	ADDSS		XMM0, XMM5				; input->x[i*d*dim + j*dim] + (input->r[jri * 2 - 1) * input->stepind[j]

	MOV		EDX, [EBP+y]
	MOVSS		[EDX+EDI*dim], XMM0

	INC			EDI
	JMP			.j_scalar

.end:
	MOV		EDX, [EBP+randIndex]		; aggiorno randIndex
	MOV		[EDX], ECX
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

	x					equ	8		; puntatore alla matrice delle posizioni dei pesci
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
	MOVUPS	XMM1, [EAX+ESI]			; input->x[i*d*dim + j*dim]
	SUBPS		XMM0, XMM1					; y[j, ..., j+p-1] - input->x[i*d*dim + j*dim]

	MOVUPS	[ECX+ESI], XMM0			; dx[i*d*dim+j*dim] = y[j, ..., j+p-1] - input->x[i*d*dim + j*dim]
	MOVUPS	[EAX+ESI], XMM2			; input->x[i*d*dim + j*dim] = y[j, ..., j+p-1]

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

	deltaX2		equ	8		; puntatore alla matrice delle variazioni di posizione
	i3					equ	12	; indice del pesce i-esimo
	d2				equ	16	; input->d

mantieniPosizionePesce:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH		EBP					; salvo il Base Pointer
	MOV		EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH		EBX					; salvo i registri da preservare
	PUSH		ESI
	PUSH		EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV		EAX, [EBP+deltaX2]		; matrice dx

	;
	; corpo della funzione
	;

	XOR			EDI, EDI					; j = 0
.j:
	ADD			EDI, p*unroll
	CMP		EDI, [EBP+d2]		; (j<d) ?
	JG			.end_j
	SUB			EDI, p*unroll

	MOV		ESI, [EBP+i3]					; i
	IMUL		ESI, [EBP+d2]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim

	XORPS		XMM0, XMM0						; 0
	MOVUPS	[EAX+ESI+p*0*dim], XMM0	; dx[i*d*dim + j*dim] = 0

	XORPS		XMM1, XMM1
	MOVUPS	[EAX+ESI+p*1*dim], XMM1

	XORPS		XMM2, XMM2
	MOVUPS	[EAX+ESI+p*2*dim], XMM2

	XORPS		XMM3, XMM3
	MOVUPS	[EAX+ESI+p*3*dim], XMM3

	ADD			EDI, p*unroll
	JMP			.j

.end_j:
	SUB			EDI, p*unroll

.j_no_unroll:
	ADD			EDI, p
	CMP		EDI, [EBP+d2]		; (j<d) ?
	JG			.end_j_no_unroll
	SUB			EDI, p

	MOV		ESI, [EBP+i3]					; i
	IMUL		ESI, [EBP+d2]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim

	XORPS		XMM0, XMM0					; 0
	MOVUPS	[EAX+ESI], XMM0			; dx[i*d*dim + j*dim] = 0

	ADD			EDI, p
	JMP			.j_no_unroll

.end_j_no_unroll:
	SUB			EDI, p

.j_scalar:
	CMP		EDI, [EBP+d2]
	JGE			.end

	MOV		ESI, [EBP+i3]					; i
	IMUL		ESI, [EBP+d2]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim

	XORPS		XMM0, XMM0					; 0
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

	MOV	EDX, [EBP+weightGain]	; 1 se il peso del banco è aumentato, 0 altrimenti

	;
	; corpo della funzione
	;

	MOVSS 	XMM5, [EAX+32]			; stepvol
	SHUFPS	XMM5, XMM5, 0				; stepvol[j, ..., j+p-1]

	MOVSS		XMM6, [EBP+distEuclidea]			; dist
	SHUFPS	XMM6, XMM6, 0				; dist[j, ..., j+p-1]

	MOVSS		XMM7, [EBP+randNum]	; randNum
	SHUFPS	XMM7, XMM7, 0				; randNum[j, ..., j+p-1]

	XOR			EDI, EDI				; j = 0
.j:
	ADD			EDI, p*unroll
	CMP		EDI, [EAX+20]		; (j<d) ?
	JG			.end_j
	SUB			EDI, p*unroll

	MOV		ESI, [EBP+i4]		; i
	IMUL		ESI, [EAX+20]		; i*d
	ADD			ESI, EDI				; i*d + j
	IMUL		ESI, ESI, dim		; i*d*dim + j*dim

	MOVUPS	XMM0, [EBX+ESI+p*0*dim]			; x[i*d+j]
	SUBPS		XMM0, [ECX+EDI*dim+p*0*dim]	; x[i*d+j] - B[j]
	MULPS		XMM0, XMM7					; (x[i*d+j]-B[j])*randNum
	MULPS		XMM0, XMM5					; (x[i*d+j]-B[j])*randNum*stepvol
	DIVPS		XMM0, XMM6					; (x[i*d+j]-B[j])*randNum*stepvol / dist (= numerator)

	MOVUPS	XMM1, [EBX+ESI+p*1*dim]
	SUBPS		XMM1, [ECX+EDI*dim+p*1*dim]
	MULPS		XMM1, XMM7
	MULPS		XMM1, XMM5
	DIVPS		XMM1, XMM6

	MOVUPS	XMM2, [EBX+ESI+p*2*dim]
	SUBPS		XMM2, [ECX+EDI*dim+p*2*dim]
	MULPS		XMM2, XMM7
	MULPS		XMM2, XMM5
	DIVPS		XMM2, XMM6

	MOVUPS	XMM3, [EBX+ESI+p*3*dim]
	SUBPS		XMM3, [ECX+EDI*dim+p*3*dim]
	MULPS		XMM3, XMM7
	MULPS		XMM3, XMM5
	DIVPS		XMM3, XMM6

	CMP		EDX, 0								; if (weightGain)
	JE				.false

	MOVUPS	XMM4, [EBX+ESI+p*0*dim]		; x[i*d+j]
	SUBPS		XMM4, XMM0								; x[i*d+j] - numerator
	MOVUPS	[EBX+ESI+p*0*dim], XMM4		; x[i*d+j] = x[i*d+j] - numerator

	MOVUPS	XMM4, [EBX+ESI+p*1*dim]
	SUBPS		XMM4, XMM1
	MOVUPS	[EBX+ESI+p*1*dim], XMM4

	MOVUPS	XMM4, [EBX+ESI+p*2*dim]
	SUBPS		XMM4, XMM2
	MOVUPS	[EBX+ESI+p*2*dim], XMM4

	MOVUPS	XMM4, [EBX+ESI+p*3*dim]
	SUBPS		XMM4, XMM3
	MOVUPS	[EBX+ESI+p*3*dim], XMM4

	JMP			.end_false

.false:
	MOVUPS	XMM4, [EBX+ESI+p*0*dim]		; x[i*d+j]
	ADDPS		XMM4, XMM0								; x[i*d+j] - numerator
	MOVUPS	[EBX+ESI+p*0*dim], XMM4		; x[i*d+j] = x[i*d+j] - numerator

	MOVUPS	XMM4, [EBX+ESI+p*1*dim]
	ADDPS		XMM4, XMM1
	MOVUPS	[EBX+ESI+p*1*dim], XMM4

	MOVUPS	XMM4, [EBX+ESI+p*2*dim]
	ADDPS		XMM4, XMM2
	MOVUPS	[EBX+ESI+p*2*dim], XMM4

	MOVUPS	XMM4, [EBX+ESI+p*3*dim]
	ADDPS		XMM4, XMM3
	MOVUPS	[EBX+ESI+p*3*dim], XMM4

.end_false:

	ADD			EDI, p*unroll
	JMP			.j

.end_j:
	SUB			EDI, p*unroll

.j_no_unroll:
	ADD			EDI, p
	CMP		EDI, [EAX+20]		; (j<d) ?
	JG			.end_j_no_unroll
	SUB			EDI, p

	MOV		ESI, [EBP+i4]		; i
	IMUL		ESI, [EAX+20]		; i*d
	ADD			ESI, EDI				; i*d + j
	IMUL		ESI, ESI, dim		; i*d*dim + j*dim

	MOVUPS	XMM0, [EBX+ESI]			; x[i*d+j]
	SUBPS		XMM0, [ECX+EDI*dim]	; x[i*d+j] - B[j]
	MULPS		XMM0, XMM7					; (x[i*d+j]-B[j])*randNum
	MULPS		XMM0, XMM5					; (x[i*d+j]-B[j])*randNum*stepvol
	DIVPS		XMM0, XMM6					; (x[i*d+j]-B[j])*randNum*stepvol / dist (= numerator)


	CMP		EDX, 0								; if (weightGain)
	JE				.false_no_unroll

	MOVUPS	XMM4, [EBX+ESI]			; x[i*d+j]
	SUBPS		XMM4, XMM0					; x[i*d+j] - numerator
	MOVUPS	[EBX+ESI], XMM4			; x[i*d+j] = x[i*d+j] - numerator

	JMP			.end_false_no_unroll

.false_no_unroll:
	MOVUPS	XMM4, [EBX+ESI]			; x[i*d+j]
	ADDPS		XMM4, XMM0					; x[i*d+j] + numerator
	MOVUPS	[EBX+ESI], XMM4			; x[i*d+j] = x[i*d+j] - numerator

.end_false_no_unroll:

	ADD			EDI, p
	JMP			.j_no_unroll

.end_j_no_unroll:
	SUB			EDI, p

.j_scalar:
	CMP		EDI, [EAX+20]
	JGE			.end

	MOV		ESI, [EBP+i4]					; i
	IMUL		ESI, [EAX+20]					; i*d
	ADD			ESI, EDI							; i*d + j
	IMUL		ESI, ESI, dim					; i*d*dim + j*dim

	MOVSS		XMM0, [EBX+ESI]			; x[i*d+j]
	SUBSS		XMM0, [ECX+EDI*dim]	; x[i*d+j] - B[j]
	MULSS		XMM0, XMM7					; (x[i*d+j]-B[j])*randNum
	MULSS 	XMM0, XMM5					; (x[i*d+j]-B[j])*randNum*stepvol
	DIVSS		XMM0, XMM6					; (x[i*d+j]-B[j])*randNum*stepvol / dist (= numerator)


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



;
; procedura che preso un vettore ne somma gli elementi
;

global sommaElementiVettore

	input		equ		8
	v				equ		12
	sumV		equ		16


sommaElementiVettore:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH		EBP					; salvo il Base Pointer
	MOV		EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH		EBX					; salvo i registri da preservare
	PUSH		ESI
	PUSH		EDI

	;
	; lettura dei parametri dal record di attivazione
	;

	MOV 		EAX, [EBP+input]	; indirizzo della struttura contenente i parametri
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

	MOV		EBX, [EBP+v]			; indirizzo al primo valore di df
	MOV		ECX, [EBP+sumV]	; indirizzo a sumdf

	;
	; corpo della funzione
	;

	XORPS		XMM0, XMM0			; sumV
	XOR			ESI, ESI					; i = 0
.i:
	ADD			ESI, p
	CMP		ESI, [EAX+16]			; (i<np) ?
	JG			.end_i
	SUB			ESI, p

	MOVAPS	XMM1, [EBX+ESI*dim]	; v[i, ..., i+p-1]
	ADDPS		XMM0, XMM1					; sumV+=v[i, ..., i+p-1]

	ADD			ESI, p
	JMP			.i

.end_i:
	SUB			ESI, p

.i_scalar:
	CMP		ESI, [EAX+16]
	JGE			.end

	MOVSS		XMM1, [EBX+ESI*dim]	; v[i]
	ADDSS		XMM0, XMM1					; sumV+=v[i]

	INC			ESI
	JMP			.i_scalar

.end:
	HADDPS	XMM0, XMM0		; effettuo le due somme orizzontali rimanenti
	HADDPS	XMM0, XMM0

	MOVSS		[ECX], XMM0		; *sumV = sumV


	;
	;	sequenza di uscita dalla funzione
	;

	POP			EDI				; ripristina i registri da preservare
	POP			ESI
	POP			EBX
	MOV		ESP, EBP		; ripristina lo Stack Pointer
	POP			EBP				; ripristina il Base Pointer
	RET								; ritorna alla funzione chiamante
