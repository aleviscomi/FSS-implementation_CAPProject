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

	indirizzoMediaP		resd		1	; locazione di memoria utile per la procedura mediaPesata
	sommaV			resd		1	; locazione che conterrà la sommatoria degli elementi del vettore vect (usata in mediaPesata)
	input_nx			resd		1	; locazione che conterrà le righe della matrice
	input_d			resd		1	; locazione che conterrà le colonne della matrice

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

	vect1	equ	8		; puntatore al vettore dei coefficienti
	vect2	equ	12		; puntatore al vettore x
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

	MOV	EAX, [EBP+vect1]		; x
	MOV	EBX, [EBP+vect2]		; y
	MOV	ECX, [EBP+n]		; n
	MOV	EDX, [EBP+ris]		; ris

	;
	; corpo della funzione
	;

	XORPS	XMM0, XMM0		; ps = 0
	XOR 	ESI, ESI			; i = 0

for_i:
	MOV 	EDI, ESI			; indTemp = i
	ADD 	EDI, p				; indTemp+=p
	CMP 	EDI, ECX			; (indTemp > n) ?
	JG		for_i_scalar		; se vero passa a lavorare con scalari, anziché vettori

	MOVAPS	XMM1, [EAX+ESI*dim]	; vect1[i, ..., i+p-1]
	MULPS 	XMM1, [EBX+ESI*dim]	; temp[i, ..., i+p-1] = vect1[i, ..., i+p-1] * vect2[i, ..., i+p-1]
	ADDPS	XMM0, XMM1			; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	ADD		ESI, p			; i+=p
	JMP		for_i

for_i_scalar:
	CMP	ESI, ECX			; (i >= n) ?
	JGE		end

	MOVSS	XMM1, [EAX+ESI*dim]	; vect1[i]
	MULSS	XMM1, [EBX+ESI*dim]	; temp[i] = vect1[i] * vect2[i]
	ADDSS	XMM0, XMM1			; ps[i, ..., i+p-1] += temp[i]

	INC		ESI				; i++
	JMP		for_i_scalar

end:
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

loopQ:
	MOV 	EDI, ESI				; temp = i
	ADD 	EDI, p				; temp += p		per controllare che siano presenti almeno p elementi nella prossima iterazione
	CMP	EDI, ECX				; temp < n
	JGE		loopR
	MOVAPS	XMM0, [EAX + ESI*dim]	; XMM0 = v1[i ... i+p-1]
	SUBPS	XMM0, [EBX + ESI*dim]	; XMM0 -= v2[i ... i+p-1]
	MULPS	XMM0, XMM0			; XMM0 = XMM0^2
	ADDPS	XMM1, XMM0			; XMM1 += XMM0
	ADD		ESI, p				; i += p
	JMP		loopQ

loopR:
	CMP	ESI, ECX				; i < n		gli elementi restanti in nÂ° < p
	JGE		endLoop
	MOVSS	XMM0, [EAX + ESI*dim]	; XMM0 = v1[i]
	SUBSS	XMM0, [EBX + ESI*dim]	; XMM0 -= v2[i]
	MULSS	XMM0, XMM0			; XMM0 = XMM0^2
	ADDSS	XMM1, XMM0			; XMM1 += XMM0
	INC		ESI					; i++
	JMP		loopR

endLoop:
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

	input	equ	8		; puntatore al vettore dei parametri
	matrix	equ	12		; puntatore alla matrice
	vect		equ	16		; puntatore al vettore vect
	mediaP	equ	20		; puntatore al vettore contenente il risultato
	sumVect	equ	24		; contiene la somma degli elementi di vect

mediaPesata:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH	EBP				; salvo il Base Pointer
	MOV	EBP, ESP			; il Base Pointer punta al record di attivazione corrente
	PUSH	ESP				; salvo i registri da preservare
	PUSH	EBX
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

	MOV	EBX, [EAX+16]	; memorizzo il numero di righe della matrice (nx)
	MOV	[input_nx], EBX

	MOV	EBX, [EAX+20]	; memorizzo il numero di colonne della matrice (d)
	MOV	[input_d], EBX

	MOV	EAX, [EBP+matrix]	; indirizzo del primo elemento della matrice

	MOV	EBX, [EBP+vect]	; indirizzo di vect

	MOV	ECX, [EBP+mediaP]		; memorizzo l'indirizzo di mediaP in memoria per usarlo dopo
	MOV	[indirizzoMediaP], ECX	; è necessario ciò perché dopo uso EBP e quindi non potrò più
								; accederci con [EBP+mediaP]

	MOV	ECX, [EBP+sumVect]	; memorizzo sumVect in quanto perdendo dopo EBP
	MOV	[sommaV], ECX		; non ci potrei più accedere

	;
	; corpo della funzione

	MOV	EBP, 0		; i = 0
for_i2:
	IMUL	ESI, EBP, dim		; i*dim
	IMUL	ESI, [input_d]		; i*d*dim

	CMP	EBP, [input_nx]	; (i<nx) ?
	JGE		end2

	MOV	ECX, 0			; j = 0
for_j2:
	ADD		ECX, p
	CMP	ECX, [input_d]		; (j<d) ?
	JG		end_for_j2
	SUB		ECX, p

	IMUL	EDI, ECX, dim		; j*dim

	MOV	EDX, ESI			; i*d*dim
	ADD		EDX, EDI			; i*d*dim + j*dim

	MOVAPS	XMM0, [EAX + EDX]	; matrix[i][j, ..., j+p-1]
	MOVSS	XMM1, [EBX+EBP*dim]	; v[i]
	SHUFPS 	XMM1, XMM1, 0		; v[i, ..., i+p-1] = v[i]
	MULPS	XMM0, XMM1			; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1]
	MOVSS	XMM1, [sommaV]		; sumV
	SHUFPS 	XMM1, XMM1, 0		; sumV[i, ..., i+p-1] = sumV
	DIVPS	XMM0, XMM1			; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]

	MOV	EDX, [indirizzoMediaP]	; indirizzo di mediaP
	ADDPS	XMM0, [EDX+EDI] 		; mediaP[j, ..., j+p-1] = matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]
	MOVAPS	[EDX+EDI], XMM0

	ADD		ECX, p
	JMP		for_j2

end_for_j2:
	SUB		ECX, p

for_j_scalar2:
	CMP	ECX, [input_d]
	JGE		update_for_i2

	IMUL	EDI, ECX, dim		; j*dim

	MOV	EDX, ESI			; i*d*dim
	ADD		EDX, EDI			; i*d*dim + j*dim

	MOVSS	XMM0, [EAX + EDX]	; matrix[i][j]
	MOVSS	XMM1, [EBX+EBP*dim]	; v[i]
	MULSS	XMM0, XMM1			; matrix[i][j] * v[i]
	DIVSS	XMM0, [sommaV]		; matrix[i][j] * v[i] / sumV

	MOV	EDX, [indirizzoMediaP]	; indirizzo di mediaP
	ADDSS	XMM0, [EDX+EDI] 		; mediaP[j] = matrix[i][j] * v[i] / sumV
	MOVSS	[EDX+EDI], XMM0

	INC		ECX
	JMP		for_j_scalar2

update_for_i2:
	INC		EBP
	JMP		for_i2

end2:
	; a questo punto in mediaP è stato già messo il risultato e posso quindi terminare


	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	POP		ESP				; ripristina lo Stack Pointer
	POP		EBP				; ripristina il Base Pointer
	RET




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
	PUSH	ESP				; salvo i registri da preservare
	PUSH	EBX
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

	MOV	EBX, [EAX+16]	; memorizzo il numero di righe della matrice (nx)
	MOV	[input_nx], EBX

	MOV	EBX, [EAX+20]	; memorizzo il numero di colonne della matrice (d)
	MOV	[input_d], EBX


	MOV	EAX, [EBP+matrix]	; indirizzo del primo elemento della matrice

	MOV	EBX, [EBP+vect]	; indirizzo di vect

	;
	; corpo della funzione
	;

	MOV	EBP, 0			; i = 0
for_i3:
	IMUL	ESI, EBP, dim		; i*dim
	IMUL	ESI, [input_d]		; i*d*dim

	CMP	EBP, [input_nx]	; (i<nx) ?
	JGE		end3

	MOV	ECX, 0			; j = 0
for_j3:
	ADD		ECX, p
	CMP	ECX, [input_d]		; (j<d) ?
	JG		end_for_j3
	SUB		ECX, p

	IMUL	EDI, ECX, dim		; j*dim

	MOV	EDX, ESI			; i*d*dim
	ADD		EDX, EDI			; i*d*dim + j*dim

	MOVAPS	XMM0, [EAX + EDX]	; matrix[i][j, ..., j+p-1]
	MOVAPS	XMM1, [EBX+ECX*dim]	; v[j, ..., j+p-1]
	ADDPS		XMM0, XMM1			; matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]
	MOVAPS	[EAX + EDX], XMM0	; matrix[i][j, ..., j+p-1]  = matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]

	ADD		ECX, p
	JMP		for_j3

end_for_j3:
	SUB		ECX, p

for_j_scalar3:
	CMP	ECX, [input_d]
	JGE		update_for_i3

	IMUL	EDI, ECX, dim		; j*dim

	MOV	EDX, ESI			; i*d*dim
	ADD		EDX, EDI			; i*d*dim + j*dim

	MOVSS	XMM0, [EAX + EDX]	; matrix[i][j]
	MOVSS	XMM1, [EBX+ECX*dim]	; v[j]
	ADDSS	XMM0, XMM1			; matrix[i][j] + v[j]
	MOVSS	[EAX + EDX], XMM0	; matrix[i][j]  = matrix[i][j] + v[j]

	INC		ECX
	JMP		for_j_scalar3

update_for_i3:
	INC		EBP
	JMP		for_i3

end3:
	; a questo punto la matrice è stata già modificata e posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	POP		EDI				; ripristina i registri da preservare
	POP		ESI
	POP		EBX
	POP		ESP				; ripristina lo Stack Pointer
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
	MOV	ECX, [EDX]

	;
	; corpo della funzione
	;

	XOR		EDI, EDI					; j = 0
for_j4:
	ADD		EDI, p
	CMP	EDI, [EAX + 20]		; (j<d) ?
	JG		end_for_j4
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
	JMP		for_j4

end_for_j4:
	SUB		EDI, p

for_j_scalar4:
	CMP	EDI, [EAX+20]
	JGE		end4

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
	JMP			for_j_scalar4

end4:
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
