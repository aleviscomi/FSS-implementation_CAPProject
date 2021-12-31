; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
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
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati

	dim			equ	8		; dimensione in byte di un singolo dato (4 se float, 8 se double)
	p				equ	4		; grado di parallelismo SIMD (8 se float, 4 se double)
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
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------


;
; procedura di calcolo del prodotto scalare
;

global prodottoScalare

prodottoScalare:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH		RBP				; salva il Base Pointer
	MOV		RBP, RSP	; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	;
	; lettura dei parametri dal record di attivazione
	;


	; RDI		; v1 	puntatore al vettore dei coefficienti
	; RSI		; v2	puntatore al vettore x
	; RDX	; n 	dimensione vettori
	; RCX	; ris	puntatore alla variabile contenente il risultato


	;
	; corpo della funzione
	;

	VXORPD	YMM0, YMM0		; ps = 0
	VXORPD	YMM1, YMM1		; ps = 0
	VXORPD	YMM2, YMM2		; ps = 0
	VXORPD	YMM3, YMM3		; ps = 0
	XOR 		RBX, RBX			; i = 0

.i:
	MOV 		RAX, RBX			; indTemp = i
	ADD 		RAX, p*unroll		; indTemp+=p
	CMP 		RAX, RDX			; (indTemp > n) ?
	JG			.i_no_unroll			; se vero passa a lavorare senza loop unrolling

	VMOVUPD	YMM4, [RDI+RBX*dim+p*0*dim]	; v1[i, ..., i+p-1]
	VMOVUPD	YMM5, [RSI+RBX*dim+p*0*dim]	; v2[i, ..., i+p-1]
	VMULPD 	YMM4, YMM5									; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	VADDPD	YMM0, YMM4									; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	VMOVUPD	YMM4, [RDI+RBX*dim+p*1*dim]	; v1[i+p, ..., i+2p-1]
	VMOVUPD	YMM5, [RSI+RBX*dim+p*1*dim]	; v2[i, ..., i+p-1]
	VMULPD 	YMM4, YMM5									; temp[i+p, ..., i+2p-1] = v1[i+p, ..., i+2p-1] * v2[i+p, ..., i+2p-1]
	VADDPD	YMM1, YMM4 								; ps[i+p, ..., i+2p-1]+= temp[i+p, ..., i+2p-1]

	VMOVUPD	YMM4, [RDI+RBX*dim+p*2*dim]	; v1[i+2p, ..., i+3p-1]
	VMOVUPD	YMM5, [RSI+RBX*dim+p*2*dim]	; v2[i, ..., i+p-1]
	VMULPD 	YMM4, YMM5									; temp[i+2p, ..., i+3p-1]= v1[i+2p, ..., i+3p-1] * v2[i+2p, ..., i+3p-1]
	VADDPD	YMM2, YMM4									; ps[i+2p, ..., i+3p-1] += temp[i+2p, ..., i+3p-1]

	VMOVUPD	YMM4, [RDI+RBX*dim+p*3*dim]	; v1[i+3p, ..., i+4p-1]
	VMOVUPD	YMM5, [RSI+RBX*dim+p*3*dim]	; v2[i, ..., i+p-1]
	VMULPD 	YMM4, YMM5									; temp[i+3p, ..., i+4p-1] = v1[i+3p, ..., i+4p-1] * v2[i+3p, ..., i+4p-1]
	VADDPD	YMM3, YMM4									; ps[i+3p, ..., i+4p-1] += temp[i+3p, ..., i+4p-1]

	ADD			RBX, p*unroll		; i+=p*unroll
	JMP			.i

.i_no_unroll:
	MOV 		RAX, RBX			; indTemp = i
	ADD 		RAX, p					; indTemp+=p
	CMP 		RAX, RDX			; (indTemp > n) ?
	JG			.end_i_no_unroll	; se vero passa a lavorare con scalari, anziché vettori

	VMOVUPD	YMM4, [RDI+RBX*dim]	; v1[i, ..., i+p-1]
	VMOVUPD YMM5, [RSI+RBX*dim]	; v2[i, ..., i+p-1]
	VMULPD 	YMM4, YMM5					; temp[i, ..., i+p-1] = v1[i, ..., i+p-1] * v2[i, ..., i+p-1]
	VADDPD	YMM0, YMM4					; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	ADD			RBX, p				; i+=p
	JMP			.i_no_unroll

.end_i_no_unroll:
	VXORPD	YMM4, YMM4		; YMM4 = 0, mantiene la somma del ciclo scalare (poiché AVX su XMM azzererebbe gli altri bit)

.i_scalar:
	CMP		RBX, RDX			; (i >= n) ?
	JGE			.end

	VMOVSD	XMM5, [RDI+RBX*dim]	; v1[i]
	VMULSD	XMM5, [RSI+RBX*dim]	; temp[i] = v1[i] * v2[i]
	VADDSD	XMM4, XMM5					; ps[i, ..., i+p-1] += temp[i]

	INC			RBX				; i++
	JMP			.i_scalar

.end:
	VADDPD	YMM0, YMM1
	VADDPD	YMM0, YMM2
	VADDPD	YMM0, YMM3

	VADDPD	YMM0, YMM4

	VHADDPD	YMM0, YMM0

	VPERM2F128	YMM1, YMM0, YMM0, 10000001b

	VADDSD	XMM0, XMM1

	VMOVSD	[RCX], XMM0		; *ris = ps

	;
	;	sequenza di uscita dalla funzione
	;

	popaq							; ripristina i registri generali
	MOV		RSP, RBP	; ripristina lo Stack Pointer
	POP			RBP				; ripristina il Base Pointer
	RET								; torna alla funzione C chiamante


;
; procedura di calcolo della distanza euclidea
;

global distanzaEuclidea

distanzaEuclidea:
	;
	; sequenza di ingresso nella funzione
	;

	PUSH		RBP				; salva il Base Pointer
	MOV		RBP, RSP	; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	;
	; lettura dei parametri dal record di attivazione
	;


	; RDI		; v1 	puntatore a v1
	; RSI		; v2	puntatore a v2
	; RDX	; n 	n° elementi
	; RCX	; ris	puntatore alla variabile contenente il risultato


	;
	; corpo della funzione
	;

	VXORPD	YMM0, YMM0		; conterrà  le somme parziali (v1[i] - v2[i])^2: sono le 4 somme parziali che porto avanti per l'unrolling
	VXORPD	YMM1, YMM1
	VXORPD	YMM2, YMM2
	VXORPD	YMM3, YMM3

	XOR 		RBX, RBX				; i = 0
.i:
	MOV 		RAX, RBX				; temp = i
	ADD 		RAX, p*unroll		; temp += p*unroll		per controllare che siano presenti almeno p*unroll elementi nella prossima iterazione
	CMP		RAX, RDX				; temp < n
	JG			.i_no_unroll

	VMOVUPD	YMM4, [RDI + RBX*dim+p*0*dim]	; XMM4 = v1[i ... i+p-1]
	VMOVUPD	YMM5, [RDI + RBX*dim+p*1*dim]
	VMOVUPD	YMM6, [RDI + RBX*dim+p*2*dim]
	VMOVUPD	YMM7, [RDI + RBX*dim+p*3*dim]

	VSUBPD	YMM4, [RSI + RBX*dim+p*0*dim] 	; XMM4 -= v2[i ... i+p-1]
	VMULPD	YMM4, YMM4 								; XMM4 = XMM4^2
	VADDPD	YMM0, YMM4									; XMM0 += XMM4

	VSUBPD	YMM5, [RSI + RBX*dim+p*1*dim]
	VMULPD	YMM5, YMM5
	VADDPD	YMM1, YMM5

	VSUBPD	YMM6, [RSI + RBX*dim+p*2*dim]
	VMULPD	YMM6, YMM6
	VADDPD	YMM2, YMM6

	VSUBPD	YMM7, [RSI + RBX*dim+p*3*dim]
	VMULPD	YMM7, YMM7
	VADDPD	YMM3, YMM7

	ADD			RBX, p*unroll									; i += p*unroll
	JMP			.i

.i_no_unroll:
	MOV 		RAX, RBX					; temp = i
	ADD 		RAX, p						; temp += p		per controllare che siano presenti almeno p elementi nella prossima iterazione
	CMP		RAX, RDX				; temp < n
	JG			.end_i_no_unroll

	VMOVUPD	YMM4, [RDI + RBX*dim]	; XMM4 = v1[i ... i+p-1]
	VSUBPD	YMM4, [RSI + RBX*dim]	; XMM4 -= v2[i ... i+p-1]
	VMULPD	YMM4, YMM4 				; XMM4 = XMM4^2
	VADDPD	YMM0, YMM4					; XMM0 += XMM4

	ADD			RBX, p								; i += p
	JMP			.i_no_unroll

.end_i_no_unroll:
	VXORPD	YMM4, YMM4		; YMM4 = 0, mantiene la somma del ciclo scalare (poiché AVX su XMM azzererebbe gli altri bit)

.i_scalar:
	CMP		RBX, RDX							; i < n		gli elementi restanti in n° < p
	JGE			.end

	VMOVSD	XMM5, [RDI + RBX*dim]	; XMM5 = v1[i ... i+p-1]
	VSUBSD	XMM5, [RSI + RBX*dim]	; XMM5 -= v2[i ... i+p-1]
	VMULSD	XMM5, XMM5 				; XMM5 = XMM5^2
	VADDSD	XMM4, XMM5					; XMM4 += XMM5

	INC			RBX									; i++
	JMP			.i_scalar

.end:
	VADDPD	YMM0, YMM1		; somma delle 4 somme parziali dovute all'unrolling
	VADDPD	YMM0, YMM2
	VADDPD	YMM0, YMM3

	VADDPD	YMM0, YMM4

	VHADDPD	YMM0, YMM0

	VPERM2F128	YMM1, YMM0, YMM0, 10000001b
	VADDSD	XMM0, XMM1

	VSQRTSD	XMM0, XMM0

	VMOVSD	[RCX], XMM0		; *ris = XMM0


	;
	;	sequenza di uscita dalla funzione
	;

	popaq							; ripristina i registri generali
	MOV		RSP, RBP	; ripristina lo Stack Pointer
	POP			RBP				; ripristina il Base Pointer
	RET								; torna alla funzione C chiamante



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

	PUSH		RBP				; salva il Base Pointer
	MOV		RBP, RSP	; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	;
	; lettura dei parametri dal record di attivazione
	;

	; RDI		; input		indirizzo della struttura contenente i parametri
		; [RDI]	input->x
		; [RDI + 8] input->xh
		; [RDI + 16] input->c
		; [RDI + 24] input->r
		; [RDI + 32] input->nx
		; [RDI + 36] input->d
		; [RDI + 40] input->iter
		; [RDI + 48] input->stepind
		; [RDI + 56] input->stepvol
		; [RDI + 64] input->wscale
		; ...
	; RSI		; matrix	indirizzo del primo elemento della matrice
	; RDX	; vect		indirizzo del primo elemento del vettore
	; RCX	; mediaP	indirizzo al vettore contenente il risultato
	; XMM0	; sumVect	contiene la somma degli elementi di vect

	MOV		R11d, [RDI+32]		; poiché nx è un numero a 32bit (non a 64) uso R11d per spostare solo quei 32 bit nella parte bassa di R11
	MOV		R12d, [RDI+36]		; stesso discorso per d
													; poiché ho bisogno di più registri uso questi a 64bit (non avrò problemi
													; poiché la mov sui 32bit bassi azzera gli alti)

	;
	; corpo della funzione
	;

	VSHUFPD			XMM0, XMM0, 0						; XMM4 = sumV[i, ..., i+p-1]
	VPERM2F128	YMM1, YMM0, YMM0, 0		; YMM4 = sumV[i, ..., i+p-1] = sumV

	MOV		RBX, 0		; i = 0
.i:
	CMP		RBX, R11	; (i<nx) ?
	JGE			.end

	VMOVSD			XMM2, [RDX+RBX*dim]			; XMM6 = v[i]
	VSHUFPD 			XMM2, XMM2, 0						; XMM6 = v[i, ..., i+p-1]
	VPERM2F128	YMM2, YMM2, YMM2, 0		; YMM6 = v[i, ..., i+p-1] = v[i]

	MOV		RAX, 0			; j = 0
.j:
	ADD			RAX, p
	CMP		RAX, R12		; (j<d) ?
	JG			.end_j
	SUB			RAX, p

	MOV		R10, RBX					; i
	IMUL		R10, R12			; i*d
	ADD			R10, RAX					; i*d + j
	IMUL		R10, dim					; i*d*dim + j*dim

	VMOVUPD	YMM0, [RSI + R10]			; matrix[i][j, ..., j+p-1]
	VMULPD	YMM0, YMM2 				; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1]
	VDIVPD		YMM0, YMM1					; matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]

	VADDPD	YMM0, [RCX+RAX*dim] 	; mediaP[j, ..., j+p-1] += matrix[i][j, ..., j+p-1] * v[i, ..., i+p-1] / sumV[i, ..., i+p-1]
	VMOVAPD	[RCX+RAX*dim], YMM0

	ADD			RAX, p					; j+=4
	JMP			.j

.end_j:
	SUB			RAX, p

.j_scalar:
	CMP		RAX, R12
	JGE			.update_i

	MOV		R10, RBX					; i
	IMUL		R10, R12					; i*d
	ADD			R10, RAX					; i*d + j
	IMUL		R10, dim					; i*d*dim + j*dim

	VMOVSD	XMM0, [RSI + R10]			; matrix[i][j]
	VMULSD	XMM0, XMM2					; matrix[i][j] * v[i]
	VDIVSD	XMM0, XMM1					; matrix[i][j] * v[i] / sumV

	VADDSD	XMM0, [RCX+RAX*dim] 	; mediaP[j] += matrix[i][j] * v[i] / sumV[i]
	VMOVSD	[RCX+RAX*dim], XMM0

	INC			RAX							; j++
	JMP			.j_scalar

.update_i:
	INC			RBX							; i++
	JMP			.i

.end:
	; a questo punto in mediaP è stato già messo il risultato e posso quindi terminare


	;
	;	sequenza di uscita dalla funzione
	;

	popaq							; ripristina i registri generali
	MOV		RSP, RBP	; ripristina lo Stack Pointer
	POP			RBP				; ripristina il Base Pointer
	RET								; torna alla funzione C chiamante





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

	PUSH		RBP				; salva il Base Pointer
	MOV		RBP, RSP	; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	;
	; lettura dei parametri dal record di attivazione
	;


	; RDI		; input		indirizzo della struttura contenente i parametri
		; [RDI]	input->x
		; [RDI + 8] input->xh
		; [RDI + 16] input->c
		; [RDI + 24] input->r
		; [RDI + 32] input->nx
		; [RDI + 36] input->d
		; [RDI + 40] input->iter
		; [RDI + 48] input->stepind
		; [RDI + 56] input->stepvol
		; [RDI + 64] input->wscale
		; ...
	; RSI		; matrix	indirizzo del primo elemento della matrice
	; RDX	; vect		indirizzo del primo elemento del vettore

	MOV		R11d, [RDI+32]		; poiché nx è un numero a 32bit (non a 64) uso R11d per spostare solo quei 32 bit nella parte bassa di R11
	MOV		R12d, [RDI+36]		; stesso discorso per d
													; poiché ho bisogno di più registri uso questi a 64bit (non avrò problemi
													; poiché la mov sui 32bit bassi azzera gli alti)

	;
	; corpo della funzione
	;

	MOV		RBX, 0						; i = 0
.i:
	CMP		RBX, R11					; (i<nx) ?
	JGE			.end

	MOV		RAX, 0						; j = 0
.j:
	ADD			RAX, p*unroll
	CMP		RAX, R12					; (j<d) ?
	JG			.end_j
	SUB			RAX, p*unroll

	MOV		R10, RBX					; i
	IMUL		R10, R12					; i*d
	ADD			R10, RAX					; i*d + j
	IMUL		R10, dim					; i*d*dim + j*dim

	VMOVUPD	YMM0, [RSI + R10+p*0*dim]		; matrix[i][j, ..., j+p-1]
	VMOVAPD	YMM4, [RDX+RAX*dim+p*0*dim]	; v[j, ..., j+p-1]
	VADDPD	YMM0, YMM4									; matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]

	VMOVUPD	YMM1, [RSI + R10+p*1*dim]
	VMOVAPD	YMM4, [RDX+RAX*dim+p*1*dim]
	VADDPD	YMM1, YMM4

	VMOVUPD	YMM2, [RSI + R10+p*2*dim]
	VMOVAPD	YMM4, [RDX+RAX*dim+p*2*dim]
	VADDPD	YMM2, YMM4

	VMOVUPD	YMM3, [RSI + R10+p*3*dim]
	VMOVAPD	YMM4, [RDX+RAX*dim+p*3*dim]
	VADDPD	YMM3, YMM4

	VMOVUPD	[RSI + R10+p*0*dim], YMM0		; matrix[i][j, ..., j+p-1]  = matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]
	VMOVUPD	[RSI + R10+p*1*dim], YMM1
	VMOVUPD	[RSI + R10+p*2*dim], YMM2
	VMOVUPD	[RSI + R10+p*3*dim], YMM3


	ADD			RAX, p*unroll						; j+=4
	JMP			.j

.end_j:
	SUB			RAX, p*unroll

.j_no_unroll:
	ADD			RAX, p
	CMP		RAX, R12					; (j<d) ?
	JG			.end_j_no_unroll
	SUB			RAX, p

	MOV		R10, RBX					; i
	IMUL		R10, R12					; i*d
	ADD			R10, RAX					; i*d + j
	IMUL		R10, dim					; i*d*dim + j*dim

	VMOVUPD	YMM0, [RSI + R10]			; matrix[i][j, ..., j+p-1]
	VMOVAPD	YMM1, [RDX+RAX*dim]	; v[j, ..., j+p-1]
	VADDPD	YMM0, YMM1					; matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]
	VMOVUPD	[RSI + R10], YMM0			; matrix[i][j, ..., j+p-1]  = matrix[i][j, ..., j+p-1] + v[j, ..., j+p-1]

	ADD			RAX, p						; j+=4
	JMP			.j_no_unroll

.end_j_no_unroll:
	SUB			RAX, p

.j_scalar:
	CMP		RAX, R12
	JGE			.update_i

	MOV		R10, RBX					; i
	IMUL		R10, R12					; i*d
	ADD			R10, RAX					; i*d + j
	IMUL		R10, dim					; i*d*dim + j*dim

	VMOVSD	XMM0, [RSI + R10]			; matrix[i][j]
	VMOVSD	XMM1, [RDX+RAX*dim]	; v[j]
	VADDSD	XMM0, XMM1					; matrix[i][j] + v[j]
	VMOVSD	[RSI + R10], XMM0			; matrix[i][j]  = matrix[i][j] + v[j]

	INC			RAX							; j++
	JMP			.j_scalar

.update_i:
	INC			RBX							; i++
	JMP			.i

.end:
	; a questo punto la matrice è stata già modificata e posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	popaq							; ripristina i registri generali
	MOV		RSP, RBP	; ripristina lo Stack Pointer
	POP			RBP				; ripristina il Base Pointer
	RET								; torna alla funzione C chiamante




;
; procedura di generazione della possibile nuova posizione dei pesci nel movimento individuale
;

global generaPossibileMovimento

	align 32
	uno				dq	1.0, 1.0, 1.0, 1.0		; 1

generaPossibileMovimento:

	;
	; sequenza di ingresso nella funzione
	;

	PUSH		RBP				; salva il Base Pointer
	MOV		RBP, RSP	; il Base Pointer punta al Record di Attivazione corrente
	pushaq						; salva i registri generali

	;
	; lettura dei parametri dal record di attivazione
	;


	; RDI		; input		indirizzo della struttura contenente i parametri
		; [RDI]	input->x
		; [RDI + 8] input->xh
		; [RDI + 16] input->c
		; [RDI + 24] input->r
		; [RDI + 32] input->nx
		; [RDI + 36] input->d
		; [RDI + 40] input->iter
		; [RDI + 48] input->stepind
		; [RDI + 56] input->stepvol
		; [RDI + 64] input->wscale
		; ...
	; RSI		; randIndex	puntatore all'indice che scorre il file dei numeri random
	; RDX	; i					indice del pesce i-esimo
	; RCX	; y				puntatore al vettore delle possibili nuove posizioni

	MOV		R10, [RDI] 				; indirizzo del primo elemento della matrice
	MOV		R11, [RDI+24]			; indirizzo del primo elemento dei numeri random

	MOV		R12d, [RDI+36]

	MOV		R13d, [RSI]
	;
	; corpo della funzione
	;


	XOR			RAX, RAX					; j = 0
.j:
	ADD			RAX, p*unroll
	CMP		RAX, R12						; (j<d) ?
	JG			.end_j
	SUB			RAX, p*unroll

	VMOVUPD	YMM0, [R11+R13*dim+p*0*dim]		; input->r[ri, ..., ri+p-1]
	VADDPD	YMM0, YMM0										; input->r[ri, ..., ri+p-1] * 2
	VSUBPD	YMM0, [uno]										; input->r[ri, ..., ri+p-1] * 2 - 1

	VMOVUPD	YMM1, [R11+R13*dim+p*1*dim]
	VADDPD	YMM1, YMM1
	VSUBPD	YMM1, [uno]

	VMOVUPD	YMM2, [R11+R13*dim+p*2*dim]
	VADDPD	YMM2, YMM2
	VSUBPD	YMM2, [uno]

	VMOVUPD	YMM3, [R11+R13*dim+p*3*dim]
	VADDPD	YMM3, YMM3
	VSUBPD	YMM3, [uno]

	ADD			R13, p*unroll

	VMOVSD			XMM4, [RDI+48]					; input->stepind
	VSHUFPD			XMM4, XMM4, 0					; input->stepind[j, ..., j+p-1]
	VPERM2F128	YMM4, YMM4, YMM4, 0	; input->stepind[j, ..., j+p-1]

	VMULPD	YMM0, YMM4				; (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]
	VMULPD	YMM1, YMM4
	VMULPD	YMM2, YMM4
	VMULPD	YMM3, YMM4

	MOV		RBX, RDX				; i
	IMUL		RBX, R12					; i*d
	ADD			RBX, RAX				; i*d + j
	IMUL		RBX, dim					; i*d*dim + j*dim

	VMOVUPD	YMM4, [R10+RBX+p*0*dim]	; input->x[i*d*dim + j*dim]
	VADDPD	YMM0, YMM4							; input->x[i*d*dim]

	VMOVUPD	YMM4, [R10+RBX+p*1*dim]
	VADDPD	YMM1, YMM4

	VMOVUPD	YMM4, [R10+RBX+p*2*dim]
	VADDPD	YMM2, YMM4

	VMOVUPD	YMM4, [R10+RBX+p*3*dim]
	VADDPD	YMM3, YMM4

	VMOVAPD	[RCX+RAX*dim+p*0*dim], YMM0
	VMOVAPD	[RCX+RAX*dim+p*1*dim], YMM1
	VMOVAPD	[RCX+RAX*dim+p*2*dim], YMM2
	VMOVAPD	[RCX+RAX*dim+p*3*dim], YMM3

	ADD			RAX, p*unroll
	JMP			.j

.end_j:
	SUB			RAX, p*unroll

.j_no_unroll:
	ADD			RAX, p
	CMP		RAX, R12								; (j<d) ?
	JG			.end_j_no_unroll
	SUB			RAX, p

	VMOVUPD	YMM0, [R11+R13*dim]		; input->r[ri, ..., ri+p-1]
	VADDPD	YMM0, YMM0						; input->r[ri, ..., ri+p-1] * 2
	VSUBPD	YMM0, [uno]						; input->r[ri, ..., ri+p-1] * 2 - 1

	ADD			R13, p

	VMOVSD			XMM4, [RDI+48]					; input->stepind
	VSHUFPD			XMM4, XMM4, 0					; input->stepind[j, ..., j+p-1]
	VPERM2F128	YMM4, YMM4, YMM4, 0	; input->stepind[j, ..., j+p-1]

	VMULPD	YMM0, YMM4				; (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]

	MOV		RBX, RDX					; i
	IMUL		RBX, R12						; i*d
	ADD			RBX, RAX					; i*d + j
	IMUL		RBX, dim						; i*d*dim + j*dim

	VMOVUPD	YMM4, [R10+RBX]		; input->x[i*d*dim + j*dim]
	VADDPD	YMM0, YMM4				; input->x[i*d*dim + j*dim] + (input->r[ri, ..., ri+p-1] * 2 - 1) * input->stepind[j, ..., j+p-1]

	VMOVAPD	[RCX+RAX*dim], YMM0

	ADD			RAX, p
	JMP			.j_no_unroll

.end_j_no_unroll:
	SUB			RAX, p

.j_scalar:
	CMP		RAX, R12
	JGE			.end

	VMOVSD	XMM0, [R11+R13*dim]		; input->r[j]
	VADDSD	XMM0, XMM0						; input->r[ri] * 2
	VSUBSD	XMM0, [uno]						; input->r[ri] * 2 - 1

	INC			R13

	VMOVSD	XMM4, [RAX+48]		; input->stepind
	VMULSD	XMM0, XMM4				; (input->r[ri] * 2 - 1) * input->stepind

	MOV		RBX, RDX					; i
	IMUL		RBX, R12						; i*d
	ADD			RBX, RAX					; i*d + j
	IMUL		RBX, dim						; i*d*dim + j*dim

	VMOVSD	XMM4, [R10+RBX]		; input->x[i*d*dim + j*dim]
	VADDSD	XMM0, XMM4				; input->x[i*d*dim + j*dim] + (input->r[jri * 2 - 1) * input->stepind[j]

	VMOVSD	[RCX+RAX*dim], XMM0

	INC			RAX
	JMP			.j_scalar

.end:
	MOV		[RSI], R13d
	; a questo punto y è stato già modificato e posso terminare

	;
	;	sequenza di uscita dalla funzione
	;

	popaq							; ripristina i registri generali
	MOV		RSP, RBP	; ripristina lo Stack Pointer
	POP			RBP				; ripristina il Base Pointer
	RET								; torna alla funzione C chiamante



