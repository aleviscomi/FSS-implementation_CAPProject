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

global prodottoScalare

	x	equ	8		; puntatore al vettore dei coefficienti
	y	equ	12		; puntatore al vettore x
	n	equ	16		; dimensione vettori
	ris	equ	20		; puntatore alla variabile contenente il risultato

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

	MOV	EAX, [EBP+x]		; x
	MOV	EBX, [EBP+y]		; y
	MOV	ECX, [EBP+n]		; n
	MOV	EDX, [EBP+ris]		; ris

	;
	; corpo della funzione
	;

	XORPS	XMM0, XMM0		; ps = 0
	MOV 	ESI, 0			; i = 0

for_i:
	MOV 	EDI, ESI			; indTemp = i
	ADD 	EDI, p			; indTemp+=p
	CMP 	EDI, ECX			; (indTemp > n) ?
	JG		for_i_scalar		; se vero passa a lavorare con scalari, anziché vettori

	MOVAPS	XMM1, [EAX+ESI*dim]	; x[i, ..., i+p-1]
	MULPS 	XMM1, [EBX+ESI*dim]	; temp[i, ..., i+p-1] = x[i, ..., i+p-1] * y[i, ..., i+p-1]
	ADDPS	XMM0, XMM1			; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	ADD		ESI, p			; i+=p
	JMP		for_i

for_i_scalar:
	CMP	ESI, ECX			; (i >= n) ?
	JGE		end

	MOVSS	XMM1, [EAX+ESI*dim]	; x[i]
	MULSS	XMM1, [EBX+ESI*dim]	; temp[i] = x[i] * y[i]
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