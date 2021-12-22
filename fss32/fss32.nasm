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

section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	ris	resd	1		; type* ris: puntatore utilizzato per i valori restituiti dalle funzioni

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

	dim	equ	4		; dimensione in byte di un singolo dato (4 se float, 8 se double)
	p	equ	4		; grado di parallelismo SIMD (4 se float, 2 se double)

	x	equ	8		; puntatore al primo vettore
	y	equ	12		; puntatore al secondo vettore
	n	equ	16		; dimensione vettori

prodottoScalare:
	;
	; sequenza di ingresso nella funzione
	;

	push	ebp				; salvo il Base Pointer
	mov		ebp, esp			; il Base Pointer punta al record di attivazione corrente
	push	ebx				; salvo i registri da preservare
	push	esi
	push	edi

	;
	; lettura dei parametri dal record di attivazione
	;

	mov		eax, [ebp+x]		; x
	mov		ebx, [ebp+y]		; y
	mov		ecx, [ebp+n]		; n

	;
	; corpo della funzione
	;

	xorps	xmm0, xmm0		; ps = 0
	mov 	esi, 0			; i = 0

for_i:
	mov 	edi, esi			; indTemp = i
	add 	edi, p			; indTemp+=p
	cmp 	edi, ecx			; (indTemp > n) ?
	jg		for_i_scalar		; se vero passa a lavorare con scalari, anziché vettori

	movaps	xmm1, [eax+esi*dim]	; x[i, ..., i+p-1]
	mulps 	xmm1, [ebx+esi*dim]	; temp[i, ..., i+p-1] = x[i, ..., i+p-1] * y[i, ..., i+p-1]
	addps	xmm0, xmm1			; ps[i, ..., i+p-1] += temp[i, ..., i+p-1]

	add		esi, p			; i+=p
	jmp		for_i

for_i_scalar:
	cmp		esi, ecx			; (i >= n) ?
	jge		end

	movss	xmm1, [eax+esi*dim]	; x[i]
	mulss	xmm1, [ebx+esi*dim]	; temp[i] = x[i] * y[i]
	addss	xmm0, xmm1			; ps[i, ..., i+p-1] += temp[i]

	inc		esi				; i++
	jmp		for_i_scalar

end:
	haddps	xmm0, xmm0		; effettuo le due somme orizzontali rimanenti
	haddps	xmm0, xmm0

	movss	[ris], xmm0		; *ris = ps

	mov		eax, ris			; return ris

	;
	;	sequenza di uscita dalla funzione
	;

	pop		edi				; ripristina i registri da preservare
	pop		esi
	pop		ebx
	mov		esp, ebp			; ripristina lo Stack Pointer
	pop		ebp				; ripristina il Base Pointer
	ret						; ritorna alla funzione chiamante
