.global asm_add
.global asm_sub
.global asm_mul
.intel_syntax noprefix

# Addition
# RDI : pointer vers dest (aussi retour)
# RSI : pointer vers src //
# RDX : longeur des nombes (index)
# RAX : Retour le reste 
asm_add:
    clc
    test rdx,rdx
    jz .fin_add

.boucle_add:
    mov rax, [rsi] 
    adc [rdi], rax #des[i] = des[i] + src[i]

    lea rdi, [rdi + 8] #Pointer au suivant 
    lea rsi, [rsi + 8] #Pointer au suivant

    dec rdx # index --
    jnz .boucle_add # while index > 0

    #Finalment o retouner Carry
    mov rax, 0
    adc rax, 0
    ret
.fin_add:
    xor rax, rax
    ret 

# Soustraction 
# Avec toujour dest > src
# RDI : pointer vers dest (aussi retour)
# RSI : pointer vers src //
# RDX : longeur des nombes (index)
# RAX : Retour le reste 
asm_sub:
    clc
    test rdx, rdx
    jz .fin_sub

.boucle_sub:
    mov rax, [rsi] 
    sbb [rdi], rax #des[i] = des[i] - src[i]

    lea rdi, [rdi + 8] #Pointer au suivant 
    lea rsi, [rsi + 8] #Pointer au suivant

    dec rdx # index --
    jnz .boucle_sub # while index > 0

    #Finalment 
    mov rax, 0
    adc rax, 0
    ret 
.fin_sub:
    xor rax, rax
    ret

# Multiplication par entier A=n*B
# RDI : pointer vers dest (aussi retour)
# RSI : pointer vers src //
# RDX : Scalaire n 
# RCX : longeur des nombes (index)
# RAX : Retour le reste 
asm_mul:
    clc
    test rcx, rcx
    jz .fin_mul

    mov r9, rdx
    xor r8, r8
    #commencer le calcul
.boucle_mul:

    # Multiplication
    mov rax, [rsi]
    mul r9

    # Ajouter le rest
    add rax, r8 
    adc rdx, 0

    # mettre les resultat en ces places
    add [rdi], rax
    adc rdx, 0
    mov r8, rdx

    # Decaller les pointer 
    lea rsi, [rsi + 8]
    lea rdi, [rdi + 8]

    #decrementer index
    dec rcx 
    jnz .boucle_mul
    mov rax, r8
    ret
.fin_mul:
    xor rax, rax
    ret






