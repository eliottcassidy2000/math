"""Push extensions + the Brouwer/Borsuk-Ulam pivot.
(A) c3 (3-cycles) = Tr(A^3)/3 = spectral moment (grid-triangle count, third face of the dictionary).
(B) The Cayley-Dickson 'i': S skew => eigenvalues i*lambda (imaginary) => the order-4 structure is INTRINSIC.
(C) p=3 mod4 <=> -1 is a QNR <=> negation swaps QR<->QNR FREELY (free Z2, Borsuk-Ulam) <=> Paley TOURNAMENT exists;
    p=1 mod4 <=> -1 is a QR <=> negation FIXES the QR partition (Brouwer-ish) <=> Paley is a GRAPH not a tournament."""
import numpy as np, itertools
# (A) c3 = Tr(A^3)/3
def c3_check(n):
    perms=[list(p) for p in itertools.permutations(range(n))]; pairs=list(itertools.combinations(range(n),2))
    ok=True
    for bits in itertools.product([0,1],repeat=len(pairs)):
        A=np.zeros((n,n),int)
        for (i,j),b in zip(pairs,bits): A[i,j]=1-b; A[j,i]=b
        s=A.sum(1); c3_score=int(sum(1 for a,b,c in itertools.combinations(range(n),3)
            if (A[a,b]+A[b,c]+A[c,a]) in (0,3)))   # cyclic triangle: all-forward or all-back among the 3
        tr=int(round(np.trace(np.linalg.matrix_power(A,3))))
        if c3_score!=tr//3 or tr%3!=0: ok=False; break
    return ok
print("(A) c3 (cyclic 3-cycles = grid-triangle count) == Tr(A^3)/3 (spectral moment) for ALL tournaments?")
for n in [3,4,5]: print(f"    n={n}: {c3_check(n)}")
# (C) QR structure: -1 QR? and negation on QR
def legendre(a,p): 
    a%=p; return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
print("\n(C) p mod4, -1 Legendre, negation on QR:")
for p in [3,5,7,11,13,17,19,23]:
    QR=set(x for x in range(1,p) if legendre(x,p)==1)
    neg_swaps = all(((p-x)%p) not in QR for x in QR)   # negation maps QR entirely OUT of QR?
    neg_fixes = all(((p-x)%p) in QR for x in QR)        # negation keeps QR
    tag = "FREE Z2 (neg swaps QR<->QNR): Borsuk-Ulam, PALEY TOURNAMENT" if (p%4==3) else "FIXED Z2 (neg keeps QR): Brouwer, Paley=GRAPH"
    print(f"    p={p:>2} (={p%4} mod4): (-1|p)={legendre(-1,p):+d}; neg swaps QR<->QNR? {neg_swaps}; keeps QR? {neg_fixes}  => {tag}")
print("\n=> p=3 mod4 <=> -1=QNR <=> negation is a FREE Z2 on QR-structure <=> Paley tournament (self-comp via negation,")
print("   1 fixed vertex=0=the real-axis Cayley eigenvalue) => the Borsuk-Ulam / three-pillars HARD regime (flip-rank obstruction).")
