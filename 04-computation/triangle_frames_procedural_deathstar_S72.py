# Procedurally analyze extension frames of the triangular numbers T_n=n(n+1)/2, along uniform axes.
from math import comb, factorial, log
from fractions import Fraction as Fr
def diff_degree(seq):
    # returns polynomial degree if seq is eventually polynomial, else None
    s=seq[:]
    for d in range(len(seq)):
        if all(x==0 for x in s[1:]): return d-1 if d>0 else 0
        if len(set(s))==1: return d
        s=[s[i+1]-s[i] for i in range(len(s)-1)]
        if len(s)<2: break
    return None
def growth(seq):
    if len(seq)<4: return "?"
    d=diff_degree(seq)
    if d is not None: return f"polynomial deg {d}"
    r=[seq[i+1]/seq[i] for i in range(len(seq)-1) if seq[i]]
    if r and abs(r[-1]-r[-2])<0.05*r[-1]: return f"exponential ratio~{r[-1]:.2f}"
    return "factorial/super-exponential (ratio ->inf)"
def parity(seq): 
    p=[x%2 for x in seq]; 
    return "all odd" if all(p) else ("all even" if not any(p) else f"mixed {p[:8]}")
N=range(1,9)
frames={
 "T triangular (T_n=C(n+1,2))":            [n*(n+1)//2 for n in N],
 "ADD polygonal s=4 (squares)":            [n*n for n in N],
 "ADD polygonal s=5 (pentagonal)":         [n*(3*n-1)//2 for n in N],
 "ADD Faulhaber p=2 (sq pyramidal Sum k^2)":[sum(k*k for k in range(1,n+1)) for n in N],
 "ADD Faulhaber p=3 (Sum k^3 = T_n^2)":    [sum(k**3 for k in range(1,n+1)) for n in N],
 "DIM simplicial d=3 (tetrahedral C(n+2,3))":[comb(n+2,3) for n in N],
 "DIM simplicial d=4 (C(n+3,4))":          [comb(n+3,4) for n in N],
 "MUL factorial (Prod k = n!)":            [factorial(n) for n in N],
 "MUL Ham-path monoid (spectrum odds\\{7,21})": [1,3,5,9,11,13,15,17],  # ordered achieved values, multiplicative
 "FAC arborescence min ((n-1)!)":          [factorial(n-1) for n in N],
 "EXP tilings (2^C(n-1,2))":               [2**comb(n-1,2) for n in N],
 "CAT triangulations (Catalan C_{n-1})":   [comb(2*(n-1),n-1)//n for n in N],
}
print(f"{'FRAME':46} {'first terms':32} {'growth':30} parity")
print("-"*130)
for name,seq in frames.items():
    print(f"{name:46} {str(seq[:6]):32} {growth([float(x) for x in seq]):30} {parity(seq)}")
print("\n--- composition-mode ladder (how each BUILDS from the triangle substrate) ---")
print(" ADD (+): sum of k terms  -> Fermat/Cauchy: every n = sum of <=s s-gonal (representability theorem)")
print(" MUL (x): product of terms -> Ham monoid: spectrum = odds\\{7,21}, closed under x (the MULT analog of Fermat)")
print(" FAC (!): branching        -> arborescence (n-1)!-bands; simplicial C(n+d-1,d) binomial branching")
print(" EXP (^): 2^substrate      -> 2^{T_{n-2}} tilings (the hypercube the tournament invariants live on)")
