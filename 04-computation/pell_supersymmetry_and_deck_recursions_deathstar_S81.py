#!/usr/bin/env python3
"""
pell_supersymmetry_and_deck_recursions_deathstar_S81.py

(1) CHASE THE PREDICTION: the a/b Pell identity E_n^2 - O_n^2 = (x^2-1)^n
    (THM-1880) has an exact GMC(2) MOMENT shadow. With a<->Z, abar<->W,
    x^2-1 = a*abar <-> ZW = s (radial), and E the Gaussian functional
    (E[Z^a W^b] = a! delta_ab):
        sym_n = (Z^n + W^n)/2   (bosonic/even),  alt_n = (Z^n - W^n)/2  (fermionic/odd)
        PREDICT:  E[sym_n^2] - E[alt_n^2] = E[(ZW)^n] = n!
    and generally E[sym(P)^2 - alt(P)^2] = E[P * Ptilde] (Ptilde = charge-conjugate,
    Z<->W). Verified here by exact Wick (symbolic charge bookkeeping).

(2) THE RECURSION MODES as subtournament decks: interpret A..G as vertex-deleted
    subtournaments; test the sign patterns ++- (n=3), +++--- (n=6), ++-+--+ (n=7)
    as signed deck-sums of invariants (signed Redei R = mac-mini THM-1936, c3, H),
    and report which relations hold / which classes are built from smaller ones.
"""
from math import comb, factorial
from itertools import combinations
from fractions import Fraction as Fr

def sep(t): print("\n"+"="*68+"\n"+t+"\n"+"="*68)

# ---------- (1) GMC(2) Wick functional on polynomials in Z,W ----------
# represent a polynomial as dict {(a,b): coeff} meaning coeff * Z^a W^b
def mul(p,q):
    r={}
    for (a,b),c in p.items():
        for (a2,b2),c2 in q.items():
            k=(a+a2,b+b2); r[k]=r.get(k,0)+c*c2
    return {k:v for k,v in r.items() if v}
def add(p,q):
    r=dict(p)
    for k,v in q.items(): r[k]=r.get(k,0)+v
    return {k:v for k,v in r.items() if v}
def scal(p,s): return {k:v*s for k,v in p.items()}
def E(p):
    # E[Z^a W^b] = a! if a==b (charge 0) else 0
    return sum(v*factorial(a) for (a,b),v in p.items() if a==b)
def conj(p):  # charge-conjugate: swap Z<->W
    return {(b,a):v for (a,b),v in p.items()}

sep("(1) PELL SUPERSYMMETRY for GMC(2): E[sym^2] - E[alt^2] = E[(ZW)^n] = n!")
print(f"{'n':>2} {'E[sym^2]':>10} {'E[alt^2]':>10} {'diff':>8} {'n!':>8} {'= E[(ZW)^n]?':>12}")
for n in range(1,8):
    Zn={(n,0):Fr(1)}; Wn={(0,n):Fr(1)}
    sym=scal(add(Zn,Wn),Fr(1,2)); alt=scal(add(Zn,scal(Wn,-1)),Fr(1,2))
    es=E(mul(sym,sym)); ea=E(mul(alt,alt)); diff=es-ea
    radial=E({(n,n):Fr(1)})   # E[(ZW)^n] = E[Z^n W^n] = n!
    ok = (diff==factorial(n)==radial)
    print(f"{n:>2} {str(es):>10} {str(ea):>10} {str(diff):>8} {factorial(n):>8} {str(ok):>12}")

sep("(1b) GENERAL: E[sym(P)^2 - alt(P)^2] = E[P * conj(P)] for random-ish P")
# P with a couple of charge modes
tests=[{(2,0):Fr(1),(0,1):Fr(3),(1,1):Fr(-2)},
       {(3,0):Fr(1),(1,2):Fr(1),(0,0):Fr(5)},
       {(2,1):Fr(1),(0,3):Fr(-1)}]
for P in tests:
    sym=scal(add(P,conj(P)),Fr(1,2)); alt=scal(add(P,scal(conj(P),-1)),Fr(1,2))
    lhs=E(mul(sym,sym))-E(mul(alt,alt)); rhs=E(mul(P,conj(P)))
    print(f"  P={P}\n     E[sym^2-alt^2]={lhs}   E[P*conj(P)]={rhs}   equal? {lhs==rhs}")
print("  => the Pell/Chebyshev E_n^2-O_n^2=(x^2-1)^n is the polynomial lift of the")
print("     GMC identity E[sym^2-alt^2]=E[P*Ptilde]; on P=Z^n it localizes on (ZW)^n=n!.")

# ---------- (2) signed subtournament deck recursions ----------
sep("(2) RECURSION MODES as signed subtournament decks (++- / +++--- / ++-+--+)")
def all_t(n):
    P=list(combinations(range(n),2))
    for bits in range(1<<len(P)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(P):
            if bits>>k&1: A[i][j]=1
            else: A[j][i]=1
        yield A
def sub(A,v,n):  # delete vertex v
    idx=[u for u in range(n) if u!=v]
    return [[A[i][j] for j in idx] for i in idx]
def c3(A,n):
    s=[sum(r) for r in A]; return comb(n,3)-sum(comb(x,2) for x in s)
def signed_redei(A,n):  # R(T)=sum_{ham paths pi} sgn(pi); via permanent-like signed DP over perms
    from itertools import permutations
    R=0
    for pi in permutations(range(n)):
        if all(A[pi[k]][pi[k+1]] for k in range(n-1)):
            # sgn of the permutation pi (as sequence)
            inv=sum(1 for i in range(n) for j in range(i+1,n) if pi[i]>pi[j])
            R+=(-1)**inv
    return R
PATT={3:[1,1,-1], 6:[1,1,1,-1,-1,-1], 7:[1,1,-1,1,-1,-1,1]}
for n in (3,6,7):
    if n==7:
        print(f"  n=7: skipped (2^21 tournaments); pattern ++-+--+ noted"); continue
    hits={'c3':0,'R':0,'H':0}; total=0
    for A in all_t(n):
        total+=1
        decks=[sub(A,v,n) for v in range(n)]
        # signed deck sum with pattern (sorted by an intrinsic order: by deck's own invariant to be canonical)
        for inv,name in [(lambda M,m: c3(M,m),'c3'),(lambda M,m: signed_redei(M,m),'R')]:
            vals=sorted(inv(M,n-1) for M in decks)   # sort for a canonical signed assignment
            s=sum(p*x for p,x in zip(PATT[n],vals))
            whole=inv(A,n)
            if s==whole or s==0: hits[name]+=1
    print(f"  n={n}: signed deck-sum (pattern {PATT[n]}) matches whole-or-vanishes: "
          f"c3 {hits['c3']}/{total}, R {hits['R']}/{total}")
print("""  READING: exploratory -- tests whether the given sign patterns, applied to the
  (sorted) vertex-deleted subtournament invariants, reproduce the whole tournament's
  invariant or vanish. mac-mini THM-1936: R is JOIN-multiplicative (build-up), the
  cleaner 'smaller classes -> larger class' law; deletion signed-sums are the dual.""")
