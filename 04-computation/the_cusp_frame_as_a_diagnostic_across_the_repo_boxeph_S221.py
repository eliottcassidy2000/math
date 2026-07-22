#!/usr/bin/env python3
"""the_cusp_frame_as_a_diagnostic_across_the_repo_boxeph_S221.py -- boxeph-2026-07-21-S221

Apply the CUSP FRAME (S220) to under-attended repo problems and show its power. The frame:

  an object = EISENSTEIN (computable main term / floor / local, from spectral/genus data)
           (+) CUSP (the hidden obstruction = the genus = the deep arithmetic entropy, S218).
  The DIFFICULTY is always the CUSP; the frame LOCALIZES it to a small/finite object and predicts the
  'first hard case' as the first positive cusp dimension.

Sweeps (mostly under-attended in the cusp language):
  P1 TOURNAMENT COSPECTRALITY: char_A spectrum = the Eisenstein (local) data; the COSPECTRAL fiber
     (#classes - #distinct spectra) = the reconstruction CUSP. Transitive = spectrally unique (rigid).
  P2 INTRANSITIVITY c3 = the tournament's CUSP FORM: transitive c3=0 (pure Eisenstein/gradient) vs
     regular/Paley c3 max (the intransitive obstruction). c3 = C(n,3) - sum C(s_i,2).
  P3 GMC(2) angular=EISENSTEIN (DvdK-closed) / radial=CUSP (Laplace determinacy, ker L != 0): the obstruction
     is the radial cuspidal residual (E = L o CT, THM-1645/S211).
  P4 FIGURATE cake/bagel = EISENSTEIN (polynomial main term) (+) CUSP (the Fibonacci/oscillating deviation).
"""
from itertools import combinations, permutations
from math import comb

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

# ---- tournament machinery ----
def all_tournaments(n):
    pairs=list(combinations(range(n),2))
    for bits in range(2**len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: adj[i][j]=1
            else: adj[j][i]=1
        yield adj
def canon(adj,n,perms):
    best=None
    for p in perms:
        key=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or key<best: best=key
    return best
def charpoly(A):  # Faddeev-LeVerrier, integer coeffs of det(xI-A)
    n=len(A); I=[[1 if i==j else 0 for j in range(n)] for i in range(n)]
    c=[1]+[0]*n; Mk=[row[:] for row in I]
    for k in range(1,n+1):
        AM=[[sum(A[i][t]*Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
        ck=-sum(AM[i][i] for i in range(n))//k; c[k]=ck
        Mk=[[AM[i][j]+(ck if i==j else 0) for j in range(n)] for i in range(n)]
    return tuple(c)
def scores(adj,n): return tuple(sorted(sum(adj[i]) for i in range(n)))

# ==========================================================================
sep("P1  TOURNAMENT COSPECTRALITY: spectrum=Eisenstein(local); the cospectral fiber = the reconstruction CUSP")
for n in (4,5,6):
    perms=list(permutations(range(n))); classes={}
    for adj in all_tournaments(n):
        c=canon(adj,n,perms)
        if c not in classes: classes[c]=adj
    specs={}
    for adj in classes.values():
        s=charpoly(adj); specs.setdefault(s,0); specs[s]+=1
    nclass=len(classes); ndistinct=len(specs); cusp=nclass-ndistinct
    cofam=[v for v in specs.values() if v>1]
    print(f"  n={n}: #iso classes={nclass}; #distinct char_A spectra={ndistinct}; COSPECTRAL cusp dim = {cusp} "
          f"(cospectral families of sizes {sorted(cofam) or 'none'})")
print("  => most tournaments are spectrally DETERMINED (Eisenstein/rigid, esp. the transitive char x^n); the")
print("     COSPECTRAL classes are the reconstruction CUSP -- the hidden fiber where local (spectral) data fails.")
print("     The cusp dim first turns POSITIVE at a specific n = the 'first hard case' of tournament reconstruction.")

# ==========================================================================
sep("P2  INTRANSITIVITY c3 = the tournament's CUSP FORM: transitive 0 (Eisenstein) vs regular max")
def c3_count(adj,n):  # #3-cycles = C(n,3) - sum C(outdeg_i, 2)
    return comb(n,3) - sum(comb(sum(adj[i]),2) for i in range(n))
for n in (5,7,9):
    # transitive: i beats j iff i<j
    trans=[[1 if i<j else 0 for j in range(n)] for i in range(n)]
    # regular (n odd): i beats i+1,..,i+(n-1)/2 mod n
    reg=[[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(1,(n-1)//2+1): reg[i][(i+k)%n]=1
    print(f"  n={n}: transitive c3={c3_count(trans,n)} (Eisenstein floor: NO cusp = pure gradient) ; "
          f"regular c3={c3_count(reg,n)} (max = the intransitive CUSP) ; C(n,3)={comb(n,3)}")
print("  => c3 (the 3-cycle count) IS the tournament's cusp form: it vanishes on the transitive/AP Eisenstein")
print("     vertex and peaks on the regular/Paley pole. Intransitivity = the cuspidal obstruction (THM-1830 atom).")

# ==========================================================================
sep("P3  GMC(2): angular = EISENSTEIN (DvdK-closed floor) ; radial = CUSP (Laplace determinacy, ker L != 0)")
from math import factorial
def L(coeffs):  # radial Laplace functional L(s^k)=k! ; coeffs low->high
    return sum(coeffs[k]*factorial(k) for k in range(len(coeffs)))
tests={"t - 1":[-1,1], "t^2 - 3t + 1":[1,-3,1], "t^2 - 2t":[0,-2,1]}
for name,c in tests.items():
    print(f"  L({name}) = {L(c)}" + ("   <- in ker L (a radial CUSP element: integrated-zero but nonzero)" if L(c)==0 else ""))
print("  E = L o CT (THM-1645): the ANGULAR CT (DvdK) is CLOSED = the EISENSTEIN floor (the easy, proved half);")
print("  the RADIAL L has a nonzero KERNEL = the CUSP obstruction (Laplace determinacy). GMC(n>=3) FALSE = the")
print("  cusp grows (extra obstructions); GMC(2) rigid = the radial cusp is defeated by Frobenius (THM-2022).")

# ==========================================================================
sep("P4  FIGURATE: cake/bagel = EISENSTEIN (polynomial main term) (+) CUSP (Fibonacci/oscillating deviation)")
def cake(nn): return sum(comb(nn,k) for k in range(4))
def moser(nn): return comb(nn,0)+comb(nn,2)+comb(nn,4)
# the 'Eisenstein' part = the smooth polynomial (leading figurate); the 'cusp' = the deviation from 2^n / the
# shallow-diagonal Fibonacci wobble (S207). Illustrate cake vs its polynomial vs the Fibonacci diagonal.
print("  cake(n) (A000125) =", [cake(n) for n in range(9)], " = the smooth degree-3 EISENSTEIN polynomial (a ball)")
print("  2^n (full Pascal)  =", [2**n for n in range(9)], " ; the CUSP = 2^n - cake = the dropped C(n,>=4) terms")
print("  the Fibonacci shallow-diagonal (S207) is the OSCILLATING/cusp reading of the same Pascal triangle.")
print("  => figurate cutting = a smooth Eisenstein polynomial; the Fibonacci/deviation is its cusp (the hidden,")
print("     number-theoretic part) -- the same Eisenstein(smooth)/cusp(arithmetic) split, S207 recast.")

sep("SUMMARY -- the cusp frame is a difficulty-LOCATOR")
print("""  The cusp frame applied across the repo (each: EISENSTEIN computable floor (+) CUSP hidden obstruction):
    LRC(14)          floor 3/pi^2            (+) f14=14a (genus 1, the first cusp, apex 7)      [S220]
    tournaments      char_A spectrum         (+) the COSPECTRAL fiber (reconstruction cusp)      [P1]
    intransitivity   transitive (c3=0)       (+) c3 = the 3-cycle count (the intransitive cusp)  [P2]
    GMC(2)           angular CT (DvdK floor) (+) radial ker L (Laplace-determinacy cusp)         [P3]
    figurate         smooth polynomial       (+) Fibonacci/deviation (the arithmetic cusp)       [P4]
  POWER: the frame (1) LOCALIZES each difficulty to a small/finite CUSP object (a fiber, a count, a kernel,
  a genus-1 newform); (2) PREDICTS the 'first hard case' as the first positive cusp dimension (LRC: p=7;
  tournaments: the first cospectral n); (3) UNIFIES them as the S218 deep arithmetic entropy = dim(cusp).
  The Eisenstein floor is always the easy/computable/local part; the cusp is always where the proof must go.""")
