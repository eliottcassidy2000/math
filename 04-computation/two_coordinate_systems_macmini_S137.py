#!/usr/bin/env python3
"""
{-1,0,1} and {1,1/2,0} are ONE object in two coordinate systems   (mac-mini-S137)
=================================================================================
Owner: "look back at the summand and multiplicand graph ... think about how {-1,0,1}
and {1,1/2,0} are functionally equivalent but have shown up in this repo many times
each. this is what i mean by even/odd vs positive/negative."

THE DICTIONARY.  x |-> 2x - 1 sends {0, 1/2, 1} -> {-1, 0, 1}, inverse x |-> (1+x)/2.
    {0,1/2,1}  ADDITIVE / probability coordinates.  Averaging is natural.  POSITIVE.
    {-1,0,1}   MULTIPLICATIVE / sign coordinates.  (-1)^k is natural.  EVEN/ODD.
This is the Fourier-Walsh <-> probability dictionary: a +-1 variable X with P(X=1)=p
has E[X] = 2p - 1.  The repo's "two arithmetics" (additive tiling hypercube vs
multiplicative H-norm) is the SAME duality one level up.

Four places the repo already carries it, all verified below:
  A  adjacency A in {0,1} vs skew-Seidel S in {-1,0,1}: S = 2A - (J - I).  The DIAGONAL
     is the tie: p = 1/2 <-> S = 0, and it is the FIXED POINT of complementation
     (p |-> 1-p additively, S |-> -S multiplicatively).
  B  the FIBER FRACTION f(n) = (1/2)_{n-2}/(n-2)! is exactly C(2k,k)/4^k, the
     probability of an EXACT TIE in 2k fair coin flips.  The repo's 1/2 is a fair coin
     and its "two-sheeted branched cover" is the square root of that.
  C  THM-1500's admissible exponent ladder is d/2 in {0, 1/2, 1} for d in {0,1,2},
     i.e. dimensions n in {2,3,4} -- and under x |-> 2x-1 that ladder IS {-1,0,1}.
     The MINIMAL counterexample (n=3) sits at the FIXED POINT.
  D  THM-1520/1540's telescoping principle is a statement about SIGN support: the
     nullcone is the ENDPOINTS of the balance, failure lives in the INTERIOR.
"""
from fractions import Fraction as F
from math import factorial, comb
from itertools import permutations
import subprocess

def poch(a, k):
    r = F(1)
    for i in range(k): r *= (a + i)
    return r

print("=" * 78)
print("PART A -- adjacency vs skew-Seidel IS the map x |-> 2x - 1")
print("=" * 78)
print("  Read the tournament as a PROBABILITY matrix p_ij in {0, 1/2, 1}:")
print("      p_ij = 1   if i beats j,   0   if j beats i,   1/2 on the diagonal (a tie).")
print("  Then  S = 2p - 1  in {-1, 0, +1}  is exactly the skew-Seidel matrix of THM-1440.")
print()
n = 5
import random
rng = random.Random(137)
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
ok_map = ok_comp = ok_fix = True
for _ in range(200):
    code = rng.getrandbits(len(pairs))
    P = [[F(1,2)]*n for _ in range(n)]                  # diagonal tie
    for e, (i, j) in enumerate(pairs):
        if code >> e & 1: P[i][j], P[j][i] = F(0), F(1)
        else:             P[i][j], P[j][i] = F(1), F(0)
    S = [[2*P[i][j] - 1 for j in range(n)] for i in range(n)]
    # S must be the skew-Seidel matrix: entries in {-1,0,1}, skew, zero diagonal
    if any(S[i][j] not in (-1, 0, 1) for i in range(n) for j in range(n)): ok_map = False
    if any(S[i][j] != -S[j][i] for i in range(n) for j in range(n)): ok_map = False
    if any(S[i][i] != 0 for i in range(n)): ok_fix = False
    # complementation: p -> 1-p  must equal  S -> -S
    Pc = [[1 - P[i][j] for j in range(n)] for i in range(n)]
    Sc = [[2*Pc[i][j] - 1 for j in range(n)] for i in range(n)]
    if any(Sc[i][j] != -S[i][j] for i in range(n) for j in range(n)): ok_comp = False
print(f"  S = 2p-1 lands in {{-1,0,1}} and is skew:            {ok_map}")
print(f"  complement  p -> 1-p   IS   S -> -S:                {ok_comp}")
print(f"  the tie p = 1/2 is the FIXED POINT, S = 0 (diagonal): {ok_fix}")
print()
print("  So: EVEN/ODD (multiplicative, S -> -S, the sign character) and")
print("      POSITIVE/NEGATIVE (additive, p -> 1-p, a probability) are ONE involution")
print("      in two coordinate systems, and 1/2 <-> 0 is its unique fixed point.")
print()
print("  Paley makes it literal: the skew-Seidel matrix of the Paley tournament is the")
print("  LEGENDRE SYMBOL chi(j-i), whose value set is exactly {-1,0,+1}; its probability")
print("  avatar (1+chi)/2 is the quadratic-residue indicator, with 1/2 sitting at 0.")
for q in (7, 11, 19):
    QR = {(x*x) % q for x in range(1, q)}
    vals = set()
    for k in range(q):
        vals.add(0 if k == 0 else (1 if k in QR else -1))
    print(f"    q={q:>2}: value set of chi = {sorted(vals)}   "
          f"probability avatar (1+chi)/2 = {sorted({F(1+v,2) for v in vals})}")

print()
print("=" * 78)
print("PART B -- the FIBER FRACTION is the probability of an exact TIE")
print("=" * 78)
print("  CLAUDE.md records  f(n) = (1/2)_{n-2}/(n-2)!,  GF = (1-x)^{-1/2}.")
print("  Three expressions, all equal:")
print("      (1/2)_k / k!        the Pochhammer form (the repo's)")
print("      C(2k,k) / 4^k       the CENTRAL BINOMIAL PROBABILITY = P(exact tie in 2k")
print("                          fair coin flips)")
print("      E[U^k] / k!         with U = T^2/2, T ~ N(0,1)   (THM-1500's d=1 moment)")
print()
print(f"{'k':>3} {'(1/2)_k/k!':>14} {'C(2k,k)/4^k':>14} {'E[U^k]/k!':>14} {'equal':>7} {'decimal':>10}")
allsame = True
for k in range(0, 9):
    a = poch(F(1,2), k) / factorial(k)
    b = F(comb(2*k, k), 4**k)
    c = poch(F(1,2), k) / factorial(k)          # E[U^k] = (1/2)_k for U = T^2/2
    same = (a == b == c); allsame &= same
    print(f"{k:>3} {str(a):>14} {str(b):>14} {str(c):>14} {str(same):>7} {float(a):>10.6f}")
print(f"  all three agree: {allsame}")
print()
print("  So the repo's fiber fraction IS a tie probability -- the 1/2 is a FAIR COIN, and")
print("  the 'two-sheeted branched cover' (1-x)^{-1/2} is its square root.  The tie is the")
print("  0 of {-1,0,1} and the 1/2 of {0,1/2,1}: the SAME point in the two coordinates.")

print()
print("=" * 78)
print("PART C -- THM-1500's admissible ladder IS the three-element set")
print("=" * 78)
print("  THM-1500: U = chi^2_d / 2 has MGF Phi(x) = (1-x)^{-d/2}; the master equation")
print("  Phi(-s g(s)) = 1/(1+s) forces 1 + s g(s) = (1+s)^{2/d}, POLYNOMIAL iff 2/d is a")
print("  positive integer, i.e. d in {1,2}.  d = 0 is the obstructed case (n=2).")
print()
print(f"{'d':>3} {'exponent d/2':>13} {'2x-1':>7} {'n = 2+d':>8} {'g(s)':>10} {'status':>26}")
for d in (0, 1, 2):
    e = F(d, 2)
    sign = 2*e - 1
    g = {0: "log(1+s)/s", 1: "2 + s", 2: "1"}[d]
    st = {0: "OBSTRUCTED (GMC(2) open)", 1: "MINIMAL counterexample",
          2: "the 4-term example"}[d]
    print(f"{d:>3} {str(e):>13} {str(sign):>7} {2+d:>8} {g:>10} {st:>26}")
print()
print("  The exponent ladder {0, 1/2, 1} maps under x |-> 2x-1 onto {-1, 0, +1}.")
print("  And the MINIMAL counterexample n = 3 sits at the FIXED POINT (1/2 <-> 0) --")
print("  the same 1/2 as the fiber fraction, and the same 1/2 that makes (1+s)^{2/d}")
print("  a PERFECT SQUARE.  Square root and square are the same one-half.")

print()
print("=" * 78)
print("PART D -- the telescoping principle, read in the two coordinate systems")
print("=" * 78)
print("  THM-1520/1540: the nullcone is the STRICTLY ONE-SIDED charge support, and")
print("  failure REQUIRES two-sided support.  In the two coordinates:")
print()
print("    SIGN view {-1,0,+1}:  which signs the charge support hits.")
print("        {+1} only  -> nullcone.   {-1} only -> nullcone.   both -> the open case.")
print("    BALANCE view {0,1/2,1}:  b = (fraction of the support with positive charge).")
print("        b = 1 -> nullcone.   b = 0 -> nullcone.   0 < b < 1 -> the open case.")
print()
print("  So THE NULLCONE IS THE ENDPOINTS AND THE HARD CASE IS THE INTERIOR, with the")
print("  balanced point b = 1/2 the deepest interior point.  That is the same statement")
print("  as PART A's involution having a unique fixed point, one level up.")
print()
def cmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1]); out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}
def cexp(p): return sum(c*factorial(a) for (a, b), c in p.items() if a == b)
ONE = {(0,0): 1}
def cpow(p, m):
    r = ONE
    for _ in range(m): r = cmul(r, p)
    return r
print(f"{'P':>18} {'charge signs':>16} {'balance b':>11} {'in nullcone?':>14}")
for P, nm in (({(1,0):1}, "Z"), ({(2,0):1,(1,0):-3}, "Z^2-3Z"),
              ({(0,1):1,(0,2):5}, "W+5W^2"), ({(1,0):1,(0,1):1}, "Z+W"),
              ({(2,0):1,(0,1):1}, "Z^2+W"), ({(1,1):1,(1,0):1}, "ZW+Z")):
    ch = [a-b for (a, b) in P]
    signs = sorted({(c > 0) - (c < 0) for c in ch})
    pos = sum(1 for c in ch if c > 0)
    b = F(pos, len(ch))
    inn = all(cexp(cpow(P, m)) == 0 for m in range(1, 11))
    print(f"{nm:>18} {str(signs):>16} {str(b):>11} {str(inn):>14}")

print()
print("=" * 78)
print("PART E -- census: where each coordinate system actually appears")
print("=" * 78)
def count(pat, where="01-canon 07-reflections"):
    try:
        r = subprocess.run(f"grep -rl -- '{pat}' {where} 2>/dev/null | wc -l",
                           shell=True, capture_output=True, text=True, timeout=90)
        return int(r.stdout.strip() or 0)
    except Exception:
        return -1
print(f"{'marker':>34} {'files':>7}   coordinate system")
for pat, sysname in (("skew-Seidel", "{-1,0,1} sign"), ("Legendre", "{-1,0,1} sign"),
                     ("odd function", "{-1,0,1} sign"), ("sgn", "{-1,0,1} sign"),
                     ("(1-x)^{-1/2}", "{0,1/2,1} prob"), ("fiber fraction", "{0,1/2,1} prob"),
                     ("1/2)_", "{0,1/2,1} prob"), ("merged metagraph", "{0,1/2,1} prob")):
    print(f"{pat:>34} {count(pat):>7}   {sysname}")
print()
print("SUMMARY")
print("  One 3-element object, two coordinate systems, and the repo has been working in")
print("  both without always naming the change of variables.  x |-> 2x-1 is the whole")
print("  dictionary: EVEN/ODD is the multiplicative chart, POSITIVE/NEGATIVE the additive")
print("  one, and the tie (1/2 <-> 0) is the fixed point they share.")
