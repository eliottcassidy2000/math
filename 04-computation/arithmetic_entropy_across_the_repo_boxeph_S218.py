#!/usr/bin/env python3
"""arithmetic_entropy_across_the_repo_boxeph_S218.py -- boxeph-2026-07-21-S218

HISTORICAL / REFUTED SYNTHESIS (MISTAKE-231).  The claimed repo-wide entropy
invariant is retracted.  The script mixes a finite-fiber Hartley count, score
statistics, a continued-fraction digit mean, and moment depth; these are not
one invariant.  It intentionally retains the original calculations so their
exact survivors and counterexamples remain reproducible.  Printed conclusions
below are historical unless the status banner identifies them as retained.

Extend "arithmetic entropy" (S217) and apply it across the repo. Unified definition:

  H_arith(X | L) = log2 | { X' : L(X') = L(X) } |   -- the bits of a GLOBAL object X HIDDEN from its
  LOCAL invariant L(X). Zero entropy = local determines global (RIGID); positive = hidden global info.
  The repo's difficulty lives exactly where this entropy is positive; its rigid extrema are the zeros.

Pillars (each an instance of one local->global information deficit):
  P1 BINARY FORMS | genus: h(D) = genera x deep; the DEEP (within-genus) part is the hidden entropy
     (invisible to congruences). Heegner (h=1: -3,-7,-11) = 0. (-15 = pure genus; -23 = pure deep.)
  P2 TOURNAMENTS | score sequence: H = log2(#iso classes sharing a score sequence). The TRANSITIVE
     tournament (score 0..n-1) is the UNIQUE realization => zero entropy = RIGID (the AP vertex); the
     regular score sequence has the largest fiber => max entropy (the kps reconstruction wall).
  P3 the DUAL entropies: score-DISTRIBUTION entropy (S217, transitive MAX spread) vs RECONSTRUCTION/fiber
     entropy (transitive MIN=0). The AP is max-spread AND zero-hidden-info: rigid.
  P4 REALS | continued-fraction prefix: golden [0;1,1,..] = min CF-entropy (worst-approx, the LRC FOIL
     S206) vs anti-golden t*=[0;13,14] = max (well-approx, the extremal).
  P5 NULLCONE | moment depth: detection depth = the depth where the entropy hits 0; LRC finite (bounded
     danger alphabet, Bonferroni depth ~5) vs GMC infinite (unbounded radial degree).
"""
from math import log2, gcd
from itertools import combinations, permutations

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

sep("STATUS: HISTORICAL / REFUTED SYNTHESIS -- READ MISTAKE-231")
print("  Retained: small reduced-form tables, n<=5 score-fiber census, and elementary arithmetic.")
print("  Retracted: one repo-wide entropy invariant, regular=max fiber, CF extrema, depth-5 classification, and LRC transfer.")

# ==========================================================================
sep("P1  BINARY FORMS | genus: h = genera x deep; the DEEP part is the hidden arithmetic entropy")
def is_reduced(a,b,c): return (-a<b<=a<=c) and not (a==c and b<0)
def reduced_forms(D):
    fs=[]; a=1
    while a*a<=-D//3+1:
        for b in range(-a,a+1):
            if (b*b-D)%(4*a)==0:
                c=(b*b-D)//(4*a)
                if a<=c and gcd(gcd(a,abs(b)),c)==1 and is_reduced(a,b,c): fs.append((a,b,c))
        a+=1
    return fs
def prime_factors(n):
    n=abs(n); s=set(); d=2
    while d*d<=n:
        while n%d==0: s.add(d); n//=d
        d+=1
    if n>1: s.add(n)
    return s
def num_genera(D):   # 2^{t-1}, t = # prime factors of D (fundamental-disc heuristic; exact for our examples)
    t=len(prime_factors(D))
    return 2**(t-1)
for D in (-3,-4,-7,-8,-11,-15,-20,-23,-24,-47,-84):
    h=len(reduced_forms(D)); g=num_genera(D)
    if h%g!=0: g=1  # guard: fall back if the heuristic mismatches
    deep=h//g
    Hgenus=log2(g); Hdeep=log2(deep)
    print(f"  disc {D:4d}: h={h:2d} = genera {g} x deep {deep:2d} ; H_genus(local)={Hgenus:.2f} + H_deep(HIDDEN)={Hdeep:.2f} bits")
print("  => Heegner (-3,-7,-11) h=1: zero entropy (rigid). -15=-3*5: pure GENUS (congruence-detectable, 0 deep).")
print("     -23,-47 (prime disc): pure DEEP -- the class group is invisible to any congruence = true hidden info.")

# ==========================================================================
sep("P2  TOURNAMENTS | score sequence: transitive = UNIQUE realization (0 entropy = rigid) vs regular (max)")
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
def scores(adj,n): return tuple(sorted(sum(adj[i]) for i in range(n)))
for n in (3,4,5):
    perms=list(permutations(range(n))); classes={}
    for adj in all_tournaments(n):
        c=canon(adj,n,perms)
        if c not in classes: classes[c]=adj
    fib={}
    for adj in classes.values():
        s=scores(adj,n); fib.setdefault(s,0); fib[s]+=1
    transitive_seq=tuple(range(n)); Htrans=log2(fib.get(transitive_seq,0)) if fib.get(transitive_seq,0) else 0.0
    maxfib=max(fib.values()); maxseq=[s for s in fib if fib[s]==maxfib]
    print(f"  n={n}: #classes={len(classes)}; transitive score {transitive_seq} realized by {fib.get(transitive_seq)} class (H={Htrans:.2f}); "
          f"max fiber={maxfib} at score(s) {maxseq[0]} (H={log2(maxfib):.2f})")
print("  => the TRANSITIVE tournament (score 0..n-1) is the UNIQUE realization: reconstruction entropy 0 = RIGID")
print("     (the AP/nullcone vertex). The largest fiber (near-regular scores) carries the hidden entropy (kps wall).")

# ==========================================================================
sep("P3  the DUAL: score-DISTRIBUTION entropy (transitive MAX) vs RECONSTRUCTION entropy (transitive 0)")
from collections import Counter
def dist_entropy(seq):
    n=len(seq); c=Counter(seq); return -sum((v/n)*log2(v/n) for v in c.values())
for n in (5,7,11):
    trans=list(range(n)); reg=[(n-1)//2]*n
    print(f"  n={n}: transitive -> score-DIST entropy {dist_entropy(trans):.2f} (MAX, spread) but RECONSTRUCTION entropy 0 (unique);")
    print(f"        regular    -> score-DIST entropy {dist_entropy(reg):.2f} (MIN, delta)  but RECONSTRUCTION entropy MAX (huge fiber).")
print("  => the two entropies are DUAL. The AP/transitive is MAX spread + ZERO hidden info = rigid; the")
print("     regular/Paley pole is MIN spread + MAX hidden info. Rigidity = max-order-yet-uniquely-determined.")

# ==========================================================================
sep("P4  REALS | continued-fraction prefix: golden (min = FOIL) vs anti-golden t*=14/183 (max = extremal)")
def cf(p,q,K=12):
    a=[]
    for _ in range(K):
        a.append(p//q); p,q=q,p-(p//q)*q
        if q==0: break
    return a
def geo_mean(a):
    prod=1.0
    for x in a: prod*=max(x,1)
    return prod**(1/len(a))
tstar=cf(14,183); golden=[1]*8  # golden ratio CF is all 1s
KHINCHIN=2.6854
print(f"  Khinchin's constant (a.e. geometric mean of partial quotients) ~= {KHINCHIN}")
print(f"  golden [0;1,1,..]: partial quotients {golden}, geo-mean {geo_mean(golden):.3f} << Khinchin = WORST-approx (min CF-entropy) = LRC FOIL (S206)")
print(f"  t* = 14/183 = [0;{','.join(map(str,tstar[1:]))}]: geo-mean {geo_mean(tstar[1:]):.3f} >> Khinchin = WELL-approx (max CF-entropy) = the extremal")
print("  => the CF 'arithmetic entropy' (info per convergent) is MINIMAL for golden (foil), MAXIMAL for t* (extremal).")

# ==========================================================================
sep("P5  NULLCONE | moment depth = certificate entropy: LRC finite (bounded alphabet) vs GMC infinite")
print("  LRC(14): danger count X in {0,..,13} is BOUNDED => the Bonferroni/moment certificate TERMINATES at")
print("           finite detection depth (~5, klein-S389/THM-671): FINITE arithmetic entropy (a finite certificate).")
print("  GMC/NC2: radial degree is UNBOUNDED => no finite moment depth determines the nullcone (Watson/Laplace")
print("           determinacy, S211): INFINITE arithmetic entropy (needs the full radial descent).")
print("  => detection depth = the moment depth at which local data pins the global object = a certificate entropy.")

sep("SUMMARY -- one information deficit across the repo")
print("""  H_arith(X | L) = log2 |{X': L(X')=L(X)}| = the GLOBAL bits HIDDEN from the LOCAL invariant. Instances:
    binary forms | genus     : DEEP class group (Heegner=0)          -- the S217 result, refined by genus theory
    tournaments  | score seq : transitive UNIQUE (0) / regular (max) -- the AP is the zero-entropy rigid vertex
    reals        | CF prefix : golden min (foil) / t* max (extremal) -- Diophantine info rate (S206)
    nullcone     | moment    : LRC finite depth / GMC infinite        -- certificate entropy (S211)
  UNIFYING: the repo's RIGID extremum -- the AP / transitive tournament / Heegner h=1 -- is the ZERO-arithmetic-
  entropy point (local determines global), and the DIFFICULTY of every thread is exactly the positive-entropy
  'hidden global' object (deep class group / cospectral fiber / CF tail / deep moment). Rigidity is why the
  extremal is unique; hidden entropy is where the proof must still go.""")
