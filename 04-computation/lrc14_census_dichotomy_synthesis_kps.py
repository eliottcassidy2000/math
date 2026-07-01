#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CENSUS SYNTHESIS: is inf meas(L_C) over 11-cores attained at a FULL-ORBIT (Z/N)* core, and is the
near-tight locus a DICHOTOMY {full unit orbit, partial two-clash}?  Synthesizing HYP-3793's modular
labels with the loop-map group structure (lrc14_loop_map_dictionary_kps.py).

kind-pasteur-2026-07-01-S8.  For a broad pool of 11-cores compute (M, primitive N of the maximizers,
#atoms, full-(Z/N)*-orbit?, meas(L_C@1/14)); report the global min, the bottom locus by family, and the
min measure per modular class N.  Test: (i) min >= 1/36, (ii) argmin is the full-orbit (Z/10)* core,
(iii) every near-tight core is full-orbit OR partial two-clash (the 11-speed THM-523 echo).
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import gcd
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
def units(N): return [a for a in range(1,N) if gcd(a,N)==1]
def M_and_atoms(S,Dmax=None):
    # M(S)=max_{a/b} min_v ||v a/b||; integer hot loop: d_v=min(va%b, b-va%b), M-candidate = (min_v d_v)/b.
    if Dmax is None: Dmax=2*max(S)+2
    dens=set()
    for v in S: dens.add(2*v)
    for v in S:
        for w in S:
            dens.add(v+w)
            if v!=w: dens.add(abs(v-w))
    bestnum=0; bestden=1; atoms=[]   # best value = bestnum/bestden
    for b in sorted(d for d in dens if 0<d<=Dmax):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            md=b
            for v in S:
                r=(v*a)%b; dv=r if r<=b-r else b-r
                if dv<md: md=dv
                if md==0: break
            if md==0: continue
            # compare md/b vs bestnum/bestden
            c=md*bestden - bestnum*b
            if c>0: bestnum,bestden=md,b; atoms=[Fr(a,b)]
            elif c==0 and md*bestden==bestnum*b: atoms.append(Fr(a,b))
    return Fr(bestnum,bestden), sorted(set(atoms))
BAND=Fr(1,14)
def safe_arcs(v): return [((Fr(k)+BAND)/v,(Fr(k+1)-BAND)/v) for k in range(v)]
def inter(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: r.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return r
def Lmeas(S):
    a=safe_arcs(S[0])
    for v in S[1:]:
        a=inter(a,safe_arcs(v))
        if not a: return Fr(0)
    return sum(h-l for l,h in a)
def classify(S):
    M,atoms=M_and_atoms(S)
    if M<=BAND: return None  # tight or covering: meas 0, not a loose core
    Ns=sorted(set(t.denominator for t in atoms))
    N=Ns[0] if len(Ns)==1 else None
    full = (N is not None and sorted(t.numerator for t in atoms)==units(N))
    return dict(S=S, M=M, N=N, na=len(atoms), full=full, meas=Lmeas(S))

pool=set()
# structured 2-drops of {1..V}
for V in [13,14,15,16]:
    base=set(range(1,V+1))
    for drop in itertools.combinations(range(1,V+1), V-11):
        C=tuple(sorted(base-set(drop)))
        if len(C)==11: pool.add(C)
# 1-drop of {1..12} + one far w
for j in range(1,13):
    for w in range(13,61):
        C=tuple(sorted(set(range(1,13))-{j}|{w}))
        if len(C)==11: pool.add(C)
# random 11-cores in [1,22]
rng=random.Random(1)
for _ in range(2500):
    pool.add(tuple(sorted(rng.sample(range(1,23),11))))
rows=[]
for C in pool:
    r=classify(list(C))
    if r is not None and r['meas']>0: rows.append(r)
rows.sort(key=lambda r: r['meas'])
print(f"POOL: {len(pool)} candidate 11-cores; {len(rows)} are LOOSE (M>1/14, meas>0).")
print(f"1/36 = {float(Fr(1,36)):.5f};  1/28 = {float(Fr(1,28)):.5f}")
print("\n=== BOTTOM 14 by meas(L_C) (the near-tight locus) ===")
print(f"  {'meas':>10} {'M':>7} {'N':>5} {'#at':>4} {'family':>10}  core")
for r in rows[:14]:
    fam = 'FULL(Z/N)*' if r['full'] else ('2-clash' if r['N'] else 'multi-N')
    print(f"  {float(r['meas']):>10.5f} {str(r['M']):>7} {str(r['N']):>5} {r['na']:>4} {fam:>10}  {r['S']}")
gmin=rows[0]
print(f"\nGLOBAL MIN over pool: meas={gmin['meas']}={float(gmin['meas']):.5f} at {gmin['S']}")
print(f"  M={gmin['M']}, atoms on (Z/{gmin['N']})*, full-orbit={gmin['full']}; >= 1/36? {gmin['meas']>=Fr(1,36)}")

print("\n=== MIN meas per MODULAR CLASS N (full-orbit cores only) ===")
byN={}
for r in rows:
    if r['full']:
        byN.setdefault(r['N'], []).append(r['meas'])
for N in sorted(byN):
    mn=min(byN[N])
    print(f"  N={N:>3} ((Z/{N})*, phi={len(units(N))}): min meas={float(mn):.5f}={mn}  #cores={len(byN[N])}  >=1/36? {mn>=Fr(1,36)}")

print("\n=== DICHOTOMY TEST on the bottom 30: full-orbit OR 2-clash? ===")
fam_counts={'FULL':0,'2CLASH':0,'OTHER':0}
for r in rows[:30]:
    if r['full']: fam_counts['FULL']+=1
    elif r['N'] is not None and r['na']<=2: fam_counts['2CLASH']+=1
    else: fam_counts['OTHER']+=1
print(f"  bottom-30 family split: {fam_counts}")
allmin=min(r['meas'] for r in rows)
print(f"\nVERDICT: global min meas = {float(allmin):.5f} (>=1/36={allmin>=Fr(1,36)}); argmin full-orbit (Z/10)*={gmin['N']==10 and gmin['full']}.")
print("  Near-tight locus = {full unit orbit, partial 2-clash} (+ any OTHER flagged above) -- the 11-speed THM-523 echo.")
print("DONE.")
