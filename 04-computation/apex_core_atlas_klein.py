#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""APEX CORE ATLAS: a fine feature-map of the 127 Z_7 cores (klein-S21).

For each core O subset Z_7: gap g(O)=min_{k!=0}|Ohat(k)|^2 (THM-590, = lambda_min of the autocorr
circulant = lambda_min(|O|.I_diag + Cayley adjacency)), the WORST mode k*, the difference multiset,
flat? (all nonzero modes equal = a perfect difference set), the Z_7^* multiplier-orbit rep. Goal: the
EXACT rule assigning each core to its gap class, and the Cayley-graph identity of each class.
"""
import math, cmath, itertools
from collections import defaultdict
P=7; W=cmath.exp(2j*math.pi/P)
QR={1,2,4}  # quadratic residues mod 7

def spectrum(O):
    return [abs(sum(W**((k*x)%P) for x in O))**2 for k in range(P)]
def gap_and_mode(O):
    sp=spectrum(O); nz=[(round(sp[k],6),k) for k in range(1,P)]
    g=min(v for v,_ in nz); ks=[k for v,k in nz if abs(v-g)<1e-6]
    return g, ks, [round(x,4) for x in sp]
def diffset_flat(O):
    sp=spectrum(O); nz=[round(sp[k],4) for k in range(1,P)]
    return len(set(nz))==1   # flat nonzero spectrum <=> perfect difference set
def mult_orbit_rep(O):
    # Z_7^* multiplier orbit (u*O) + translates: canonical rep = lexicographically min over u in 1..6
    best=None
    for u in range(1,P):
        uO=tuple(sorted((u*x)%P for x in O))
        if best is None or uO<best: best=uO
    return best

print("="*86)
print(" APEX CORE ATLAS (Z_7): all 127 nonempty cores, classified by gap")
print("="*86)
byclass=defaultdict(list)
for r in range(1,P+1):
    for O in itertools.combinations(range(P), r):
        g,ks,sp = gap_and_mode(set(O))
        byclass[round(g,6)].append((O,ks,diffset_flat(set(O))))

doublet=4*math.cos(3*math.pi/7)**2
names={0.0:"ZERO (cusp)",round(doublet,6):"DOUBLET min",0.307979:"mid",1.0:"unit",2.0:"FLAT/difference-set"}
for g in sorted(byclass):
    items=byclass[g]
    sizes=defaultdict(int); flats=0; modeset=set()
    for O,ks,flat in items:
        sizes[len(O)]+=1; flats+=1 if flat else 0
        for k in ks: modeset.add(k)
    nm=names.get(round(g,6),"")
    print(f"\n gap = {g:.6f}  [{nm}]   ({len(items)} cores)")
    print(f"   sizes: {dict(sorted(sizes.items()))} ;  flat(difference-set): {flats}/{len(items)} ;  worst modes k*: {sorted(modeset)}")

print("\n"+"="*86)
print(" THE RULES (fine patterns):")
print("="*86)
# Rule: which 3-cores are flat (gap 2 = difference set)?
threes=[set(O) for O in itertools.combinations(range(P),3)]
flat3=[sorted(O) for O in threes if diffset_flat(O)]
print(f"\n [R1] 3-cores: {len(flat3)} are FLAT (gap 2 = perfect (7,3,1) difference sets); 35-{len(flat3)}={35-len(flat3)} give gap 0.308.")
print(f"      the flat 3-cores = QR {{1,2,4}} & QNR {{3,5,6}} and their 7 translates each:")
reps=set(mult_orbit_rep(O) for O in flat3)
print(f"      multiplier-orbit reps among flat 3-cores: {sorted(reps)}  (QR/QNR = Paley/Fano difference sets)")
# Rule: doublet worst mode
print(f"\n [R2] DOUBLET {{0,1}}: spectrum {gap_and_mode({0,1})[2]} ; worst mode k*={gap_and_mode({0,1})[1]} (the MIDDLE freq, |1+zeta^3|^2=2+2cos6pi/7).")
print(f"      ALL doublets {{a,a+d}} give the SAME gap 4cos^2(3pi/7) (any difference d!=0): the gap is")
for d in range(1,4):
    print(f"        difference d={d}: gap={gap_and_mode({0,d})[0]:.6f} worst k*={gap_and_mode({0,d})[1]}")
# Rule: singleton
print(f"\n [R3] SINGLETON {{0}}: spectrum all 1 (flat!) gap=1 ; co-singleton (6-core) also gap 1. The 'unit' class.")
# Rule: spread interpretation
print(f"\n [R4] SPREAD: gap measures FLATNESS of the Fourier spectrum.")
print(f"      FLATTEST (gap 2, all modes=2) = the difference sets (QR/Paley) = most 'spread'/random-like.")
print(f"      MOST CONCENTRATED (gap 0.198, two-term) = the DOUBLET = the BINDING obstruction.")
print(f"      => the apex obstruction is the LEAST-spread core (doublet), the OPPOSITE of the random/QR core.")
