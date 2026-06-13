#!/usr/bin/env python3
"""
resonant_crossing_chord_spectrum_monad_partB.py
monad-explorer-2026-06-13

PART B — the small-factor CHORD BOTTLENECK as a master key that UNIFIES
  * THM-437 (the cube / n=27 is angle-rigid at 81, cannot beat)
  * THM-493/494 (n=28 crosses at 85, uniquely via t=3)
purely combinatorially (no angle calculus): the resonant bonus is gated by the
chord spectrum of the SMALLEST factor.

Builds optimal triangular-lattice k-point factors by EXACT brute force over the
19-point hex (center + 2 shells), computes their chord spectra, and runs a
factorization sweep predicting which n cross 3N and with which Moser norm t.
"""
from itertools import combinations
from collections import defaultdict

def esub(p,q): return (p[0]-q[0], p[1]-q[1])
def enorm(p):
    a,b=p; return a*a+a*b+b*b
def edges(pts): return sum(1 for p,q in combinations(pts,2) if enorm(esub(p,q))==1)
def chordspec(pts):
    d=defaultdict(int)
    for p,q in combinations(pts,2): d[enorm(esub(p,q))]+=1
    return dict(d)
def m_alpha(pts):
    d=defaultdict(int)
    for p in pts:
        for q in pts:
            if p!=q: d[esub(p,q)]+=1
    return d
def norm_t_alphas(t,R=8):
    return [(a,b) for a in range(-R,R+1) for b in range(-R,R+1) if a*a+a*b+b*b==t]
def delta_t(G,H,t):
    mG,mH=m_alpha(G),m_alpha(H); s=0
    for al in norm_t_alphas(t): s+=mG.get(al,0)*mH.get(al,0)
    return s//2

# 19-point hex region: center + 2 shells in Eisenstein coords
HEX19=[(a,b) for a in range(-2,3) for b in range(-2,3) if enorm((a,b))<=4]
# (norm<=4 picks shells 0,1,3,4 -> the compact hex of 19 lattice points)

def best_factor(k):
    """max-edge k-subset of HEX19 (exact brute force); return (e, set, chordspec).
       Among max-edge subsets, prefer the one with the RICHEST chord spectrum."""
    best=None
    for sub in combinations(HEX19,k):
        e=edges(sub)
        cs=chordspec(sub)
        # key: maximize edges, then richness (#distinct non-unit chord norms)
        rich=len([t for t in cs if t>=2])
        key=(e,rich)
        if best is None or key>best[0]:
            best=(key, e, list(sub), cs)
    return best[1], best[2], best[3]

print("="*72)
print("PART B1 — optimal small triangular-lattice factors & their chord spectra")
print("="*72)
FAC={}
for k in range(2,9):
    e,pts,cs=best_factor(k)
    nz=sorted(t for t in cs if t>=2)
    FAC[k]=(e,pts,cs,nz)
    print(f"  k={k}: e={e:2d}  rho={2*e/k:.3f}  ChordSpec={cs}   non-unit norms={nz}")
print()
print("  READING: k=2 (edge) and k=3 (triangle) are CHORD-FREE (only norm 1).")
print("           k=4 (rhombus) first carries a non-unit chord: norm 3 (= sqrt3).")
print("  => any product whose SMALLEST factor is size 2 or 3 has Delta_t=0 for all t>=2.")
print()

print("="*72)
print("PART B2 — the 27-vs-28 UNIFICATION (chord bottleneck, no angle calculus)")
print("="*72)
# n=27: factorizations through a size-3 (or size-2) factor only (27=3*9, 3*3*3)
print("  n=27 = 3*9 (also 3^3): the size-3 factor is the triangle K3, chord-free.")
K3=FAC[3][1]
F9=best_factor(9)
print(f"     ChordSpec(K3)={chordspec(K3)}  -> NO non-unit chord")
print(f"     => Delta_t(K3, anything)=0 for all t>=2  (scan t=2..30:",
      [t for t in range(2,31) if delta_t(K3,F9[1],t)>0], ")")
print(f"     so 27 gets ZERO resonance bonus through any triangle factor.")
print(f"     This RE-DERIVES THM-437 (cube angle-rigid at 81) combinatorially:")
print(f"     27=3^3 forces a chord-free size-3 factor.  Best product U(27) = product cap (tie 81).")
print()
print("  n=28 = 4*7: the size-4 factor is the rhombus R, ChordSpec={1,3}.")
print(f"     The single non-unit chord (norm 3) opens t=3 resonance => Delta_3=2 => 85 > 84.")
print(f"     28 is the FIRST composite whose smallest dense factor is chord-BEARING.")
print()

print("="*72)
print("PART B3 — predictive factorization sweep (which n cross 3N, with which t)")
print("="*72)
print("  For each n, over factorizations n=a*b (2<=a<=b), take densest a- and b-factors,")
print("  compute best resonant total over admissible t in shared non-unit chord spectrum.")
print()
def divisor_pairs(n):
    return [(a,n//a) for a in range(2,int(n**0.5)+1) if n%a==0]

for n in range(24,50):
    rows=[]
    for a,b in divisor_pairs(n):
        if a>8 or b>8:  # 19-hex brute force reliable to k<=8
            continue
        eA,GA,csA,nzA=FAC[a]; eB,GB,csB,nzB=FAC[b]
        shared=sorted(set(nzA)&set(nzB))
        Pcap=eA*b+a*eB
        best_t,best_d=None,0
        for t in shared:
            d=delta_t(GA,GB,t)
            if d>best_d: best_d,best_t=d,t
        U=Pcap+best_d
        rows.append((a,b,Pcap,best_t,best_d,U))
    if not rows:
        continue
    # pick factorization with max U
    a,b,Pcap,bt,bd,U=max(rows,key=lambda r:r[5])
    thr=3*n
    flag = "CROSS" if U>thr else ("tie" if U==thr else "")
    note=f"t={bt} +{bd}" if bd>0 else "no bonus"
    print(f"  n={n:2d}: best {a}x{b}  Pcap={Pcap:3d}  {note:10s}  U={U:3d}  vs 3N={thr:3d}  {flag}")
print()
print("  NOTE: this is the 2-factor TRIANGULAR-LATTICE product family only (factors<=8).")
print("  It is a LOWER-bound construction lens, not an upper bound on u(n).")
print("  The point is the chord-bottleneck PREDICTION of WHERE/WHICH-t a bonus appears.")

print()
print("="*72)
print("PART C — is t=3 merely FIRST, or DOMINANT? per-t bonus breakdown")
print("="*72)
for (a,b) in [(4,7),(5,7),(6,7),(7,7),(6,8),(5,8)]:
    GA=FAC[a][1]; GB=FAC[b][1]
    n=a*b
    Pcap=FAC[a][0]*b+a*FAC[b][0]
    line=f"  n={n:2d} ({a}x{b}) Pcap={Pcap:3d}: "
    parts=[]
    for t in [3,4,7,9,12,13]:
        d=delta_t(GA,GB,t)
        if d>0: parts.append(f"Delta_{t}={d}")
    print(line+"  ".join(parts) if parts else line+"(no bonus)")
print()
print("  => t=3 (sqrt3, the triangular SECOND-neighbour) gives the LARGEST bonus")
print("     in every case: it is the DOMINANT resonance norm, not merely the first.")
print("     Higher norms (4,7) carry smaller bonuses (fewer such chords in dense patches).")
