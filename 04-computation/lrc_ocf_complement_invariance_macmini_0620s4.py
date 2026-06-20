#!/usr/bin/env python3
"""
lrc_ocf_complement_invariance_macmini_0620s4.py  (mac-mini-2026-06-20-S4)

Final consolidation: the OCF / odd-cycle structure (Claim A's object) is COMPLEMENT-INVARIANT,
so the central computation halves on the half-tiling.  Also cross-check half(n)=floor((n-1)^2/4)
by a SECOND independent route (Burnside orbit count of the tournament complement on tilings).
"""
import itertools, sys
sys.stdout.reconfigure(line_buffering=True)
def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in itertools.product([0,1],repeat=len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for (b,(i,j)) in zip(bits,pairs):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        yield adj
def comp(adj,n): return [[adj[j][i] for j in range(n)] for i in range(n)]
def odd_cycle_counts(adj,n):
    """vector of (#directed cycles of length L) for L=3,5,7 (odd)."""
    out={}
    for L in (3,5,7):
        if L>n: break
        c=0
        for verts in itertools.permutations(range(n),L):
            if verts[0]!=min(verts): continue   # canonical start
            # count directed cycles on this ordered tuple as a cycle
            if all(adj[verts[t]][verts[(t+1)%L]]==1 for t in range(L)): c+=1
        out[L]=c
    return out
print("OCF / odd-cycle complement-invariance (reverse all arcs):")
for n in [4,5,6]:
    tot=0; inv=0
    for adj in all_tournaments(n):
        tot+=1
        if odd_cycle_counts(adj,n)==odd_cycle_counts(comp(adj,n),n): inv+=1
    print(f"  n={n} ({tot} tournaments): odd-cycle-count vector invariant in {inv}/{tot}  {'ALL' if inv==tot else 'SOME FAIL'}")
print("  => OCF (all odd-cycle counts) complement-invariant => Claim-A / mu computation halves on the half-tiling.")

# independent half(n): Burnside orbits of T->T^op on the 2^m tilings (fixed base path) via tile-reflection
def tiles(n): return [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x-y>=2]
def sigma(t,n): x,y=t; return (n+1-y,n+1-x)
print("\nIndependent half(n) via Burnside (#orbits of sigma on 2^m tilings = (2^m + 2^fix)/2... NO:")
print("  half = #orbits of sigma on the m TILES = (m+#fixed-tiles)/2):")
for n in range(2,13):
    T=tiles(n); m=len(T); fix=sum(1 for t in T if sigma(t,n)==t)
    print(f"  n={n}: (m+fix)/2 = ({m}+{fix})/2 = {(m+fix)//2}  vs floor((n-1)^2/4)={(n-1)**2//4}  {'OK' if (m+fix)//2==(n-1)**2//4 else 'FAIL'}")
