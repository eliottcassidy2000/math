"""
Poke connections of the odd-n covering-min to repo structures:
 (1) the SAFE BAND = an independent set in the resonance/Cayley graph -> LRC = fractional-chromatic search.
 (2) the off-cusp odd core (descent) is a Z_p subset = a circulant; apex gap = its Cayley spectral gap.
 (3) circulant TOURNAMENTS (QR/Paley) on Z_p, their apex gaps, odd cycles (OCF) -- the tournament side.
"""
import cmath, math
import numpy as np
def apex_gap(O,p):
    w=cmath.exp(2j*cmath.pi/p)
    return min(abs(sum(w**((k*x)%p) for x in O))**2 for k in range(1,p))
p=7
print("(2)/(3) APEX GAPS of Z_7 subsets (odd cores / tournament connection sets) = Cayley spectral gaps:")
sets={"doublet {0,1}":[0,1],"QR/Paley {1,2,4}":[1,2,4],"covmin-7 odd core {0,1,5}":[0,1,5],
      "NQR {3,5,6}":[3,5,6],"{1,3,5}":[1,3,5],"full Z_7 (cusp)":[0,1,2,3,4,5,6]}
for name,O in sets.items():
    print(f"   {name:>26}: apex gap g(O) = {apex_gap(O,p):.5f}  |O|={len(O)}")
print("   => QR/Paley {1,2,4} = the Paley TOURNAMENT connection set; its apex gap is a tournament invariant.")
print("      the covmin-7 odd core {0,1,5} is OFF-cusp (proper, gap>0): a circulant on Z_7.")
print()
# (3) Paley tournament on Z_7: i->j iff j-i in QR; count its 3-cycles (the OCF odd cycles)
QR={1,2,4}
def paley_3cycles(p,QR):
    c=0
    for a in range(p):
        for b in range(p):
            for cc in range(p):
                if len({a,b,cc})<3: continue
                if ((b-a)%p in QR) and ((cc-b)%p in QR) and ((a-cc)%p in QR): c+=1
    return c//3
print(f"(3) Paley TOURNAMENT on Z_7 (QR={sorted(QR)}): directed 3-cycles = {paley_3cycles(p,QR)} (the OCF's odd cycles)")
print(f"    H(Paley_7) is odd (Redei); the apex-7 = the LRC C_7 = the Paley/circulant tournament's odd-cycle core.")
print()
# (1) the safe band as an independent set / the danger graph
print("(1) SMARTER SEARCH = fractional-chromatic on the resonance graph:")
print("    M(S)>=m/q  <=>  the speeds avoid the m-ball of 0 at witness a/q  <=> {sa mod q} is an INDEPENDENT")
print("    SET in the circulant 'danger graph' C(q; +-1..+-(m-1)) (vertices Z/q, edges within distance m).")
print("    => covering-min = min over covering S of (1 - fractional-clique / q)-type bound; search the")
print("       DANGER-GRAPH INDEPENDENT SETS that are simultaneously COVERING -- a constraint-satisfaction, not brute force.")
print()
print("SMARTER SEARCH IDEAS (poked):")
print("  A. WITNESS-FIRST/band: fix (q,a,m), allowed speed-residues = a^{-1}*{m..q-m}; ENUMERATE covering sets")
print("     with residues in the band (not greedy) -> densest = covmin. Band restricts residues mod q.")
print("  B. DESCENT/APEX: covmin off-cusp has a proper odd core O (circulant); minimize over O (apex gap g(O)>0)")
print("     then lift. The apex gap = Cayley spectral gap (THM-590) bounds the per-level density.")
print("  C. CIRCULANT-TOURNAMENT: the odd core as a tournament connection set; the OCF odd-cycle count is an")
print("     invariant; search connection sets with the right spectral/odd-cycle profile.")
print("  D. FRACTIONAL-CHROMATIC: covering-min = a circular-chromatic / independent-set problem on the danger")
print("     circulant; use LP relaxation + covering constraints (IP), far smarter than enumerating sets.")
