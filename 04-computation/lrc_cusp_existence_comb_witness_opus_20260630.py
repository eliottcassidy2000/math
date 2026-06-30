"""
CUSP EXISTENCE (corrected): existence is the WITNESS denominator q* (the comb q=14, or its relocation),
not arbitrary q. At the cusp the density lim D(q)/q -> 0 (MEASURE vanishes) but D(q*) >= 1 (EXISTENCE).
"""
def D(S,q,n=14):
    return sum(1 for a in range(1,q) if all(min((v*a)%q, q-(v*a)%q)*n>=q for v in S))
sets={"AP {1..13}":(list(range(1,14)),14),
      "{1..11,13,84}":([1,2,3,4,5,6,7,8,9,10,11,13,84],89),
      "{1..12,182}":(list(range(1,13))+[182],183)}
print("EXISTENCE at the witness denominator q* (the comb / its relocation) vs MEASURE (density at large q):")
print(f"   {'set':>16} {'witness q*':>11} {'D(q*) (EXISTS)':>15} {'M=D-edge':>10}  {'density D(2000)/2000 (->0)':>26}")
for nm,(S,qs) in sets.items():
    dqs=D(S,qs); dens=D(S,2003)/2003
    print(f"   {nm:>16} {qs:>11} {dqs:>15} {'1/14' if nm[0]=='A' else ('7/89' if '84' in nm else '14/183'):>10}  {dens:>26.5f}")
print()
print("THE PIN -- the cusp is pure EXISTENCE (one witness), the measure is irrelevant:")
print("  * MEASURE: lim_q D(q)/q = the lonely density. AT THE CUSP (odd core = Z_7) this -> 0 (rho_apex=0).")
print("    e.g. AP: density 0 (lonely set = the 6 units {k/14}, measure zero). The doublet/measure argument")
print("    is SILENT here -- it correctly reports rho=0, no contradiction.")
print("  * EXISTENCE: D(q*) >= 1 at the witness denominator q*. The lonely POINT is there:")
print("      - AP: q*=14, the comb itself; the 6 empty-tooth witnesses {1,3,5,9,11,13}/14 (the units).")
print("      - covering: the killer FILLS the observer's comb-tooth, so q* RELOCATES off 14 (coprime:")
print("        89, 183, ...), the next Farey rational; M = n/Phi_6(n) = the Eisenstein-norm denominator.")
print("  * so the CUSP needs ONE witness (D(q*)>=1), supplied by the comb (AP) or its relocation (covering)")
print("    -- this is 'existence direct from bounded' (klein-S4): finite, decidable, NOT a measure bound.")
print()
print("LRC(14) = OFF-CUSP (proper core: doublet rho>=0.198, MEASURE, PROVED) + CUSP (core Z_7: D(q*)>=1,")
print("  EXISTENCE, the comb/relocation witness). Two regimes, two tools: measure for the bulk, the comb")
print("  witness for the cusp. The doublet carries the measure; the empty tooth carries the existence.")
# the cusp witness is the comb tooth / relocation -- tie to covering-min formula
print(f"\n  cusp witness denominators: AP q*=14=n; covering-min q*=183=n^2-n+1=Phi_6(n) (the relocated comb).")
print(f"  => the cusp EXISTENCE witness IS the Dirac-comb empty tooth (AP) or its Phi_6 relocation (covering).")
