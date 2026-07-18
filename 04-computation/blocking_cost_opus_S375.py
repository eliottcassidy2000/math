# opus-2026-07-17-S375 -- HISTORICAL UNIT-SPEED BLOCKING-COST TABLE.
#
# CORRECTION (codex-S67): every assertion below involving k_q is valid only
# when every speed is coprime to q.  A nonunit speed belongs to a larger gcd
# stratum and can kill more reduced numerators.  The general exact count is
# K(q,g)=phi(q)/phi(q/g) * #{r in W_q : gcd(r,q)=g}; see THM-1110 and
# extended_sieve_gcd_stratified_referee_codex_S67.py.  At q=90 the unrestricted
# reduced-numerator blocking cost is 3, while the primitive cost is exactly 4.
#
# A speed v coprime to q kills the numerators p = v^{-1} w for w in W_q\{0}.
# Only w coprime to q yields a p coprime to q, so v kills exactly
#     k_q := #{ w in W_q, w != 0, gcd(w,q) = 1 }
# of the phi(q) available numerators.  Hence, AMONG UNIT SPEEDS,
#     blocking q costs at least  ceil(phi(q) / k_q)  speeds.
# Blocking VIA a multiple (q | v) costs ONE speed and kills every p at once --
# that is the cheap route, and the lcm construction of THM-1105 is exactly it.
# This script does not price arbitrary nonzero residues.
from math import gcd
def W(q): return set([0]+[r%q for j in range(1,(q-1)//14+1) for r in (j,-j)])
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def k_q(q): return sum(1 for w in W(q) if w!=0 and gcd(w,q)==1)

print("CORRECTION: this is the UNIT-SPEED cost only; see THM-1110 for gcd strata.")
print("(4) UNIT-SPEED COST ceil(phi(q)/k_q)")
print("    q   phi(q)  |W_q|  k_q   cost   affordable with 13 speeds?")
rows=[]
for q in range(15,61):
    k=k_q(q)
    if k==0:
        print(f"    {q:3d} {phi(q):6d} {len(W(q)):6d} {k:4d}    INF    NO -- unblockable by residues!")
        rows.append((q,None)); continue
    c=-(-phi(q)//k)
    rows.append((q,c))
    print(f"    {q:3d} {phi(q):6d} {len(W(q)):6d} {k:4d} {c:6d}   {'yes' if c<=13 else 'NO -- must use a multiple'}")

print()
print("(5) UNIT-SPEED TABLE ONLY -- no conclusion for arbitrary nonunit speeds")
exp=[(q,c) for q,c in rows if c is None or c>13]
print(f"    moduli in [15,60] with blocking cost > 13 (or infinite): {[q for q,c in exp]}")
print(f"    count: {len(exp)} of {len(rows)}")
print()
print("    These rows do not force divisibility: a proper nonunit gcd stratum")
print("    may be cheaper than the displayed unit-speed cost.")

print()
print("(6) HISTORICAL q=15 DIAGNOSTIC (not a unit-speed lower-bound test)")
V=[11,70,77,137,144,156,175,213,226,232,246,262,281]
q=15
Wq=W(q); killers=[v for v in V if any((v*p)%q in Wq for p in range(1,q) if gcd(p,q)==1)]
print(f"    q=15: unit-speed cost ceil(phi/k) = {-(-phi(15)//k_q(15))} speeds")
print(f"    speeds actually doing the blocking: {sorted(set(killers))[:6]}")
print(f"    number of distinct blockers: {len(set(killers))}")
