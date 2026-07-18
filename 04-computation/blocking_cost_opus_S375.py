# opus-2026-07-17-S375 -- THE BLOCKING COST, the sharp replacement for the lemma.
#
# A speed v coprime to q kills the numerators p = v^{-1} w for w in W_q\{0}.
# Only w coprime to q yields a p coprime to q, so v kills exactly
#     k_q := #{ w in W_q, w != 0, gcd(w,q) = 1 }
# of the phi(q) available numerators.  Hence, WITHOUT using a multiple of q,
#     blocking q costs at least  ceil(phi(q) / k_q)  speeds.
# Blocking VIA a multiple (q | v) costs ONE speed and kills every p at once --
# that is the cheap route, and the lcm construction of THM-1105 is exactly it.
# So the adversary faces a real budget: 13 speeds, and each q > 14 must be paid
# for either by one multiple or by ceil(phi(q)/k_q) residues.
from math import gcd
def W(q): return set([0]+[r%q for j in range(1,(q-1)//14+1) for r in (j,-j)])
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def k_q(q): return sum(1 for w in W(q) if w!=0 and gcd(w,q)==1)

print("(4) BLOCKING COST ceil(phi(q)/k_q) -- speeds needed to block q by RESIDUES")
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
print("(5) THE EXPENSIVE MODULI -- these force the adversary to spend a multiple")
exp=[(q,c) for q,c in rows if c is None or c>13]
print(f"    moduli in [15,60] with blocking cost > 13 (or infinite): {[q for q,c in exp]}")
print(f"    count: {len(exp)} of {len(rows)}")
print()
print("    Any family blocking ALL of these must contain, for each such q,")
print("    a speed divisible by q.  With only 13 speeds that is a strong")
print("    simultaneous-divisibility constraint -- and it is exactly the")
print("    structure the THM-1105 lcm construction exploits.")

print()
print("(6) SANITY: does the cost prediction match the S375 counterexample?")
V=[11,70,77,137,144,156,175,213,226,232,246,262,281]
q=15
Wq=W(q); killers=[v for v in V if any((v*p)%q in Wq for p in range(1,q) if gcd(p,q)==1)]
print(f"    q=15: predicted cost ceil(phi/k) = {-(-phi(15)//k_q(15))} speeds")
print(f"    speeds actually doing the blocking: {sorted(set(killers))[:6]}")
print(f"    number of distinct blockers: {len(set(killers))}  (prediction was a LOWER bound)")
