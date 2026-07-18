# opus-2026-07-17-S375 -- HYP-7580: THE EXTENDED SIEVE LEMMA FOR q > 14.
#
# TARGET AS POSED: "q divides no speed => some p/q is lonely", for q > 14.
# My own S374 data already says this fails 34.8% of the time when q0 > 14, so
# before attempting a proof I verify the failure explicitly and extract the
# structural reason.
#
# THE FORBIDDEN WINDOW.  p/q is lonely iff for every speed v,
#     v*p mod q  NOT IN  W_q := { r : min(r, q-r)*14 < q }.
# W_q = {0, +-1, ..., +-floor((q-1)/14)}, so
#     |W_q| = 2*floor((q-1)/14) + 1.
# |W_q| = 1  <=>  q <= 14.  THAT is exactly why 14 is the classical boundary:
# for q <= 14 the only forbidden residue is 0, so "q divides no speed" IS the
# whole hypothesis.  For q > 14 the window widens and the hypothesis is no
# longer sufficient.
from math import gcd
from functools import reduce

def W(q): return set([0] + [r%q for j in range(1,(q-1)//14+1) for r in (j,-j)])
def lonely_at(V,p,q):
    Wq=W(q)
    return all((v*p)%q not in Wq for v in V)
def works(V,q): return any(lonely_at(V,p,q) for p in range(1,q) if gcd(p,q)==1)

print("(1) THE WINDOW |W_q| = 2*floor((q-1)/14)+1 -- why 14 is the boundary")
for q in [7,13,14,15,16,28,29,30,42,43]:
    print(f"    q={q:3d}: |W_q| = {len(W(q)):3d}   W_q = {sorted(W(q))[:7]}{'...' if len(W(q))>7 else ''}")
print("    |W_q| = 1 exactly for q <= 14  =>  classical hypothesis 'q divides no")
print("    speed' is EXACTLY 'no speed lands in W_q', and is sufficient only there.")

print()
print("(2) EXPLICIT COUNTEREXAMPLE to the extended lemma (from the S374 run)")
V=[11,70,77,137,144,156,175,213,226,232,246,262,281]
q=15
print(f"    V = {V}")
print(f"    q = {q}; does q divide any speed?  {[v for v in V if v%q==0]}  (empty => hypothesis HOLDS)")
print(f"    gcd(V) = {reduce(gcd,V)} (primitive)")
print(f"    W_15 = {sorted(W(15))}")
print("    p : residues v*p mod 15, and which speed lands in W_15")
for p in range(1,15):
    if gcd(p,15)!=1: continue
    hit=[(v,(v*p)%15) for v in V if (v*p)%15 in W(15)]
    print(f"      p={p:2d}: blocked by v={hit[0][0]} (v*p = {hit[0][1]} mod 15)" if hit else f"      p={p:2d}: LONELY")
print(f"    => no p works at q=15, yet 15 divides no speed.  THE LEMMA IS FALSE.")

print()
print("(3) WHY NO COUNTING PROOF CAN EXIST AT A SINGLE q")
print("    each speed v kills at most |W_q|-1 numerators (v*p = w has one")
print("    solution p per nonzero w), so  #bad p <= 13*(|W_q|-1) = 26*floor((q-1)/14).")
print("    q    phi(q)   13*(|W_q|-1)   does the bound fire?")
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
for q in [15,16,20,23,25,29,31,37,41,43,53,71]:
    bad=13*(len(W(q))-1)
    print(f"    {q:4d}  {phi(q):6d}   {bad:12d}   {'YES' if phi(q)>bad else 'no'}")
print("    26*floor((q-1)/14) ~ 1.857*q > q > phi(q), so the union bound NEVER")
print("    fires -- the same 13/7 > 1 obstruction as everywhere else in this program.")
