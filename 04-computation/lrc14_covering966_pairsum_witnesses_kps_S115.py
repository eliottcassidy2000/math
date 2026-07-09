from itertools import combinations
from math import gcd

def is_cov(S): return all(any(s % q == 0 for s in S) for q in range(2, 15))
def prim(S):
    g = 0
    for s in S: g = gcd(g, s)
    return g == 1

def pairsum_moduli(v):
    return sorted({v[i]+v[j] for i in range(len(v)) for j in range(i,len(v))})

def band_witness(v):
    """(p,q): q a pair-sum modulus, all residues (v_l*p mod q) in [q/14,13q/14] (i.e. q<=14r<=13q)."""
    for q in pairsum_moduli(v):
        if q < 2: continue
        for p in range(1, q):
            if all(q <= 14*((vl*p)%q) <= 13*q for vl in v):
                return (p, q)
    return None

# enumerate the 966
covering = [S for S in combinations(range(1,19),13) if prim(S) and is_cov(S)]
print(f"#covering primitive [1,18] 13-sets = {len(covering)}")

# witnesses
wit = {}
missing = []
qmax = 0; pmax = 0
for S in covering:
    w = band_witness(list(S))
    if w is None:
        missing.append(S)
    else:
        wit[S] = w
        qmax = max(qmax, w[1]); pmax = max(pmax, w[0])
print(f"#with pair-sum band witness = {len(wit)} / {len(covering)}")
print(f"#missing = {len(missing)}")
if missing:
    print("MISSING (no pair-sum band witness):", missing[:5])
print(f"max q used = {qmax}, max p used = {pmax}")
# sanity: verify each witness truly clears (min residue >= q/14)
bad=0
for S,(p,q) in wit.items():
    if not all(q <= 14*((s*p)%q) <= 13*q for s in S): bad+=1
print(f"witness re-verification failures = {bad}")
# save witnesses for the Lean generator
import json
with open("/tmp/wit966.json","w") as f:
    json.dump([[list(S), wit[S][0], wit[S][1]] for S in covering if S in wit], f)
print("saved /tmp/wit966.json")
