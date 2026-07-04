"""kps-2026-07-04-S5: exact witness t*=p/q + residue table for each drop-j one-swap family
F(j,X)={1..13}\{j} u {X}, to Lean-prove the whole one-swap covering stratum (klein-S127 HYP-4082).
For each j and covering ladder X=X0*k, find the max-min witness and its residue table."""
from math import gcd
from fractions import Fraction

def resdist(a, q):  # distance of a to nearest multiple of q, as integer in [0, q/2]
    r = a % q
    return min(r, q - r)

def M_and_witness(V, Qmax):
    """exact max over t=p/q (q<=Qmax) of min_i ||v_i t||; returns (M:Fraction, p, q, minres)."""
    best = (Fraction(0), 0, 1, 0)
    for q in range(2, Qmax + 1):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            md = min(resdist(v * p, q) for v in V)
            val = Fraction(md, q)
            if val > best[0]:
                best = (val, p, q, md)
    return best

def covering(V):
    return all(any(v % q == 0 for v in V) for q in range(2, 15))

# X0 per drop: smallest X making F(j,X) covering with X the "spread" runner (multiple of lcm(j,14)/...)
# Empirically the ladder base. Let's just scan X and find covering ones + the witness.
print("drop-j one-swap families: witness t*=p/q, M, residue table (V at t*)")
print("="*70)
ladders = {}
for j in range(1, 14):
    base = [v for v in range(1, 14) if v != j]
    # find X0 = smallest covering X (X>13), then ladder X0*k... but covering X may not be multiples.
    # Instead: for k=1,2,3 use X = (candidate). First find all covering X up to 300.
    covX = [X for X in range(14, 400) if covering(base + [X])]
    if not covX:
        print(f"j={j:2d}: NO covering X<=400 (base {base} already covers 14? {any(v%14==0 for v in base)})")
        # base has no mult of 14 unless j avoids... {1..13} has none; need X mult of 14
        continue
    X0 = covX[0]
    row = []
    for k in range(1, 4):
        X = X0 * k
        if not covering(base + [X]):
            # find the k-th covering multiple-ish; use covX index
            X = covX[k-1] if k-1 < len(covX) else X0*k
        V = base + [X]
        M, p, q, md = M_and_witness(V, min(400*k + 100, 2500))
        row.append((k, X, M, p, q, md))
    ladders[j] = (X0, row)
    k1 = row[0]
    print(f"j={j:2d}: X0={X0:3d} | k=1: X={k1[1]:3d} M={k1[2]}={float(k1[2]):.5f} t*={k1[3]}/{k1[4]}  (14/183={14/183:.5f})")

print("\n=== LADDER FITS (p,q,M vs k) + drop-12 (residue-liar) cross-check ===")
for j in [12, 13, 11, 9, 5, 1]:
    if j not in ladders: continue
    X0, row = ladders[j]
    print(f"\ndrop-{j} (X0={X0}):")
    for (k, X, M, p, q, md) in row:
        V = [v for v in range(1,14) if v!=j] + [X]
        residues = sorted([(v, (v*p) % q, resdist(v*p, q)) for v in V], key=lambda t: t[2])
        binding = [(v, r, d) for (v,r,d) in residues if d == md]
        print(f"  k={k}: X={X:4d} M={M} t*={p}/{q} md={md}  binding runners(d={md}): {[(v,r) for v,r,d in binding]}")
