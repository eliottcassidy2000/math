# Seven-vertex census: base span plus quadratic-twisted polynomial parts q(q-1) (D/q) q^j (j <= 2), D in {2, 3, -3, 5,
# 7, -7, 6, -6}, alone and together with the weight-3/5/6 candidate orbits.  mac-mini-2026-10-25-pyth.
import sys
sys.argv = [sys.argv[0]]
exec(open('tt7_candidate_fit.py').read().split("print(f\"{len(ps)} values")[0])
from sympy import jacobi_symbol
def kron(D, p): return int(jacobi_symbol(D % p, p)) if D % p else 0
for D in (2, 3, -3, 5, 7, -7, 6, -6):
    cols = dict(base)
    for j in range(3): cols[f'({D}/q)q^{j}'] = {p: p * (p - 1) * kron(D, p) * p ** j for p in ps}
    sol = in_span(cols)
    print(f"twist ({D}/q):", 'EXACT ' + str({k: str(v) for k, v in sol.items() if v}) if sol else 'NOT in span')
cols = dict(base)
for D in (2, 3, -3, 7, -7):
    for j in range(3): cols[f'({D}/q)q^{j}'] = {p: p * (p - 1) * kron(D, p) * p ** j for p in ps}
sol = in_span(cols)
print(f"all twists together ({len(cols)} columns):", 'EXACT' if sol else 'NOT in span')
for k in (3, 5, 6):
    c2 = dict(cols)
    for lab in forms:
        if lab_k(lab) == k: c2.update(orbit_cols(lab))
    if len(c2) < len(ps) - 3:
        sol = in_span(c2)
        print(f"twists + weight {k} orbits ({len(c2)} columns):", 'EXACT ' + str({kk: str(v) for kk, v in sol.items() if v}) if sol else 'NOT in span')
