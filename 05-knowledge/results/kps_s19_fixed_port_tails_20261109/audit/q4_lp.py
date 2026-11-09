"""S19 audit, Q4: LP lower bound for the minimum core cover, with an exactly verified dual certificate.

Input: the sparse file of q4.c (each line: target value, then the indices of the core states whose tails reach it).
Dropping targets only weakens the covering constraints, so the optimum over any subset of targets is a lower bound
for the full minimum. For a subset of targets (those with at most K covering cores) solve the LP relaxation
  min sum x_c  s.t.  sum_{c covers t} x_c >= 1,  0 <= x_c,
take its dual y >= 0, rescale it so that every core column sum is <= 1 (checked in exact rational arithmetic), and
report ceil(sum y) as a rigorous lower bound (weak LP duality). Optionally run a time-limited MILP for comparison.
Usage: python q4_lp.py SPARSE K [MILP_SECONDS [CERT_OUT]]
"""
import sys, math
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds
from scipy.sparse import csr_matrix

path, K = sys.argv[1], int(sys.argv[2])
milp_s = float(sys.argv[3]) if len(sys.argv) > 3 else 0
cert_out = sys.argv[4] if len(sys.argv) > 4 else None
rows, targets = [], []
for line in open(path):
    f = list(map(int, line.split()))
    if len(f) - 1 <= K:
        targets.append(f[0]); rows.append(f[1:])
cols = sorted(set(c for r in rows for c in r))
cidx = {c: i for i, c in enumerate(cols)}
data, ri, ci = [], [], []
for i, r in enumerate(rows):
    for c in r:
        ri.append(i); ci.append(cidx[c]); data.append(1.0)
A = csr_matrix((data, (ri, ci)), shape=(len(rows), len(cols)))
print(f"targets with <= {K} covering cores: {len(rows)}; core states involved: {len(cols)}; nonzeros: {len(data)}")
res = linprog(np.ones(len(cols)), A_ub=-A, b_ub=-np.ones(len(rows)), bounds=(0, None), method="highs")
print(f"LP status: {res.status} ({res.message}); LP optimum {res.fun:.6f}")
y = np.maximum(-res.ineqlin.marginals, 0.0)     # dual of A x >= 1
# exact verification of a feasible dual: y_t as fractions (rounded down to 1e-9), column sums <= 1 after scaling
yq = [Fraction(int(math.floor(v * 1e9)), 10 ** 9) for v in y]
colsum = [Fraction(0)] * len(cols)
for i, r in enumerate(rows):
    if yq[i]:
        for c in r: colsum[cidx[c]] += yq[i]
mx = max(colsum)
scale = max(Fraction(1), mx)
total = sum(yq) / scale
print(f"dual certificate: sum y = {float(sum(yq)):.6f}, max core column sum = {float(mx):.9f}, "
      f"rescaled certified bound = {float(total):.6f} -> every cover needs at least {math.ceil(total)} core states")
if cert_out:
    with open(cert_out, "w") as fh:
        fh.write(f"# dual certificate: target value, y_t (exact fraction); sum = {total} ; max core column sum before scaling = {mx}" + "\n")
        for t, v in zip(targets, yq):
            if v: fh.write(f"{t} {v}" + "\n")
if milp_s > 0:
    r = milp(np.ones(len(cols)), constraints=LinearConstraint(A, lb=np.ones(len(rows)), ub=np.inf),
             integrality=np.ones(len(cols)), bounds=Bounds(0, 1), options={"time_limit": milp_s, "disp": False})
    db = getattr(r, "mip_dual_bound", None)
    print(f"MILP (time limit {milp_s}s): status {r.status} ({r.message}); best cover found {r.fun}; dual bound {db}")
