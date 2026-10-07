import numpy as np
from fractions import Fraction as F
from itertools import combinations
import importlib.util, sys
spec = importlib.util.spec_from_file_location("ov", "/private/tmp/claude-501/-Users-e-Documents-GitHub-math/f1b41b3c-5f00-4184-8c88-059bb9352d19/scratchpad/oai2/lds/calc/c4_overlap.py")
src = open(spec.origin).read().split("dl = F(14, 183)")[0]
ns = {}; exec(src, ns)
pair = ns['pair']
dl = F(14, 183)
t = (np.arange(4_000_000) + 0.5)/4_000_000
def bad(v): x = (v*t) % 1.0; return np.minimum(x, 1-x) < float(dl)
for q, r in ((12, 182), (6, 4), (1, 7), (13, 182)):
    mc = np.mean(bad(q) & bad(r))
    print(q, r, "exact", pair(q, r, dl), float(pair(q, r, dl)), "grid", mc)
S = list(range(1, 13)) + [182]
P = sum(pair(q, r, dl) for q, r in combinations(S, 2)); sB = 2*dl*13
print("deep well: sum pair overlaps =", P, "=", float(P))
CE = sB**2/(sB + 2*P)
print("Chung-Erdos lower bound for |union| =", CE, "=", float(CE))
# grid check of union and of the sum of pair overlaps
Bs = [bad(v) for v in S]
cnt = np.sum(Bs, axis=0)
print("grid: |union| =", np.mean(cnt > 0), " sum pairs =", np.mean(cnt*(cnt-1)/2), " mean multiplicity =", np.mean(cnt))
# multiplicity distribution: covering is 'balanced' (higher-order structure)
vals, counts = np.unique(cnt, return_counts=True)
print("multiplicity distribution:", dict(zip(vals.tolist(), np.round(counts/len(t), 4).tolist())))
import random
random.seed(1); R = sorted(random.sample(range(1, 200), 13)); cntR = np.sum([bad(v) for v in R], axis=0)
vals, counts = np.unique(cntR, return_counts=True)
print("random 13-set multiplicity distribution:", dict(zip(vals.tolist(), np.round(counts/len(t), 4).tolist())))
print("Poisson(1.989) reference:", {k: round(np.exp(-1.9891)*1.9891**k/np.math.factorial(k), 4) for k in range(6)} if hasattr(np, 'math') else None)
