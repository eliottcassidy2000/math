# Doubly regular 23-tournaments: do H differences track the first census level at which two classes differ?
# For each nauty class: H (hk64), the 5-vertex census (key = sorted scores + H of the 5-subtournament).
# Usage: drt_census_hierarchy.py HK codes_file [codes_file ...].  mac-mini-2026-10-22-sierp.
import sys, subprocess, itertools
from collections import Counter
import numpy as np
sys.path.insert(0, '.')
from drt_switch import code_to_S, d6
HK = sys.argv[1]
codes = []
for f in sys.argv[2:]:
    codes += [l.strip() for l in open(f) if l.strip()]
q = 23; sq = {(t * t) % q for t in range(1, q)}
codes.append(''.join('1' if (j - i) % q in sq else '0' for i in range(q) for j in range(i + 1, q)))   # P_23 last
canon = subprocess.run(['labelg', '-q'], input='\n'.join(d6(code_to_S(c)) for c in codes) + '\n', capture_output=True, text=True).stdout.split()
uniq = {}
for c, cf in zip(codes, canon):
    if cf not in uniq: uniq[cf] = c
reps = list(uniq.values())
Hs = [int(x) for x in subprocess.run([HK], input='\n'.join(reps) + '\n', capture_output=True, text=True).stdout.split()]
perms5 = list(itertools.permutations(range(5)))
def census5(S):
    n = S.shape[0]; cnt = Counter()
    for X in itertools.combinations(range(n), 5):
        sub = S[np.ix_(X, X)]
        sc = tuple(sorted(int((sub[i] == 1).sum()) for i in range(5)))
        h = sum(1 for p in perms5 if all(sub[p[i], p[i + 1]] == 1 for i in range(4)))
        cnt[(sc, h)] += 1
    return cnt
cens = [census5(code_to_S(c)) for c in reps]
HP = Hs[[k for k, c in enumerate(reps) if c == codes[-1]][0]] if codes[-1] in reps else None
print(f"{len(reps)} distinct classes (P_23 included); H(P_23) = {HP}")
for i in range(len(reps)):
    for j in range(i + 1, len(reps)):
        same5 = cens[i] == cens[j]
        print(f"  classes {i},{j}: same 5-census {same5}; |H_i/H_j - 1| = {abs(Hs[i] / Hs[j] - 1):.3e}")

# Is H/H(P) - 1 linear in the 5-census difference?  Rank of the census-difference matrix and a least-squares fit.
keys = sorted(set(k for c in cens for k in c))
ip = len(reps) - 1 if reps[-1] == codes[-1] else [k for k, c in enumerate(reps) if c == codes[-1]][0]
D = np.array([[cens[i][k] - cens[ip][k] for k in keys] for i in range(len(reps)) if i != ip], float)
y = np.array([Hs[i] / Hs[ip] - 1 for i in range(len(reps)) if i != ip])
print('5-census types:', len(keys), '; rank of census differences:', np.linalg.matrix_rank(D))
U, sv, Vt = np.linalg.svd(D, full_matrices=False)
r = int(np.sum(sv > 1e-9 * sv[0]))
coords = D @ Vt[:r].T
for k in range(1, r + 1):
    X = coords[:, :k]
    c, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ c
    print(f"  fit H/H(P) - 1 on the first {k} census direction(s): max |resid| = {np.max(np.abs(y - pred)):.3e} (data range {np.min(y):.3e}..{np.max(y):.3e})")
print('  census coordinates of each class (first direction):', np.round(coords[:, 0], 2), ' y:', [f"{v:.3e}" for v in y])

# Linear model with intercept in the transitive-5-set count TT5 (the free 5-census parameter), non-Paley classes only.
tt = np.array([cens[i][((0, 1, 2, 3, 4), 1)] for i in range(len(reps))], float)
Hr = np.array(Hs, float)
nonP = [i for i in range(len(reps)) if i != ip]
A1 = np.column_stack([np.ones(len(nonP)), tt[nonP]])
c, *_ = np.linalg.lstsq(A1, Hr[nonP], rcond=None)
pred = A1 @ c
print('TT5 counts:', tt.astype(int).tolist(), ' (Paley index', ip, ')')
print(f"  non-Paley fit H = {c[0]:.6e} + {c[1]:.6e} TT5: max relative residual {np.max(np.abs(pred / Hr[nonP] - 1)):.3e} "
      f"(spread of H/H(P): {np.ptp(Hr[nonP] / Hr[ip]):.3e})")
print(f"  Paley off the line by {(Hr[ip] - (c[0] + c[1] * tt[ip])) / Hr[ip]:+.3e} (relative)")
