# AUDIT task 3: hostile numerical tests of the Section-2 Proposition of HYP-9371 (X_3 of a doubly regular T is
# affine in the 5-vertex census) on the DRT(23) classes.
#  (A) DRT algebra: S^2 = J - nI (so Shat^2 = -P/n), S 1 = 0; all 6-vertex excess-3 cores are cubic.
#  (B) Parity: split X_3(z) = E(z) + O(z) by the parity of #Shat in the {P, Shat} labelling; O must vanish identically
#      (X_3 is even in z).  Evaluated at z = 1, 0.7, 1.3 (complex z = 0.9 e^(0.3i) too) for all 13 classes, with the
#      scale of the individual terms; per-core O_core(1) for the 6-vertex cores.
#  (C) Injective 6-sums inj6[t] (Moebius inversion over the 203 partitions of the 6 core vertices) for every labelling
#      t in {P, Shat}^9 of every 6-vertex core: D23_1 vs D23_2 (equal 5-census, different 6-census) must agree for every
#      t with #Shat <= 8 (step 3: a vertex of Shat-degree <= 2 can be summed out) and may differ only for t = all-Shat.
#  (D) The all-Shat coefficient of each 6-vertex core (reduced weight tensor entry W2[1,...,1](z)); the genuinely
#      6-vertex part sum_(V=6) <W2(1), inj6>/|G| for all 13 classes: affine in TT5 with Paley on the line?
import sys, itertools, time, json
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE)
import loop_engine4 as LE4
from drt_switch import code_to_S

t0 = time.time()
cls = json.load(open(f'{HERE}/drt23_classes.json'))
n = 23
cores = LE4.load_cores(3)
six = [i for i, (V, e, a) in enumerate(cores) if V == 6]
print(f"(A) {len(cores)} excess-3 cores; 6-vertex cores: {len(six)}; all cubic: "
      f"{all(sorted(np.bincount(np.array(cores[i][1]).ravel(), minlength=6)) == [3] * 6 for i in six)}")
for c in cls:
    S = code_to_S(c['code'])
    assert np.array_equal(S @ S, np.ones((n, n), dtype=int) - n * np.eye(n, dtype=int)) and not S.sum(axis=1).any()
print("    all 13 classes: S^2 = J - nI and S 1 = 0 (so Shat^2 = -P/n, Shat P = Shat): True")

def popcount_mask(E):
    idx = np.indices((2,) * E).reshape(E, -1).sum(axis=0).reshape((2,) * E)
    return idx

zlist = [1.0, 0.7, 1.3, 0.9 * np.exp(0.3j)]
W2 = {}   # (i, z) -> reduced weight tensor
for i in range(len(cores)):
    W = LE4.weight_tensor(3, i)
    for z in zlist:
        W2[(i, z)] = LE4.reduce_W(W, LE4.drt_coeffs(n, z))
print(f"    reduced weight tensors at {len(zlist)} z values: {time.time() - t0:.0f}s", flush=True)

def gamma2_all(S):
    Sh = S / n; P = np.eye(n) - np.ones((n, n)) / n
    B = np.array([P, Sh])
    return {i: LE4.network(3, i).contract(B) for i in range(len(cores))}

G2 = {}
for c in cls:
    G2[c['class']] = gamma2_all(code_to_S(c['code']))
print(f"    Gamma2 tensors for 13 classes: {time.time() - t0:.0f}s", flush=True)

print("\n(B) parity split X_3 = E + O (O must vanish):")
for c in cls:
    k = c['class']; line = []
    for z in zlist:
        Ev = 0; Ov = 0; scale = 0
        for i, (V, edges, aut) in enumerate(cores):
            T = W2[(i, z)] * G2[k][i] / aut
            m = popcount_mask(len(edges)) % 2
            Ev += T[m == 0].sum(); Ov += T[m == 1].sum(); scale = max(scale, np.abs(T).max())
        line.append(f"z={z:.3g}: E={Ev:+.10e} O={abs(Ov):.1e} (max |term| {scale:.1e})")
    print(f"   {c['label']:7s} " + '; '.join(line), flush=True)
k0 = 0
print("   per 6-vertex core at z = 1 (class D23_0): O_core = odd-#Shat part, A_core = all-Shat term W2[1..1] Gamma2[1..1]/|G|, W2[1..1](1):")
for i in six:
    V, edges, aut = cores[i]
    T = W2[(i, 1.0)] * G2[k0][i] / aut
    m = popcount_mask(len(edges)) % 2
    print(f"      core {i}: |G| = {aut:4d}, O_core = {T[m == 1].sum():+.3e}, E_core = {T[m == 0].sum():+.3e}, A_core = {T[(1,) * 9]:+.3e}, W2[all-Shat](1) = {W2[(i, 1.0)][(1,) * 9]:+.6e}")

# (C) injective 6-sums by Moebius inversion over set partitions of the 6 core vertices
def set_partitions(items):
    if not items: yield []; return
    first, rest = items[0], items[1:]
    for part in set_partitions(rest):
        for j in range(len(part)):
            yield part[:j] + [[first] + part[j]] + part[j + 1:]
        yield [[first]] + part
from math import factorial
PARTS = list(set_partitions(list(range(6))))
assert len(PARTS) == 203
NETS = {}
def inj6(i, S):
    V, edges, aut = cores[i]
    Sh = S / n; P = np.eye(n) - np.ones((n, n)) / n
    B = np.array([P, Sh])
    if i not in NETS:
        lst = []
        for part in PARTS:
            mu = 1
            for b in part: mu *= (-1) ** (len(b) - 1) * factorial(len(b) - 1)
            lab = {}
            for bi, b in enumerate(part):
                for v in b: lab[v] = bi
            lst.append((mu, LE4.Network(len(part), [(lab[a], lab[b]) for a, b in edges])))
        NETS[i] = lst
    tot = np.zeros((2,) * len(edges))
    for mu, net in NETS[i]:
        tot += mu * net.contract(B)
    return tot

print("\n(C) injective 6-vertex sums, D23_1 vs D23_2 (identical 5-census):")
I1 = {i: inj6(i, code_to_S(cls[1]['code'])) for i in six}
I2 = {i: inj6(i, code_to_S(cls[2]['code'])) for i in six}
print(f"    computed in {time.time() - t0:.0f}s", flush=True)
for i in six:
    E = len(cores[i][1]); cnt = popcount_mask(E)
    d = np.abs(I1[i] - I2[i]); sc = np.abs(I1[i]).max()
    lo = d[cnt <= 8].max(); hi = d[cnt == 9].max()
    print(f"    core {i}: max |inj6(D23_1) - inj6(D23_2)| over the 511 labellings with #Shat <= 8: {lo:.1e}; all-Shat labelling: "
          f"{hi:.3e} (inj6 all-Shat = {I1[i][(1,) * 9]:+.6e} vs {I2[i][(1,) * 9]:+.6e}); scale {sc:.1e}")

print("\n(D) genuinely 6-vertex part sum_(V=6) <W2(1), inj6>/|G| over all 13 classes (and its all-Shat piece):")
rows = []
for c in cls:
    S = code_to_S(c['code'])
    I = I1 if c['class'] == 1 else I2 if c['class'] == 2 else {i: inj6(i, S) for i in six}
    tot = sum(np.sum(W2[(i, 1.0)] * I[i]) / cores[i][2] for i in six)
    allS = sum(W2[(i, 1.0)][(1,) * 9] * I[i][(1,) * 9] / cores[i][2] for i in six)
    rows.append((c['label'], c['TT5'], tot, allS))
    print(f"    {c['label']:7s} TT5 {c['TT5']}: six-part {tot:+.12e}; all-Shat piece {allS:+.3e}", flush=True)
tt = np.array([r[1] for r in rows if r[0] != 'P_23'], float); y = np.array([r[2] for r in rows if r[0] != 'P_23'])
A = np.column_stack([np.ones(len(tt)), tt]); (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
yP = [r[2] for r in rows if r[0] == 'P_23'][0]
print(f"    fit over 12 non-Paley classes: slope {b:+.6e} per TT5, max residual {np.max(np.abs(y - A @ [a, b])):.1e}; "
      f"P_23 minus line at TT5 = 1518: {yP - (a + b * 1518):+.1e}")
print(f"done [{time.time() - t0:.0f}s]")
