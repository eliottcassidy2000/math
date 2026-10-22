# Four-loop prediction of the Hamiltonian-path differences among doubly regular 23-tournaments (HYP-9371).
#
# For a doubly regular T (S S^T = nI - J) every propagator is alpha P + gamma Shat (P = I - J/n), so each core
# contributes <W x_e C2(z), Gamma2>/|G| with Gamma2 the 2^E tensor of {P, Shat} homomorphism numbers of the core
# (loop_engine4).  Cores with <= 4 vertices see only the 4-vertex census, which is fixed by n: X_0, X_1, X_2 and the
# <= 4-vertex part of X_3 are universal.  The differences Delta X_3 = X_3(D) - X_3(P_23) split exactly into
#   five-vertex part : 5-vertex cores, plus every 6-vertex core with one pair of vertices merged (one P-edge's identity
#                      part, or one coincidence of the homomorphism), whose non-universal part is the injective 5-sum;
#   six-vertex part  : 6-vertex cores minus their 15 pair merges (the injective 6-sums).
# The five-vertex part is a linear function of the 5-vertex census, which for these classes has rank 1 (TT5).
#
# Predictions of Delta log H = log H(D) - log H(P_23), no fitted parameters:
#   strict third order : Delta X_3(1)  (all other terms through n^-3 are universal);
#   + n^-4 response    : sum_j [1 + (j D_1 + j(j-1)/2)/n] [z^j] Delta X_3, D_1 = 2n/(n+1) (one-loop mean length):
#                        the finite-size (g(L) = 1 + C(L,2)/n + ...) dressing of the four-loop term, which is the
#                        part of the n^-4 difference that four loops determine (Delta X_4(1) is not included);
#   negative control   : the 'resummed' response sum_L g(L) [z^L](Phi_2 Delta X_3) / sum_L g(L) [z^L] Phi_2 with
#                        Phi_2 = exp(X_0 + X_1 + X_2) fails: Phi_2 has no zero near z^2 = -2.2 and a heavy tail at
#                        L ~ n where g(L) = n^L/(n)_L ~ 10^7 (see drt23_anatomy.py for the exact measure).
# Output: drt23_four_loop.out, drt23_four_loop.json.  mac-mini-2026-10-22-loop4.
import sys, json, time, itertools
from math import factorial
import numpy as np
sys.path.insert(0, '.')
import loop_engine4 as LE4
from drt_switch import code_to_S
from mpmath import mp, mpf
mp.dps = 40

cls = json.load(open('drt23_classes.json'))
iP = [c['class'] for c in cls if c['label'] == 'P_23'][0]
n = 23
K, R = 256, 1.5                           # Cauchy contour for the z-series (radius of convergence sqrt(n) = 4.8)
zs = R * np.exp(2j * np.pi * np.arange(K) / K)
LMAX = 40

def coeffs(vals):
    return np.array([(np.mean(vals * np.exp(-2j * np.pi * np.arange(K) * L / K)) / R ** L).real for L in range(LMAX + 1)])

# ---- per-core reduced weight tensors W2(z) on the contour and at z = 1 (tournament-independent)
t0 = time.time()
C2z = [LE4.drt_coeffs(n, z) for z in zs]
C21 = LE4.drt_coeffs(n, 1.0)
def merged_graph(V, edges, u, v):
    m = {}; k = 0
    for x in range(V):
        if x == v: continue
        m[x] = k; k += 1
    m[v] = m[u]
    return V - 1, [(m[a], m[b]) for a, b in edges]

red = {}            # (excess, i) -> (W2 at z = 1, [W2 on contour])
for e in (1, 2, 3):
    for i, (V, edges, aut) in enumerate(LE4.load_cores(e)):
        W = LE4.weight_tensor(e, i)
        red[(e, i)] = (LE4.reduce_W(W, C21), np.array([LE4.reduce_W(W, C) for C in C2z]))
print(f"reduced weight tensors for {len(red)} cores in {time.time() - t0:.0f}s", flush=True)

# ---- per-class homomorphism tensors
def class_data(S):
    P = np.eye(n) - np.ones((n, n)) / n; Sh = S / n
    B = np.array([P, Sh])
    out = {}
    for e in (1, 2, 3):
        for i, (V, edges, aut) in enumerate(LE4.load_cores(e)):
            G = LE4.network(e, i).contract(B)
            merged = []
            if e == 3 and V == 6:
                for u, v in itertools.combinations(range(V), 2):
                    Vm, em = merged_graph(V, edges, u, v)
                    merged.append(LE4.Network(Vm, em).contract(B))
            out[(e, i)] = (G, np.sum(merged, axis=0) if merged else None)
    return out

def series_and_value(data, keys, part='all'):
    """returns (value at z = 1, coefficients) of sum over the given cores; part: 'all', 'five' (pair merges of
    6-vertex cores + all cores with V < 6) or 'six' (6-vertex cores minus their pair merges)"""
    v1 = 0.0; vz = np.zeros(K, dtype=complex)
    for key in keys:
        e, i = key; V, edges, aut = LE4.load_cores(e)[i]
        G, Gm = data[key]
        if part == 'five': T = Gm if V == 6 else G
        elif part == 'six': T = (G - Gm) if V == 6 else None
        else: T = G
        if T is None: continue
        W1, Wz = red[key]
        v1 += np.sum(W1 * T) / aut
        vz += np.tensordot(Wz, T, axes=T.ndim) / aut
    return v1, coeffs(vz)

t0 = time.time()
D = {}
for c in cls:
    S = code_to_S(c['code']); assert LE4.check_drt(S)
    D[c['class']] = class_data(S)
print(f"homomorphism tensors for {len(cls)} classes in {time.time() - t0:.0f}s", flush=True)

keys = {e: [(e, i) for i in range(len(LE4.load_cores(e)))] for e in (1, 2, 3)}
res = {}
for c in cls:
    k = c['class']
    r = {}
    for e in (1, 2):
        r[f'X{e}'] = series_and_value(D[k], keys[e])
    r['X3'] = series_and_value(D[k], keys[3])
    r['X3_five'] = series_and_value(D[k], keys[3], 'five')
    r['X3_six'] = series_and_value(D[k], keys[3], 'six')
    r['X3_small'] = series_and_value(D[k], [key for key in keys[3] if LE4.load_cores(3)[key[1]][0] <= 4])
    res[k] = r

# ---- universality checks
print("\nUniversality (max over classes of |X(D) - X(P_23)|, value at z = 1 and largest series coefficient):")
for name in ('X1', 'X2', 'X3_small'):
    dv = max(abs(res[k][name][0] - res[iP][name][0]) for k in res)
    ds = max(np.max(np.abs(res[k][name][1] - res[iP][name][1])) for k in res)
    print(f"   {name:8s}: value at z=1 {res[iP][name][0]:+.10e}, max class difference {dv:.1e} (series {ds:.1e})")
print(f"   X3 (P_23) = {res[iP]['X3'][0]:+.10e};  n^3 X3 = {n ** 3 * res[iP]['X3'][0]:+.4f}")

# ---- universal second-order measure Phi_2 = exp(X_0 + X_1 + X_2)
X0c = np.zeros(LMAX + 1)
for m in range(1, LMAX // 2 + 1, 2): X0c[2 * m] = (n - 1) / (m * n ** m)
X12c = res[iP]['X1'][1] + res[iP]['X2'][1]
def exp_series(a):
    b = np.zeros_like(a); b[0] = np.exp(a[0])
    for L in range(1, len(a)):
        b[L] = sum(k * a[k] * b[L - k] for k in range(1, L + 1)) / L
    return b
phi2 = exp_series(X0c + X12c)
g = np.array([float(mpf(n) ** L / mp.ff(n, L)) for L in range(n)])           # L <= n - 1
kappa2 = float(np.sum(phi2[:n] * g))
kapP = float(mpf(cls[iP]['H']) * mpf(2) ** (n - 1) / mp.factorial(n))
print(f"\nResummed second order: sum_(L<=22) g(L) phi_2(L) = {kappa2:.8f}; exact kappa(P_23) = {kapP:.8f} (ratio - 1 = {kappa2 / kapP - 1:+.3e})")
w = np.array([np.sum(phi2[:n - j] * g[j:]) / kappa2 for j in range(n)])
print("weights w_j = E[g(M + j)]/E[g(M)] under Phi_2:", ' '.join(f"w_{j}={w[j]:.3f}" for j in range(0, 17, 2)))

def resummed(co):
    return float(np.sum(w * co[:n]))
D1 = 2 * n / (n + 1)
wnext = np.array([1 + (j * D1 + j * (j - 1) / 2) / n for j in range(LMAX + 1)])
def dressed(co):
    return float(np.sum(wnext * co))
print("negative control: the Phi_2 resummation weights above are unusable (sum g phi_2 is 80 times kappa).")

# ---- table
HP = cls[iP]['H']
print("\nPer class (relative to P_23):  exact = log H(D) - log H(P_23);  strict = Delta X_3(1) = five + six;  dressed = strict + its n^-4 finite-size response")
print(f"{'class':8s} {'TT5':>5s} {'exact':>12s} {'strict':>12s} {'five':>12s} {'six':>12s} {'dressed':>12s} {'exact/strict':>12s}   lowest Delta coefficients z^8, z^10, z^12")
rows = []
for c in cls:
    k = c['class']
    ex = float(mp.log(mpf(c['H']) / mpf(HP)))
    d = {name: (res[k][name][0] - res[iP][name][0], res[k][name][1] - res[iP][name][1]) for name in ('X3', 'X3_five', 'X3_six')}
    row = {'class': k, 'label': c['label'], 'TT5': c['TT5'], 'dTT5': c['TT5'] - cls[iP]['TT5'], 'exact': ex,
           'strict': d['X3'][0], 'five': d['X3_five'][0], 'six': d['X3_six'][0],
           'resummed': resummed(d['X3'][1]), 'dressed': dressed(d['X3'][1]),
           'coeffs': d['X3'][1][:24].tolist(), 'coeffs_five': d['X3_five'][1][:24].tolist(), 'coeffs_six': d['X3_six'][1][:24].tolist()}
    rows.append(row)
    print(f"{c['label']:8s} {c['TT5']:5d} {ex:+.5e} {row['strict']:+.5e} {row['five']:+.5e} {row['six']:+.5e} {row['dressed']:+.5e} {ex / row['strict'] if row['strict'] else float('nan'):12.3f}   "
          + ' '.join(f"{x:+.3e}" for x in d['X3'][1][8:13:2]))

# ---- TT5 laws over the non-Paley classes
nonP = [r for r in rows if r['class'] != iP]
tt = np.array([r['TT5'] for r in nonP], float)
def fit(y):
    A = np.column_stack([np.ones(len(tt)), tt]); c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c, np.max(np.abs(y - A @ c))
print("\nLinear laws in TT5 over the 12 non-Paley classes (slope per transitive 5-set; max residual):")
for name in ('exact', 'strict', 'five', 'six', 'dressed', 'resummed'):
    y = np.array([r[name] for r in nonP])
    (a, b), mr = fit(y)
    print(f"   {name:9s}: slope {b:+.4e} (x H(P_23) = {b * HP:.4e} per TT5), max residual {mr:.2e}, line at Paley's TT5 = 1518: {a + b * 1518:+.3e}")
json.dump({'n': n, 'rows': rows, 'w': w.tolist(), 'kappa2': kappa2, 'kappa_P23': kapP,
           'X3_P23': res[iP]['X3'][0], 'X3_P23_coeffs': res[iP]['X3'][1][:24].tolist()}, open('drt23_four_loop.json', 'w'), indent=1)
