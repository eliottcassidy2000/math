# AUDIT of the claims added to HYP-9371 during the audit (doubled-Paley DRT(23) class; (C1) at n = 31):
#  (a) rebuild the doubled-Paley DRT(23) independently (bordered P_11 -> skew-Hadamard K of order 12 ->
#      [[K, K], [-K^T, K^T]] -> normalise at every row, delete it), certify DRT and nauty class = drt23_doubled.json;
#  (b) own exact forest numbers (aud_fwdpoly, 2 primes): H, Delta N_8 vs 7168 Delta TT5, per-L shares, relative H;
#      the annealed TT5 law's out-of-sample error; own 5-census (true iso classes) direction;
#  (c) n = 31: twice-doubled Paley skew-Hadamard matrix of order 32 from P_7; own nauty classes, TT5 and 5-census
#      direction against P_31.
import sys, json, subprocess, itertools
from fractions import Fraction
from math import factorial
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE + '/audit_loop4')
from aud_task1 import F_mod, lagrange_coeffs, crt2, forest_N, P1, P2, code_adj, tt

def paley_S(q):
    sq = {(t * t) % q for t in range(1, q)}
    return np.array([[0 if i == j else (1 if (j - i) % q in sq else -1) for j in range(q)] for i in range(q)])
def skew_hadamard_from_drt(S):
    n = S.shape[0]; B = np.zeros((n + 1, n + 1), dtype=int); B[1:, 1:] = S; B[0, 1:] = 1; B[1:, 0] = -1
    K = np.eye(n + 1, dtype=int) + B
    assert (K @ K.T == (n + 1) * np.eye(n + 1, dtype=int)).all()
    return K
def doubled(K):
    m = K.shape[0]; K2 = np.block([[K, K], [-K.T, K.T]])
    assert (K2 @ K2.T == 2 * m * np.eye(2 * m, dtype=int)).all() and (K2 + K2.T == 2 * np.eye(2 * m, dtype=int)).all()
    return K2
def drts_from(K):
    m = K.shape[0]; B = K - np.eye(m, dtype=int); out = []
    for r in range(m):
        d = np.array([B[r, j] if j != r else 1 for j in range(m)])
        Br = d[:, None] * B * d[None, :]
        keep = [j for j in range(m) if j != r]
        D = Br[np.ix_(keep, keep)]
        assert (D @ D.T == (m - 1) * np.eye(m - 1, dtype=int) - 1).all()
        out.append(D)
    return out
def d6(S):
    n = S.shape[0]; bits = [1 if S[i, j] == 1 else 0 for i in range(n) for j in range(n)]
    while len(bits) % 6: bits.append(0)
    s = '&' + chr(63 + n)
    for k in range(0, len(bits), 6):
        v = 0
        for b in bits[k:k + 6]: v = 2 * v + b
        s += chr(63 + v)
    return s
def canon(Ss):
    return subprocess.run(['labelg', '-q'], input='\n'.join(d6(S) for S in Ss) + '\n', capture_output=True, text=True).stdout.split()
def S_to_code(S):
    n = S.shape[0]; return ''.join('1' if S[i, j] == 1 else '0' for i in range(n) for j in range(i + 1, n))
# own 5-census by isomorphism class (as in aud_census5.py)
PAIRS = [(i, j) for i in range(5) for j in range(i + 1, 5)]
perms = list(itertools.permutations(range(5)))
def adj10(c):
    A = [[0] * 5 for _ in range(5)]
    for b, (i, j) in enumerate(PAIRS):
        if c >> b & 1: A[i][j] = 1
        else: A[j][i] = 1
    return A
def code10(A): return sum(1 << b for b, (i, j) in enumerate(PAIRS) if A[i][j])
CAN = {c: min(code10([[adj10(c)[p[i]][p[j]] for j in range(5)] for i in range(5)]) for p in perms) for c in range(1024)}
CLS = sorted(set(CAN.values()))
def census5(S):
    n = S.shape[0]; cnt = {c: 0 for c in CLS}
    for X in itertools.combinations(range(n), 5):
        c = 0
        for b, (i, j) in enumerate(PAIRS):
            if S[X[i], X[j]] == 1: c |= 1 << b
        cnt[CAN[c]] += 1
    return cnt
TTc = 0   # canonical code of the transitive 5-tournament
VEC = None

cls = json.load(open(f'{HERE}/drt23_classes.json'))
dj = json.load(open(f'{HERE}/drt23_doubled.json'))
SP23 = paley_S(23)
fam = drts_from(doubled(skew_hadamard_from_drt(paley_S(11))))
cf = canon([SP23] + fam + [np.array(SP23)])
cands = {c: D for c, D in zip(cf[1:-1], fam) if c != cf[0]}
print(f"(a) doubled-Paley family: {len(fam)} row normalisations, {len(set(cf[1:-1]))} nauty classes, non-Paley classes: {len(cands)}")
Dd = list(cands.values())[0]
from_json = np.zeros((23, 23), dtype=int); k = 0
for i in range(23):
    for j in range(i + 1, 23):
        s = 1 if dj['code'][k] == '1' else -1; from_json[i, j] = s; from_json[j, i] = -s; k += 1
cJ = canon([from_json, Dd])
print(f"    rebuilt class equals the class of drt23_doubled.json: {cJ[0] == cJ[1]}; JSON tournament doubly regular: "
      f"{(from_json @ from_json.T == 23 * np.eye(23, dtype=int) - 1).all()}; not among the 13 classes: {cJ[0] not in {c['canon'] for c in cls}}")
code = dj['code']; n = 23
v1 = F_mod(code, P1, range(n)); v2 = F_mod(code, P2, range(n))
f = [crt2(a, b) for a, b in zip(lagrange_coeffs([v1[x] for x in range(n)], P1), lagrange_coeffs([v2[x] for x in range(n)], P2))]
N = forest_N(f, n); H = f[0]
nn, out = code_adj(code); T5 = tt(5, (1 << n) - 1, out)
own = json.load(open(f'{HERE}/audit_loop4/aud_task1_N.json'))
iP = [k for k, v in own.items() if v['label'] == 'P_23'][0]
NP = [int(x) for x in own[iP]['N']]; TP = own[iP]['TT5']
HP = int([c for c in cls if c['label'] == 'P_23'][0]['H'])
print(f"(b) own: H = {H} (json {dj['H']}, equal {str(H) == str(dj['H'])}), f palindromic {f == f[::-1]}, sum f = 23! {sum(f) == factorial(23)}; TT5 = {T5}; "
      f"own N equals json N: {[str(x) for x in N] == [str(x) for x in (dj['N'] if isinstance(dj['N'], list) else [])]}")
dT = T5 - TP
print(f"    Delta TT5 = {dT}; Delta N_2,4,6 = {[N[L] - NP[L] for L in (2, 4, 6)]}; Delta N_8 = {N[8] - NP[8]} = 7168 x {Fraction(N[8] - NP[8], 7168)}")
rel = Fraction(H, HP) - 1
print(f"    H/H(P_23) - 1 = {float(rel):+.6e}")
g = [Fraction(n ** L, factorial(n) // factorial(n - L)) for L in range(n)]
kapP = sum(Fraction(NP[L], n ** L) * g[L] for L in range(n))
pl = [float(Fraction(N[L] - NP[L], n ** L) * g[L] / kapP) for L in range(n)]
print("    shares: " + ', '.join(f"L{L} {pl[L] / float(rel):.1%}" for L in (8, 10, 12, 14)) + f", L>=16 {sum(pl[16:]) / float(rel):.2%}")
nonP = [c for c in cls if c['label'] != 'P_23']
A = np.array([[1.0, c['TT5']] for c in nonP]); y = np.array([float(c['H']) for c in nonP])
(a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
print(f"    annealed TT5 law (H-units) predicts {a + b * T5:.6e}: relative error {(a + b * T5) / H - 1:+.3e}")
# log-units law (as in Section 3): exact log ratio vs line
yl = np.array([float(np.log(int(c['H']) / HP)) for c in nonP]); (al, bl), *_ = np.linalg.lstsq(A, yl, rcond=None)
print(f"    annealed TT5 law in log-units predicts log H/H(P) = {al + bl * T5:+.5e}; exact {float(np.log(H / HP)):+.5e}; difference {al + bl * T5 - float(np.log(H / HP)):+.2e}")
print(f"    strict and dressed from the affine laws (slopes 1.00612e-7, 2.99720e-7 per TT5, Paley on the line): {1.00612e-7 * dT:+.4e}, {2.99720e-7 * dT:+.4e}; "
      f"exact/strict {float(np.log(H / HP)) / (1.00612e-7 * dT):.3f}, dressed/exact {2.99720e-7 * dT / float(np.log(H / HP)):.3f}")
cP = census5(SP23); cD = census5(from_json)
dirD = {c: Fraction(cD[c] - cP[c], cD[TTc] - cP[TTc]) for c in CLS}
print(f"    own 5-census direction of the doubled class (true iso classes, canonical codes {CLS}): {[str(dirD[c]) for c in CLS]}")

print("(c) n = 31: twice-doubled Paley skew-Hadamard matrix of order 32 (from P_7)")
SP31 = paley_S(31)
fam31 = drts_from(doubled(doubled(skew_hadamard_from_drt(paley_S(7)))))
cf31 = canon([SP31] + fam31)
reps = {}
for c, D in zip(cf31[1:], fam31): reps.setdefault(c, D)
print(f"    {len(fam31)} row normalisations -> {len(reps)} nauty classes; Paley among them: {cf31[0] in reps}")
cP31 = census5(SP31)
for c, D in reps.items():
    if c == cf31[0]: continue
    cD = census5(D); dt = cD[TTc] - cP31[TTc]
    dirv = [str(Fraction(cD[x] - cP31[x], dt)) if dt else 'nan' for x in CLS]
    conv = canon([D, -D]);
    print(f"    class {c[:14]}...: Delta TT5 = {dt}; direction {dirv}; self-converse: {conv[0] == conv[1]}", flush=True)
print(f"    reference direction (n = 19, 23; aud_census5.out): ['1', '-1', '1', '-1', '-3', '1', '5', '-1', '1', '3', '-7', '1']")
