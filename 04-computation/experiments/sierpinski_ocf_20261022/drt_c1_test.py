# Test of HYP-9371 (C1) and (C2) on doubly regular tournaments not used to find them.  Paley switching families
# (drt_switch.family) of P_27 (over GF(27) = GF(3)[x]/(x^3 - x - 1)) and P_31 are single classes, so non-Paley DRTs
# come from the skew-Hadamard doubling K -> [[K, K], [-K^T, K^T]]: Paley order 12 doubled gives order 24 (DRT(23)),
# Paley order 8 doubled twice gives order 32 (DRT(31)).  Each skew-Hadamard K = I + B is normalised at every row r
# (r made a source, then deleted), one
# representative per nauty class.  For each non-Paley class against the Paley tournament of the same order:
#   (C1) is Delta(5-census)/Delta TT5 the integer direction found at n = 19, 23 (drt_census5_direction.out)?
#   (C2) is n^8 [z^8] Delta X_3 = 7168 Delta TT5 (four-loop engine, two-matrix reduction)?
# Census by complete isomorphism type (drt_census5_direction.census5_full), after audit_loop4/REPORT.md.
# Output: drt_c1_test.out.  mac-mini-2026-10-22-loop4.
import sys, itertools, subprocess, json
from collections import Counter
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
import loop_engine4 as LE4
from drt_switch import family, d6

from drt_census5_direction import census5_full
# per-TT5 direction of the 5-census (12 isomorphism types), from the n = 19 and 23 data (drt_census5_direction.py)
V = {k: Fraction(v) for k, v in json.load(open('drt_census5_direction.json')).items()}

def gf27():
    els = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]      # a + b x + c x^2
    def mul(u, v):
        p = [0] * 5
        for i in range(3):
            for j in range(3): p[i + j] += u[i] * v[j]
        for d in (4, 3):                                                          # x^3 = x + 1
            c = p[d]; p[d] = 0; p[d - 3] += c; p[d - 2] += c
        return tuple(x % 3 for x in p[:3])
    sq = {mul(u, u) for u in els if u != (0, 0, 0)}
    n = 27; S = np.zeros((n, n), dtype=int)
    for i, u in enumerate(els):
        for j, v in enumerate(els):
            if i != j: S[i, j] = 1 if tuple((b - a) % 3 for a, b in zip(u, v)) in sq else -1
    return S

def paley_prime(q):
    sq = {(t * t) % q for t in range(1, q)}
    return np.array([[0 if i == j else (1 if (j - i) % q in sq else -1) for j in range(q)] for i in range(q)])

census5 = census5_full                                  # complete isomorphism-type census

def delta8(S, K=64, R=1.0):
    n = S.shape[0]; cores = LE4.load_cores(3)
    G2 = [LE4.drt_gamma2(S, 3, i) for i in range(len(cores))]
    zz = R * np.exp(2j * np.pi * np.arange(K) / K)
    vals = np.array([sum(np.sum(LE4.reduce_W(LE4.weight_tensor(3, i), LE4.drt_coeffs(n, z)) * G2[i]) / cores[i][2]
                         for i in range(len(cores))) for z in zz])
    return (np.mean(vals * np.exp(-2j * np.pi * np.arange(K) * 8 / K)) / R ** 8).real

def bordered(S):
    n = S.shape[0]; B = np.zeros((n + 1, n + 1), dtype=int)
    B[1:, 1:] = S; B[0, 1:] = 1; B[1:, 0] = -1
    return B
def double(B):
    m = B.shape[0]; K = np.eye(m, dtype=int) + B
    K2 = np.block([[K, K], [-K.T, K.T]])
    assert (K2 @ K2.T == 2 * m * np.eye(2 * m, dtype=int)).all() and (K2 + K2.T == 2 * np.eye(2 * m, dtype=int)).all()
    return K2 - np.eye(2 * m, dtype=int)
def normalisations(B):
    m = B.shape[0]; out = []
    for r in range(m):
        d = np.array([B[r, j] if j != r else 1 for j in range(m)])
        Br = d[:, None] * B * d[None, :]
        keep = [j for j in range(m) if j != r]
        out.append(Br[np.ix_(keep, keep)])
    return out

if __name__ == "__main__":
    for name, S in (('P_27', gf27()), ('P_31', paley_prime(31))):
        assert LE4.check_drt(S), name
        c = subprocess.run(['labelg', '-q'], input='\n'.join(d6(D) for D in [S] + family(S)) + '\n', capture_output=True, text=True).stdout.split()
        print(f"{name}: its own switching family has {len(set(c))} nauty class(es)", flush=True)
    B23 = double(bordered(paley_prime(11)))                # order 24
    B31 = double(double(bordered(paley_prime(7))))         # order 32
    for name, S, Bd in (('P_23', paley_prime(23), B23), ('P_31', paley_prime(31), B31)):
        n = S.shape[0]
        assert LE4.check_drt(S), name
        fam = normalisations(Bd)
        canon = subprocess.run(['labelg', '-q'], input='\n'.join(d6(D) for D in [S] + fam) + '\n', capture_output=True, text=True).stdout.split()
        reps = {}
        for D, cf in zip([S] + fam, canon): reps.setdefault(cf, D)
        print(f"{name} and the row normalisations of the doubled Paley skew-Hadamard matrix of order {n + 1}: {len(reps)} nauty classes (Paley included)", flush=True)
        cP = census5(S); d8P = delta8(S); tP = cP['(0, 1, 2, 3, 4)|1']
        for cf, D in reps.items():
            if cf == canon[0]: continue
            assert LE4.check_drt(D)
            cD = census5(D); dt = cD['(0, 1, 2, 3, 4)|1'] - tP
            keys = set(cD) | set(cP) | set(V)
            ratio = {k: Fraction(cD.get(k, 0) - cP.get(k, 0), dt) if dt else None for k in keys}
            c1 = dt != 0 and all(ratio[k] == V.get(k, 0) for k in keys)
            d8 = delta8(D) - d8P
            print(f"   class {cf[:12]}...: Delta TT5 = {dt}; (C1) census direction equal to the n = 19/23 vector: {c1}; "
                  f"n^8 Delta delta_8 / Delta TT5 = {n ** 8 * d8 / dt if dt else float('nan'):.6f} (C2 predicts 7168)", flush=True)
            if not c1 and dt:
                print('      direction:', {k: str(v) for k, v in sorted(ratio.items())})
