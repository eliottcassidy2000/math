#!/usr/bin/env python3
"""
t511_witness_lift_kps2.py — kind-pasteur-2026-06-09-S2, Branch I: T511 lower bound

Direct has_TT(T511, 31) died (memory). Structural alternative:
T511 contains D(D(T127)) induced (two-level border-twin identity, verified in
tower_step_structure_kps2). D is monotone on induced subtournaments, so for any
X <= T127 induced, D(D(X)) <= D(D(T127)) <= T511 induced.

Take C15 = a maximum transitive chain of T127 (witness reconstruction from the
trans memo), and X = C15 (+ extra vertices). Compute trans(D(D(X))) exactly on
60-68 vertices. Any value >= 31 certifies trans(T511) >= 31 (the prediction of
the two-step law trans(level k) = 2*trans(level k-2) + 1).

Output: 05-knowledge/results/t511_witness_lift_kps2.out
"""
import sys, time
import numpy as np

sys.setrecursionlimit(100000)
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\04-computation')
from skew_doubling_core_kps1 import D_skew, normalize_first_row, core_tournament, is_skew_hadamard
from erdos_moser_trans_tower_kps2 import TransSolver, masks_of, trans_of, double_S

OUT = open(r'C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\t511_witness_lift_kps2.out',
           'w', encoding='utf-8')
def w(s=''):
    OUT.write(s + '\n'); OUT.flush(); print(s, flush=True)

def witness(sv, S):
    """Reconstruct one maximum transitive chain in S from a solved TransSolver."""
    chain = []
    while S:
        target = sv.trans(S) - 1
        T = S
        found = None
        while T:
            b = T & -T
            v = b.bit_length() - 1
            T ^= b
            if sv.trans(sv.masks[v] & S) == target:
                found = v
                break
        assert found is not None
        chain.append(found)
        S = sv.masks[found] & S
    return chain

def main():
    t_start = time.time()
    w('=== t511_witness_lift_kps2 — certifying trans(T511) >= 31 via D(D(X)) ===')
    S = np.array([[1]], dtype=np.int64)
    while S.shape[0] < 128:
        S = double_S(S)
    T127 = core_tournament(normalize_first_row(S))
    mk, n = masks_of(T127)
    sv = TransSolver(mk)
    full = (1 << n) - 1
    t15 = sv.trans(full)
    C15 = witness(sv, full)
    w(f'trans(T127) = {t15};  witness chain C15 = {C15}')
    # verify chain transitivity
    ok = all(T127[C15[i], C15[j]] for i in range(len(C15)) for j in range(i + 1, len(C15)))
    w(f'witness verifies as transitive chain: {ok}')
    w('')

    def DD_of(X):
        return D_skew(D_skew(X)[0])[0]

    A0 = T127[np.ix_(C15, C15)]
    t0 = time.time()
    tDD, _ = trans_of(DD_of(A0))
    w(f'trans(D(D(TT15 chain))) [60 vertices] = {tDD}   ({time.time()-t0:.1f}s)')
    best = tDD
    w('')

    # X = C15 + one extra vertex (keep trans(X) = 15)
    w('--- X = C15 + v (one extra, trans(X) must stay 15) ---')
    rng = np.random.default_rng(20260609)
    others = [v for v in range(127) if v not in C15]
    rng.shuffle(others)
    tried = 0
    for v in others:
        if tried >= 12 or time.time() - t_start > 1500 or best >= 31:
            break
        Xv = sorted(C15 + [v])
        AX = T127[np.ix_(Xv, Xv)]
        tX, _ = trans_of(AX)
        if tX != 15:
            continue
        tried += 1
        t0 = time.time()
        tDD, _ = trans_of(DD_of(AX))
        w(f'v={v:>3}: trans(X)=15, trans(D(D(X))) [64 vertices] = {tDD}   ({time.time()-t0:.1f}s)')
        best = max(best, tDD)
    w('')

    if best < 31:
        # two extras
        w('--- X = C15 + {u, v} (two extras, trans(X) = 15) ---')
        tried = 0
        for i in range(len(others)):
            if tried >= 8 or time.time() - t_start > 2400 or best >= 31:
                break
            for j in range(i + 1, len(others)):
                if tried >= 8 or best >= 31:
                    break
                Xv = sorted(C15 + [others[i], others[j]])
                AX = T127[np.ix_(Xv, Xv)]
                tX, _ = trans_of(AX)
                if tX != 15:
                    continue
                tried += 1
                t0 = time.time()
                tDD, _ = trans_of(DD_of(AX))
                w(f'u={others[i]:>3} v={others[j]:>3}: trans(D(D(X))) [68 v.] = {tDD} '
                  f'({time.time()-t0:.1f}s)')
                best = max(best, tDD)
    w('')
    w(f'BEST trans(D(D(X))) found = {best}')
    if best >= 31:
        w('=> trans(T511) >= 31 CERTIFIED (D(D(X)) embeds induced in T511).')
        w('   Two-step-law prediction trans(T511) = 31: lower bound confirmed;')
        w('   upper side (no TT32) remains UNVERIFIED (direct search OOM).')
    else:
        w(f'=> certified lower bound trans(T511) >= {best} only (plus >= 24 via D(T255)+1).')
    w('')
    w(f'=== done in {time.time()-t_start:.1f}s ===')
    OUT.close()

if __name__ == '__main__':
    main()
