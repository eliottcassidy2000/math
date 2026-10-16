"""Referee check (D): witness rules (HYP-9280 post-audit 'Which in-neighbour is the witness?'; results section 11.4).

Independent of the session's subset DP: explicit DFS over Hamiltonian paths (pa_core.enumerate_paths_dfs, which also
returns L[w, u] = # paths with u immediately before w) and brute force over all n! permutations for the explicit examples.
Rules at a vertex w with in-neighbourhood U (as defined in the session script):
  W  predecessor-weighted average  sum_u Pr[pred(w) = u] (2 p(u, w) - 1) = sum_u L[w,u] (2 C[u,w] - H) / H^2  (> 0 ?)
  A  in-neighbour of minimum out-degree (in T);   D  minimum score inside T[U];   E  maximum score inside T[U];
  G  beats the fewest out-neighbours of w;          P  the most frequent immediate predecessor of w (max L[w, u]).
'some-fail' = NO tied choice is a witness (p > 1/2); 'all-fail' = some tied choice is not a witness.
Modes:  w10       the n = 10 tournament 'transitive on 0..9, higher beats lower, arcs 0 -> 4 and 1 -> 9': exact W at
                  every vertex (integers / Fractions), DFS and brute force over 10! permutations.
        lab NMAX  all labelled strong tournaments n = 4..NMAX: counts for every rule (to compare with the session) and
                  the first explicit examples (smallest code) of A-some-fail, P-some-fail, D-some-fail, E-some-fail,
                  G-all-fail; each example re-verified by brute force over n! permutations with exact fractions.
        near n k ...  all strong 'transitive + k reversed arcs (hi - lo >= 2)': counts for every rule; first example of
                  G-some-fail and of W <= 0 (re-verified by brute force).
"""
import sys, time, itertools
from fractions import Fraction as Fr
import numpy as np
import numba as nb
from pa_core import enumerate_paths_dfs, brute_force_perms, strong_landau, transitive_plus

NRULE = 5
NAMES = ['A min out-degree', 'D min score in T[U]', 'E max score in T[U]', 'G beats fewest out-nbrs of w',
         'P most frequent predecessor']


@nb.njit(cache=False)
def rule_key(r, u, n, A, w, L):
    if r == 0:
        s = 0
        for x in range(n):
            s += A[u, x]
        return s
    if r == 1 or r == 2:
        s = 0
        for x in range(n):
            if x != w and A[x, w] and A[u, x]:
                s += 1
        return s if r == 1 else -s
    if r == 3:
        s = 0
        for x in range(n):
            if A[w, x] and A[u, x]:
                s += 1
        return s
    return -L[w, u]


@nb.njit(cache=False)
def eval_rules(n, A, C, L, H, res, ex_some, ex_all, code):
    """res[r] some-fail, res[NRULE + r] all-fail, res[2 NRULE] W <= 0, res[2 NRULE + 1] pairs, res[2 NRULE + 2] W == 0.
    ex_some[r] / ex_all[r] / ex_some[NRULE] (W): first (code, w) seen, stored as code * 16 + w (kept minimal)."""
    for w in range(n):
        res[2 * NRULE + 1] += 1
        Wsum = 0
        for u in range(n):
            if u != w and A[u, w]:
                Wsum += L[w, u] * (2 * C[u, w] - H)
        if Wsum <= 0:
            res[2 * NRULE] += 1
            key = code * 16 + w
            if ex_some[NRULE] < 0 or key < ex_some[NRULE]:
                ex_some[NRULE] = key
        if Wsum == 0:
            res[2 * NRULE + 2] += 1
        for r in range(NRULE):
            best = 1 << 62
            for u in range(n):
                if u != w and A[u, w]:
                    k = rule_key(r, u, n, A, w, L)
                    if k < best:
                        best = k
            some = False
            allw = True
            for u in range(n):
                if u != w and A[u, w] and rule_key(r, u, n, A, w, L) == best:
                    if 2 * C[u, w] > H:
                        some = True
                    else:
                        allw = False
            if not some:
                res[r] += 1
                key = code * 16 + w
                if ex_some[r] < 0 or key < ex_some[r]:
                    ex_some[r] = key
            if not allw:
                res[NRULE + r] += 1
                key = code * 16 + w
                if ex_all[r] < 0 or key < ex_all[r]:
                    ex_all[r] = key


@nb.njit(parallel=True, cache=False)
def sweep_lab(n, pu, pv, nchunk):
    E = pu.shape[0]
    total = 1 << E
    res = np.zeros((nchunk, 2 * NRULE + 3), dtype=np.int64)
    exs = np.full((nchunk, NRULE + 1), -1, dtype=np.int64)
    exa = np.full((nchunk, NRULE), -1, dtype=np.int64)
    nstrong = np.zeros(nchunk, dtype=np.int64)
    chunk = (total + nchunk - 1) // nchunk
    for c in nb.prange(nchunk):
        A = np.zeros((n, n), dtype=np.int8)
        C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
        for code in range(c * chunk, min(total, (c + 1) * chunk)):
            for a in range(n):
                for b in range(n):
                    A[a, b] = 0
            for e in range(E):
                if (code >> e) & 1:
                    A[pu[e], pv[e]] = 1
                else:
                    A[pv[e], pu[e]] = 1
            if not strong_landau(n, A):
                continue
            nstrong[c] += 1
            H = enumerate_paths_dfs(n, A, C, Nd, L)
            eval_rules(n, A, C, L, H, res[c], exs[c], exa[c], code)
    tot = np.zeros(2 * NRULE + 3, dtype=np.int64)
    for c in range(nchunk):
        for i in range(2 * NRULE + 3):
            tot[i] += res[c, i]
    ex_s = np.full(NRULE + 1, -1, dtype=np.int64)
    ex_a = np.full(NRULE, -1, dtype=np.int64)
    for c in range(nchunk):
        for r in range(NRULE + 1):
            if exs[c, r] >= 0 and (ex_s[r] < 0 or exs[c, r] < ex_s[r]):
                ex_s[r] = exs[c, r]
        for r in range(NRULE):
            if exa[c, r] >= 0 and (ex_a[r] < 0 or exa[c, r] < ex_a[r]):
                ex_a[r] = exa[c, r]
    return tot, ex_s, ex_a, nstrong.sum()


@nb.njit(parallel=True, cache=False)
def sweep_near(n, arcsets, nthreads):
    m = arcsets.shape[0]; k = arcsets.shape[1]
    res = np.zeros((nthreads, 2 * NRULE + 3), dtype=np.int64)
    exs = np.full((nthreads, NRULE + 1), -1, dtype=np.int64)
    exa = np.full((nthreads, NRULE), -1, dtype=np.int64)
    chunk = (m + nthreads - 1) // nthreads
    for t in nb.prange(nthreads):
        A = np.zeros((n, n), dtype=np.int8)
        C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
        for idx in range(t * chunk, min(m, (t + 1) * chunk)):
            for a in range(n):
                for b in range(n):
                    A[a, b] = 1 if b < a else 0
            for e in range(k):
                lo = arcsets[idx, e, 0]; hi = arcsets[idx, e, 1]
                A[hi, lo] = 0; A[lo, hi] = 1
            if not strong_landau(n, A):
                continue
            H = enumerate_paths_dfs(n, A, C, Nd, L)
            eval_rules(n, A, C, L, H, res[t], exs[t], exa[t], idx)
    tot = np.zeros(2 * NRULE + 3, dtype=np.int64)
    for c in range(nthreads):
        for i in range(2 * NRULE + 3):
            tot[i] += res[c, i]
    ex_s = np.full(NRULE + 1, -1, dtype=np.int64)
    ex_a = np.full(NRULE, -1, dtype=np.int64)
    for c in range(nthreads):
        for r in range(NRULE + 1):
            if exs[c, r] >= 0 and (ex_s[r] < 0 or exs[c, r] < ex_s[r]):
                ex_s[r] = exs[c, r]
        for r in range(NRULE):
            if exa[c, r] >= 0 and (ex_a[r] < 0 or exa[c, r] < ex_a[r]):
                ex_a[r] = exa[c, r]
    return tot, ex_s, ex_a


def report(tag, r, t0):
    print(f'{tag}: (T, w) pairs {int(r[2 * NRULE + 1])}; predecessor-weighted average <= 0: {int(r[2 * NRULE])} (exactly 0: {int(r[2 * NRULE + 2])}); ' +
          '; '.join(f'{NAMES[i]}: some-fail {int(r[i])}, all-fail {int(r[NRULE + i])}' for i in range(NRULE)) +
          f'  [{time.time() - t0:.1f} s]', flush=True)


def outsets(n, A):
    return [sorted(int(x) for x in range(n) if A[v, x]) for v in range(n)]


def describe(n, A, w, tag):
    """Exact re-verification of one (T, w) by brute force over all n! permutations; prints the data."""
    H, C, Nd, L = brute_force_perms(n, A)
    H = int(H)
    ins = [u for u in range(n) if u != w and A[u, w]]
    outs = [x for x in range(n) if A[w, x]]
    print(f'  [{tag}] n = {n}, w = {w}; out-neighbour sets (v: [x with v -> x]): {outsets(n, A)}', flush=True)
    print(f'     H = {H} (brute force over {n}! permutations); in-neighbours of w: {ins}; out-neighbours of w: {outs}')
    Wsum = 0
    for u in ins:
        p = Fr(int(C[u, w]), H)
        outdeg = int(A[u].sum())
        scU = sum(1 for x in ins if A[u, x])
        g = sum(1 for x in outs if A[u, x])
        Wsum += int(L[w, u]) * (2 * int(C[u, w]) - H)
        print(f'     u = {u}: p(u before w) = {p} ({"witness" if p > Fr(1, 2) else "NOT a witness"}); out-degree {outdeg}; '
              f'score in T[U] {scU}; beats {g} out-nbrs of w; Pr[pred(w) = u] = {Fr(int(L[w, u]), H)}')
    print(f'     predecessor-weighted average sum_u Pr[pred = u](2p - 1) = {Fr(Wsum, H * H)}; '
          f'best in-neighbour p = {max(Fr(int(C[u, w]), H) for u in ins)}', flush=True)
    return H, C, L


def decode_lab(n, code):
    prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    A = np.zeros((n, n), dtype=np.int8)
    for e, (i, j) in enumerate(prs):
        if (code >> e) & 1:
            A[i, j] = 1
        else:
            A[j, i] = 1
    return A


def mode_w10():
    n = 10
    A = transitive_plus(n, [(0, 4), (1, 9)])
    print('MODE w10: transitive on 0..9 (higher beats lower) with 0 -> 4 and 1 -> 9', flush=True)
    print(f'  strong (Landau): {strong_landau(n, A)}; out-neighbour sets: {outsets(n, A)}')
    C = np.zeros((n, n), dtype=np.int64); Nd = np.zeros((n, n + 1), dtype=np.int64); L = np.zeros((n, n), dtype=np.int64)
    H = int(enumerate_paths_dfs(n, A, C, Nd, L))
    t0 = time.time()
    H2, C2, Nd2, L2 = brute_force_perms(n, A)
    same = int(H2) == H and (C2 == C).all() and (L2 == L).all()
    print(f'  H = {H}; DFS and brute force over 10! permutations agree on H, C, L: {same}  [{time.time() - t0:.1f} s]')
    for w in range(n):
        ins = [u for u in range(n) if u != w and A[u, w]]
        Wsum = sum(int(L[w, u]) * (2 * int(C[u, w]) - H) for u in ins)
        best = max(Fr(int(C[u, w]), H) for u in ins)
        print(f'  w = {w}: in-degree {len(ins)}; W(w) = {Fr(Wsum, H * H)} = {Wsum / H / H:+.6f}  ({"<= 0" if Wsum <= 0 else "> 0"}); '
              f'best in-neighbour p = {best}')
    for w in range(n):
        ins = [u for u in range(n) if u != w and A[u, w]]
        Wsum = sum(int(L[w, u]) * (2 * int(C[u, w]) - H) for u in ins)
        if Wsum <= 0:
            describe(n, A, w, 'W <= 0')


def mode_lab(nmax):
    print('MODE lab: all labelled strong tournaments', flush=True)
    for n in range(4, nmax + 1):
        t0 = time.time()
        prs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pu = np.array([p[0] for p in prs], dtype=np.int64); pv = np.array([p[1] for p in prs], dtype=np.int64)
        r, ex_s, ex_a, ns = sweep_lab(n, pu, pv, 512)
        report(f'n = {n} labelled ({int(ns)} strong)', r, t0)
        if n >= 6:
            for rr, tag in ((0, 'A some-fail'), (4, 'P some-fail'), (1, 'D some-fail'), (2, 'E some-fail')):
                if ex_s[rr] >= 0:
                    code, w = int(ex_s[rr]) // 16, int(ex_s[rr]) % 16
                    describe(n, decode_lab(n, code), w, f'{tag}, labelled code {code}')
            if ex_a[3] >= 0:
                code, w = int(ex_a[3]) // 16, int(ex_a[3]) % 16
                describe(n, decode_lab(n, code), w, f'G all-fail (a tied G-choice is not a witness), labelled code {code}')


def mode_near(jobs):
    print('MODE near: all strong transitive + k reversed arcs (hi - lo >= 2)', flush=True)
    for n, k in jobs:
        t0 = time.time()
        cand = [(lo, hi) for lo in range(n) for hi in range(lo + 2, n)]
        combos = np.array(list(itertools.combinations(cand, k)), dtype=np.int64).reshape(-1, k, 2)
        r, ex_s, ex_a = sweep_near(n, combos, 12)
        report(f'n = {n}, transitive + {k} reversed arcs', r, t0)
        for rr, tag in ((3, 'G some-fail'), (NRULE, 'W <= 0')):
            if ex_s[rr] >= 0 and n <= 10:
                idx, w = int(ex_s[rr]) // 16, int(ex_s[rr]) % 16
                arcs = [tuple(map(int, a)) for a in combos[idx]]
                print(f'  first {tag} example: reversed arcs {arcs}, w = {w}')
                describe(n, transitive_plus(n, arcs), w, tag)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'w10':
        mode_w10()
    elif mode == 'lab':
        mode_lab(int(sys.argv[2]))
    elif mode == 'near':
        a = [int(x) for x in sys.argv[2:]]
        mode_near(list(zip(a[0::2], a[1::2])))
