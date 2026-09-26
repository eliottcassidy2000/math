#!/usr/bin/env python3
"""
procgen_cubedist_20260925_controls.py -- controls for the cube-distance lane (2026-09-25, HYP-9138).

SHEET.  nu(x) = -x conjugates T_sigma to T_(nu sigma), (nu sigma)(m) = -sigma(-m) (THM-4474 Theorem E).
        It maps Collatz with flips R to 3n-1 with flips -R, so the distances of Collatz and 3n-1 to class (i)
        are equal.  Checked independently: the IHS is re-run from the all-minus strategy (no seeds), and the
        nu-images of the optimal sets and of the undecided construction are verified.
DRIFT.  The 5n+-1 cube: T(n) = n/2, (5n + sigma(n mod 2^k))/2; a cycle is expanding iff 5^a > 2^p, i.e.
        odd density > log_5 2 = 0.4307.  Class (i) is decided exactly as in Theorem A.  We decide for each
        small k whether class (i) is empty (an infeasible no-good set is an UNSAT certificate, re-checked by an
        independent SAT solver and, for k <= 5, by exhaustive enumeration), and compute the distance of
        5n+1 (all +) to class (i) when it is nonempty.
"""
import sys
import os
import time
import itertools
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_cubedist_20260925_lib as L    # noqa: E402
import procgen_cubedist_20260925_exact as E  # noqa: E402


def neg_res(r, k):
    return (-r) % (1 << k)


def nu_changes(changes, k):
    """changes relative to Collatz (all +) -> changes relative to 3n-1 (all -): the set -R"""
    return sorted(neg_res(r, k) for r in changes)


def undecided_set(k, mul, sign):
    """residues r (odd) on which the all-`sign` map (m n + sign)/2 has no descent within k steps"""
    out = []
    for r in range(1, 1 << k, 2):
        n, a, ok = r, 0, True
        for j in range(1, k + 1):
            if n & 1:
                a += 1
                n = (mul * n + sign) // 2
            else:
                n //= 2
            if not mul ** a > 2 ** j:
                ok = False
                break
        if ok:
            out.append(r)
    return out


def karp_small(k, sigma_minus, mul):
    """pure-Python exact max odd density (Karp) for small k: sigma_minus = set of residues with sigma = -1."""
    K = 1 << k
    H = K >> 1
    succ = []
    for s in range(K):
        t = s // 2 if s % 2 == 0 else (mul * s + (-1 if s in sigma_minus else 1)) // 2
        t0 = t % H
        succ.append((t0, t0 + H))
    NEG = -10 ** 9
    D = [[0] * K]
    for j in range(1, K + 1):
        prev = D[-1]
        cur = [NEG] * K
        for u in range(K):
            if prev[u] > NEG:
                val = prev[u] + (u & 1)
                for v in succ[u]:
                    if val > cur[v]:
                        cur[v] = val
        D.append(cur)
    best = None
    for v in range(K):
        if D[K][v] <= NEG:
            continue
        worst = None
        for j in range(K):
            if D[j][v] > NEG:
                f = Fraction(D[K][v] - D[j][v], K - j)
                if worst is None or f < worst:
                    worst = f
        if best is None or worst > best:
            best = worst
    return best


def sat_unsat(nogoods):
    """independent UNSAT check of a no-good set with a CDCL SAT solver (pysat Glucose4); variables y_v = [sigma(v) = -];
    the no-good (cyc, minus_odd) says: some odd v of cyc has y_v different from the recorded value."""
    from pysat.solvers import Glucose4
    idx = {}
    cls = []
    for cyc, mo in nogoods:
        c = []
        for v in sorted(set(u for u in cyc if u & 1)):
            j = idx.setdefault(v, len(idx) + 1)
            c.append(-j if v in mo else j)
        cls.append(c)
    with Glucose4() as s:
        for c in cls:
            s.add_clause(c)
        return not s.solve()


def stationary_min_flip_mass(k, minus, mul):
    """for every closed class of G_sigma: the stationary law pi of the uniform-lift chain; returns
    min over closed classes of (pi(flipped), pi(odd)) -- flipped = residues where sigma = - (differs from all +)."""
    import numpy as np
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    K = 1 << k
    H = K >> 1
    rows, cols = [], []
    for s in range(K):
        t = s // 2 if s % 2 == 0 else (mul * s + (-1 if s in minus else 1)) // 2
        t0 = t % H
        rows += [s, s]
        cols += [t0, t0 + H]
    A = csr_matrix((np.ones(2 * K), (rows, cols)), shape=(K, K))
    ncomp, lab = connected_components(A, directed=True, connection='strong')
    # closed classes: no edge leaves the component
    closed = []
    for c in range(ncomp):
        nodes = np.nonzero(lab == c)[0]
        ok = True
        for s in nodes:
            for t in A.indices[A.indptr[s]:A.indptr[s + 1]]:
                if lab[t] != c:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            closed.append(nodes)
    out = []
    for nodes in closed:
        pos = {int(s): i for i, s in enumerate(nodes)}
        n = len(nodes)
        pi = np.ones(n) / n
        P = np.zeros((n, n)) if n <= 4096 else None
        L.check(P is not None, "closed class too large for the dense solve")
        for s in nodes:
            t = s // 2 if s % 2 == 0 else (mul * s + (-1 if s in minus else 1)) // 2
            t0 = t % H
            for u in (t0, t0 + H):
                P[pos[int(s)], pos[int(u)]] += 0.5
        # solve pi (P - I) = 0, sum pi = 1
        Mt = (P - np.eye(n)).T
        Mt[-1, :] = 1.0
        b = np.zeros(n)
        b[-1] = 1.0
        pi = np.linalg.solve(Mt, b)
        L.check(np.all(pi > -1e-12), "negative stationary mass")
        pf = sum(pi[pos[int(s)]] for s in nodes if s in minus)
        po = sum(pi[pos[int(s)]] for s in nodes if s & 1)
        out.append((n, pf, po))
    return out


def main(log=print, sheet_kmax=9, drift_kmax=8, collatz_opt=None, drift_prune_kmax=14):
    log("S. SHEET control: Collatz vs 3n-1 (nu = negation)")
    for k in range(2, sheet_kmax + 1):
        rc = E.exact_delta(k, mul=3, base='plus', verbose=False, shortP=14)
        rm = E.exact_delta(k, mul=3, base='minus', verbose=False)
        L.check(rc['status'] == 'optimal' and rm['status'] == 'optimal', "IHS failed")
        L.check(rc['delta'] == rm['delta'], f"sheet distances differ at k={k}")
        if collatz_opt is not None and k in collatz_opt:
            L.check(rc['delta'] == collatz_opt[k], "Collatz delta differs from the reference")
        # nu-image of the Collatz optimum is a class-(i) 3n-1 neighbour at the same distance
        img = nu_changes(rc['changed'], k)
        minus_set = set(range(1, 1 << k, 2)) - set(img)       # sigma = - except at the changed residues
        d = L.karp_density(k, L.mask_of(minus_set))
        L.check(3 ** d.numerator < 2 ** d.denominator, "nu-image not class (i)")
        L.check(d == rc['rho_max'], "nu does not preserve rho_max")
        # nu-image of the undecided construction = the 3n-1 undecided construction
        badc = undecided_set(k, 3, +1)
        badm = undecided_set(k, 3, -1)
        L.check(sorted(neg_res(r, k) for r in badc) == badm, "-Bad_k(3n+1) != Bad_k(3n-1)")
        ms = set(range(1, 1 << k, 2)) - set(badm)
        d2 = L.karp_density(k, L.mask_of(ms))
        L.check(3 ** d2.numerator < 2 ** d2.denominator, "3n-1 undecided construction not class (i)")
        log(f"  k={k}: delta(Collatz) = delta(3n-1) = {rc['delta']} (independent IHS runs, {rm['iterations']} "
            f"iterations from all-minus); nu(optimum) class (i) with rho_max {d}; -Bad_k(3n+1) = Bad_k(3n-1) "
            f"({len(badm)} residues), whose flip is class (i) (rho_max {d2})")

    log("")
    log("D. DRIFT control: the 5n+-1 cube (expanding iff 5^a > 2^p; threshold log_5 2 = 0.430677)")
    res = {}
    for k in range(2, drift_kmax + 1):
        t0 = time.time()
        r = E.exact_delta(k, mul=5, base='plus', P=k + 2, verbose=False)
        res[k] = r
        if r['status'] == 'infeasible':
            S = r['ihs']
            for cyc, mo in S.nogoods:
                L.check(E.nogood_ok(k, 5, cyc, mo), "stored 5n+-1 no-good is not a genuine expanding cycle")
            L.check(sat_unsat(S.nogoods), f"SAT solver finds the 5n+-1 no-good set satisfiable at k={k}")
            exhaust = ''
            if k <= 5:
                cnt = 0
                odd = list(range(1, 1 << k, 2))
                for bits in itertools.product((0, 1), repeat=len(odd)):
                    minus = set(o for o, b in zip(odd, bits) if b)
                    d = karp_small(k, minus, 5)
                    if 5 ** d.numerator < 2 ** d.denominator:
                        cnt += 1
                L.check(cnt == 0, f"exhaustive search finds a 5n+-1 class-(i) strategy at k={k}")
                exhaust = f"; exhaustive: 0 of {2 ** len(odd)} strategies in class (i)"
            log(f"  k={k}: class (i) EMPTY -- {r['clauses']} expanding-cycle no-goods admit no hitting set "
                f"(MIP infeasible; UNSAT re-confirmed by Glucose4){exhaust}  [{time.time() - t0:.1f}s]")
        elif r['status'] == 'optimal':
            dd = r['rho_max']
            log(f"  k={k}: class (i) nonempty; distance of 5n+1 = {r['delta']} of {1 << (k - 1)} residues "
                f"(Haar {r['delta'] / (1 << (k - 1)):.4f}); optimum rho_max {dd} = {float(dd):.4f} "
                f"[{r['iterations']} iterations, {r['clauses']} no-goods, {time.time() - t0:.0f}s]")
            log(f"        changed residues: {r['changed']}")
        else:
            log(f"  k={k}: {r['status']} after {r.get('iterations')} iterations")

    # class (i) is closed under lifting (the lift has the same map), so it is nonempty at every level >= the first
    first = min(k for k, r in res.items() if r['status'] == 'optimal')
    L.check(all(res[k]['status'] == 'infeasible' for k in res if k < first), "non-monotone emptiness")
    log(f"  => the 5n+-1 provable class is empty at levels 2..{first - 1} and nonempty at every level >= {first} "
        f"(lifting preserves class (i)).")

    log("")
    log("D2. 5n+1: lower bound (necklaces, the proof of Theorem 2 verbatim), heuristic upper bound (greedy pruning")
    log("    of the lifted level-%d optimum, every set verified class (i)), and the stationary flip-mass lemma" % first)
    log("   k   N_k(5) lower  Haar     pruned upper  Haar     k*Haar(upper)  closed classes: min pi(flipped) >= (1/2-log_5 2)/k ?")
    import random
    cur = set(res[first]['changed'])
    thr_c = 0.5 - __import__('math').log(2) / __import__('math').log(5)
    for k in range(first, drift_prune_kmax + 1):
        t0 = time.time()
        if k > first:
            cur = set(r for r in range(1, 1 << k, 2) if (r % (1 << (k - 1))) in cur)
        thr = L.expanding_threshold(k, 5)
        st, _ = L.certificate(k, L.mask_of(cur), thr.numerator, thr.denominator, mul=5)
        L.check(st == 'OK', "lifted 5n+-1 set not class (i)")
        order = sorted(cur)
        random.Random(k).shuffle(order)
        for r in order:
            cur.discard(r)
            st, _ = L.certificate(k, L.mask_of(cur), thr.numerator, thr.denominator, mul=5)
            if st != 'OK':
                cur.add(r)
        F = L.best_lower_approx(1 << k, 5)
        st, psi = L.certificate(k, L.mask_of(cur), F.numerator, F.denominator, mul=5)
        L.check(st == 'OK' and L.verify_certificate(k, cur, psi, F.numerator, F.denominator, mul=5),
                "pruned 5n+-1 set: certificate fails")
        N5 = L.necklace_lower_bound(k, 5)
        stat = stationary_min_flip_mass(k, cur, 5) if k <= 11 else None
        stxt = '-'
        if stat is not None:
            mn = min(pf for _, pf, _ in stat)
            L.check(mn >= thr_c / k - 1e-12, "stationary flip-mass lemma violated")
            L.check(all(po < __import__('math').log(2) / __import__('math').log(5) for _, _, po in stat), "sandwich")
            stxt = f"{len(stat)} class(es), min pi(flipped) = {mn:.4f} >= {thr_c / k:.4f}: yes"
        log(f"  {k:2d}   {N5:8d}   {N5 / 2 ** (k - 1):.4f}   {len(cur):8d}    {len(cur) / 2 ** (k - 1):.4f}   "
            f"{k * len(cur) / 2 ** (k - 1):6.3f}        {stxt}  [{time.time() - t0:.0f}s]")
    return res


if __name__ == '__main__':
    main()