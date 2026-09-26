#!/usr/bin/env python3
"""
procgen_cubedist_20260925_run.py -- full pipeline of the cube-distance lane (session collatz-procgen-20260922,
2026-09-25): Collatz's distance delta_k to the provable class (i) of the strategy cube (HYP-9138).

    python3 04-computation/experiments/procgen_cubedist_20260925_run.py > 05-knowledge/results/procgen_cubedist_20260925.out

Sections
  A  Theorem 1 (the undecided-residue construction sigma_k) and Theorem 2 (necklace lower bound): checks.
  B  exact delta_k: IHS for k <= 9; k = 10 bounds from the stored certificate of the long IHS run
     (05-knowledge/results/procgen_cubedist_20260925_k10_certificate.json.gz), re-verified from scratch: every
     stored no-good is re-derived as a genuine expanding cycle and the hitting-set MIP below the lower bound is
     re-proved infeasible; the upper-bound set is re-verified class (i).  Structure of the optima.
  C  smaller explicit sets: greedy pruning inside Bad_k, k <= 14; one-flip-per-necklace and run-based sets.
  D  controls: SHEET (3n-1) and DRIFT (5n+-1).
Every claim is a check() that raises on failure.  Engine: procgen_cubedist_20260925_engine.c, compiled into
scratch/procgen_cubedist/ (caches there are not to be committed).
"""
import os
import sys
import time
import math
import json
import gzip
import platform
import resource

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_cubedist_20260925_lib as L        # noqa: E402
import procgen_cubedist_20260925_bad as B        # noqa: E402
import procgen_cubedist_20260925_exact as E      # noqa: E402
import procgen_cubedist_20260925_prune as P      # noqa: E402
import procgen_cubedist_20260925_controls as C   # noqa: E402

A_ = math.log(1.5)
B_ = math.log(2)
RESULTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '05-knowledge', 'results'))


def log(*a):
    print(*a, flush=True)


def necklace_profile(k, flips):
    """for the flip set of an optimum: how many rotations of each expanding necklace are flipped, and how many
    flipped residues lie on non-expanding necklaces"""
    K = 1 << k
    word = {r: tuple(L.parity_word(r, k)) for r in range(K)}
    canon = {}
    for r in range(K):
        w = word[r]
        canon[r] = min(w[i:] + w[:i] for i in range(k))
    flipped = set(flips)
    per = {}
    off = 0
    for r in flipped:
        c = canon[r]
        if 3 ** sum(c) > 2 ** k:
            per[c] = per.get(c, 0) + 1
        else:
            off += 1
    exp_neck = set(canon[r] for r in range(K) if 3 ** sum(word[r]) > 2 ** k)
    hist = {}
    for c in exp_neck:
        m = per.get(c, 0)
        hist[m] = hist.get(m, 0) + 1
    return dict(sorted(hist.items())), off


def main():
    t_all = time.time()
    log("procgen_cubedist_20260925 -- Collatz's distance to provability in the strategy cube (HYP-9138)")
    log(f"python {platform.python_version()}; engine {L.ensure_engine()}")
    log("")
    log("=" * 100)
    log("A. THE UNDECIDED-RESIDUE CONSTRUCTION (Theorem 1) AND THE NECKLACE LOWER BOUND (Theorem 2)")
    log("=" * 100)
    B.main(log=log, kmax_cert=20, kmax_karp=15, kmax_lemma=13)

    log("")
    log("=" * 100)
    log("B. EXACT delta_k (implicit hitting set, HiGHS MIP; every optimum re-verified by exact Karp)")
    log("=" * 100)
    ref = {2: 1, 3: 2, 4: 2, 5: 4, 6: 5, 7: 9, 8: 14, 9: 23}      # THM-4474 section 6 (cube lane, MaxSAT)
    opt = {}
    log("   k  delta_k  N_k  |Bad_k|  Haar delta_k/2^(k-1)  rho_max  iterations  no-goods(seed)  seconds")
    for k in range(2, 10):
        r = E.exact_delta(k, verbose=False, shortP=14)
        L.check(r['status'] == 'optimal', f"IHS failed at k={k}")
        L.check(r['delta'] == ref[k], f"delta_{k} differs from the cube lane's value")
        N = L.necklace_lower_bound(k)
        Bk = L.ballot_count(k)
        L.check(N <= r['delta'] <= Bk, "sandwich N_k <= delta_k <= |Bad_k| fails")
        opt[k] = r['changed']
        log(f"  {k:2d}  {r['delta']:5d}  {N:4d}  {Bk:6d}     {r['delta'] / 2 ** (k - 1):.4f}            "
            f"{str(r['rho_max']):>5}    {r['iterations']:4d}      {r['clauses']:5d}({r['seed']})     {r['seconds']:.0f}")
        log(f"        optimum: {r['changed']}")
    log("  delta_k for k <= 9 agree with the cube lane's MaxSAT values (independent solver, independent code).")
    cert = os.path.join(RESULTS, 'procgen_cubedist_20260925_k10_certificate.json.gz')
    if os.path.exists(cert):
        t0 = time.time()
        data = json.load(gzip.open(cert, 'rt'))
        k = data['k']
        S = E.IHS(k, 3, 'plus')
        for c, m in data['nogoods']:
            L.check(E.nogood_ok(k, 3, c, frozenset(m)), f"k={k} stored no-good is not a genuine expanding cycle")
            S.add(c, frozenset(m))
        lb = data['LB']
        x, st = E.solve_hitting(len(S.res), S.clauses, cutoff=lb - 1)
        L.check(st == 'cutoff', f"k={k}: a hitting set of size {lb - 1} exists; stored lower bound not certified")
        best = data['best']
        minus = set(best)
        d = L.karp_density(k, L.mask_of(minus))
        L.check(3 ** d.numerator < 2 ** d.denominator, f"k={k} upper-bound set not class (i)")
        F = L.best_lower_approx(1 << k)
        stc, psi = L.certificate(k, L.mask_of(minus), F.numerator, F.denominator)
        L.check(stc == 'OK' and L.verify_certificate(k, minus, psi, F.numerator, F.denominator), "k=10 UB certificate")
        opt[k] = best
        tag = 'EXACT' if lb == len(best) else 'bounds'
        log(f"  {k:2d}  {lb} <= delta_{k} <= {len(best)} ({tag})  N_k = {L.necklace_lower_bound(k)}, |Bad_k| = "
            f"{L.ballot_count(k)}; Haar in [{lb / 2 ** (k - 1):.4f}, {len(best) / 2 ** (k - 1):.4f}]")
        log(f"        lower bound: {len(S.clauses) - S.nseed} stored no-goods + {S.nseed} seeds, each re-derived as a genuine "
            f"expanding cycle; the hitting-set MIP with sum x <= {lb - 1} is infeasible  [{time.time() - t0:.0f}s]")
        log(f"        upper bound: a class-(i) set of size {len(best)} (Karp rho_max {d}; edge-checked certificate): {best}")
    else:
        log("  k=10: certificate file not present -- skipped")
    log(f"  11  {L.necklace_lower_bound(11)} = N_11 <= delta_11 <= 72 (the pruned set of section C)")
    log("")
    log("  structure of the optima (k = 10: the best known set, not proved optimal; u->d: residue 3 mod 4 flipped;")
    log("  d->u: 1 mod 4 flipped; necklace profile =")
    log("  number of expanding k-necklaces with m flipped rotations, {m: count}; 'off' = flips on non-expanding necklaces)")
    for k in sorted(opt):
        fl = opt[k]
        bad = set(int(r) for r in B.bad_numpy(k)[0])
        hist, off = necklace_profile(k, fl)
        n3 = sum(1 for r in fl if r % 4 == 3)
        log(f"  k={k:2d}: {len(fl)} flips; in Bad_k {sum(1 for r in fl if r in bad)}; u->d {n3}, d->u {len(fl) - n3}; "
            f"-1 flipped: {(1 << k) - 1 in fl}; -5 flipped: {(1 << k) - 5 in fl}; necklace profile {hist}, off {off}")
        L.check(0 not in hist, "an expanding necklace without a flip (contradicts Theorem 2's argument)")

    log("")
    log("=" * 100)
    log("C. SMALLER EXPLICIT SETS INSIDE Bad_k (greedy pruning; upper bounds on delta_k)")
    log("=" * 100)
    P.main(log=log, kmax=14)

    log("")
    log("=" * 100)
    log("D. CONTROLS")
    log("=" * 100)
    C.main(log=log, sheet_kmax=9, drift_kmax=7, collatz_opt={k: len(v) for k, v in opt.items() if k <= 9},
           drift_prune_kmax=14)

    log("")
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024 if sys.platform == 'darwin' else 1024)
    log(f"ALL CHECKS PASSED.  wall time {time.time() - t_all:.0f} s; peak RSS of this process {peak:.0f} MB")


if __name__ == '__main__':
    main()
