#!/usr/bin/env python3
"""
beta2_rank_formula.py - Search for combinatorial formula for beta_2 vanishing

Key insight: beta_2 = 0 iff rk(d_3) = dim(Om_2) - rk(d_2)
where rk(d_2) = (n-1)(n-2)/2 - beta_1.

So beta_2 = 0 iff rk(d_3) = dim(Om_2) - (n-1)(n-2)/2 + beta_1.

This script computes dim(Om_2), rk(d_2), rk(d_3), beta_1, beta_2
and several tournament invariants to search for patterns.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(pairs):
            if bits & (1 << k):
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A, bits

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    c = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if {i,j,k} forms a 3-cycle
                s = A[i][j] + A[j][k] + A[k][i]
                if s == 3 or s == 0:  # cyclic or anti-cyclic
                    c += 1
    return c

def score_sequence(A, n):
    """Return sorted score sequence."""
    return tuple(sorted(sum(A[i]) for i in range(n)))

def transitive_triples(A, n):
    """Count transitive triples (a->b->c with a->c)."""
    tt = 0
    for a in range(n):
        for b in range(n):
            if a == b or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    tt += 1
    return tt

def compute_beta12(A, n, max_dim=3):
    """Compute beta_1 and beta_2 (and beta_3 if max_dim>=3)."""
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
            if not allowed[p]:
                break

    omega_basis = {}
    for p in range(max_dim + 2):
        if p not in allowed or not allowed[p]:
            omega_basis[p] = np.zeros((0, 0))
            continue
        basis = compute_omega_basis(A, n, p, allowed[p],
                                     allowed[p-1] if p-1 in allowed else [])
        omega_basis[p] = basis

    bd_omega = {}
    for p in range(1, max_dim + 2):
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        if dim_p == 0:
            continue
        bd = build_full_boundary_matrix(allowed[p], allowed[p-1] if p-1 in allowed else [])
        bd_omega[p] = bd @ omega_basis[p]

    betti = []
    dims = []
    ranks_down = []  # rk(d_p)
    ranks_up = []    # rk(d_{p+1})

    for p in range(max_dim + 1):
        dim_p = omega_basis[p].shape[1] if omega_basis[p].ndim == 2 else 0
        dims.append(dim_p)
        if dim_p == 0:
            betti.append(0)
            ranks_down.append(0)
            ranks_up.append(0)
            continue

        if p not in bd_omega:
            ker = dim_p
            rk_d = 0
        else:
            S_vals = np.linalg.svd(bd_omega[p], compute_uv=False)
            rk_d = sum(s > 1e-8 for s in S_vals)
            ker = dim_p - rk_d
        ranks_down.append(rk_d)

        if p+1 not in bd_omega:
            im = 0
        else:
            S_vals = np.linalg.svd(bd_omega[p+1], compute_uv=False)
            im = sum(s > 1e-8 for s in S_vals)
        ranks_up.append(im)

        betti.append(ker - im)

    return betti, dims, ranks_down, ranks_up


print("=" * 70)
print("BETA_2 RANK FORMULA SEARCH")
print("=" * 70)

for n in [4, 5, 6]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    t0 = time.time()

    # Collect data
    data = []
    count = 0

    if n <= 5:
        gen = all_tournaments(n)
    else:
        # Sample for n=6
        import random
        random.seed(42)
        pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
        m = len(pairs)
        def sample_tournaments(n, num):
            for _ in range(num):
                A = [[0]*n for _ in range(n)]
                for i,j in pairs:
                    if random.random() < 0.5:
                        A[i][j] = 1
                    else:
                        A[j][i] = 1
                yield A, 0
        gen = all_tournaments(n) if n <= 6 else sample_tournaments(n, 1000)

    results = {}  # (t3, score, tt) -> (beta, dims, rk_down, rk_up, count)

    for A, bits in gen:
        count += 1
        if count % 5000 == 0:
            print(f"  ... {count} ({time.time()-t0:.0f}s)")

        t3 = count_3cycles(A, n)
        sc = score_sequence(A, n)
        tt = transitive_triples(A, n)

        betti, dims, rk_down, rk_up = compute_beta12(A, n, max_dim=3)

        key = (t3, sc)
        if key not in results:
            results[key] = {
                'beta': betti,
                'dims': dims,
                'rk_down': rk_down,
                'rk_up': rk_up,
                'tt': tt,
                'count': 0,
                'beta2_nonzero': 0
            }
        results[key]['count'] += 1
        if len(betti) > 2 and betti[2] != 0:
            results[key]['beta2_nonzero'] += 1

    print(f"\n  Total: {count} tournaments in {time.time()-t0:.1f}s")
    print(f"\n  Key: (t3, score) -> Om dims, rk_d2, rk_d3, beta_1, beta_2")
    print(f"  {'t3':>3} {'score':<25} {'Om0':>4} {'Om1':>4} {'Om2':>4} {'Om3':>4} {'rkd2':>5} {'rkd3':>5} {'b1':>3} {'b2':>3} {'cnt':>5}")

    for key in sorted(results.keys()):
        r = results[key]
        t3, sc = key
        d = r['dims']
        rd = r['rk_down']
        ru = r['rk_up']
        b = r['beta']
        # rk_d2 = rd[2] (rank of d_2: Om_2 -> Om_1)
        # rk_d3 = ru[2] = rd[3] projected (rank of d_3: Om_3 -> Om_2, i.e. im in Om_2)
        # Actually rk_up[p] = im(d_{p+1}) landing in Om_p
        rk_d2 = rd[2] if len(rd) > 2 else 0
        rk_d3 = ru[2] if len(ru) > 2 else 0
        b1 = b[1] if len(b) > 1 else 0
        b2 = b[2] if len(b) > 2 else 0
        d0 = d[0] if len(d) > 0 else 0
        d1 = d[1] if len(d) > 1 else 0
        d2 = d[2] if len(d) > 2 else 0
        d3 = d[3] if len(d) > 3 else 0

        # Check: rk_d2 + rk_d3 should equal d2 if beta_2=0
        check = "OK" if rk_d2 + rk_d3 == d2 else "FAIL"

        print(f"  {t3:>3} {str(sc):<25} {d0:>4} {d1:>4} {d2:>4} {d3:>4} {rk_d2:>5} {rk_d3:>5} {b1:>3} {b2:>3} {r['count']:>5}  {check}")

    # Formula search
    print(f"\n  --- Formula analysis ---")
    print(f"  (n-1)(n-2)/2 = {(n-1)*(n-2)//2}")
    print(f"  C(n,3) = {n*(n-1)*(n-2)//6}")
    for key in sorted(results.keys()):
        r = results[key]
        t3, sc = key
        d2 = r['dims'][2] if len(r['dims']) > 2 else 0
        tt_val = r['tt']
        cn3 = n*(n-1)*(n-2)//6
        # dim(Om_2) = C(n,3) - t3 + cancellation_dim
        cancel = d2 - (cn3 - t3)
        print(f"    t3={t3}, score={sc}: Om_2={d2}, TT={cn3-t3}, cancel={cancel}, tt_directed={tt_val}")

print(f"\nDone.")
