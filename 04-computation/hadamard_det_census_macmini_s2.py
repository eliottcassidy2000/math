#!/usr/bin/env python3
"""
hadamard_det_census_macmini_s2.py — mac-mini-2026-06-10-S2

THE DETERMINANT LENS ON TOURNAMENTS (THM-466 / HYP-2383..2388 / T777).

For each tournament iso class at n=3..8 (via nauty gentourng for n>=5),
compute the tournament determinant det(I+S), S = A - A^T, and test:

  1. FLOOR (HYP-2384):  det(I+S) = 2^(n-1)  <=>  locally transitive
                        <=>  switching of a transitive tournament.
     Equality counts should be (2n-2)!! labeled.
  2. CEILING (HYP-2385): odd n: det(I+S) <= (n+1)^((n-1)/2), equality iff
     switching class of a DRT.  even n: det <= n^(n/2), equality iff
     SS^T = (n-1)I (skew conference).
  3. AVERAGE (HYP-2383): E_labeled[det(I+S)] = involution numbers I(n).
     Exact via class weights n!/|Aut| at n<=7; Monte Carlo at n=8.
  4. HERMITE (HYP-2387): E[det(xI+S)] = matching polynomial of K_n
     (exact labeled at n<=6).
  5. ALIGNMENT (HYP-2386): correlation of d(T) with H(T) per class;
     argmax_d vs argmax_H; complement-invariance of d; |Pf(S)| cross-links
     (THM-174/442 world) at even n.

Provenance: mac-mini-2026-06-10-S2. Reuses gentourng parsing pattern from
h63_n8_isoclass_census_s10.py (opus); H bitmask DP standard.
"""

import sys, os, subprocess, itertools, time
from math import comb, factorial
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------- utilities

def pairs_of(n):
    return list(itertools.combinations(range(n), 2))

def adj_from_gentourng_bits(bits, n, m, PAIRS):
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(PAIRS):
        if bits & (1 << (m - 1 - k)):   # nauty: MSB first
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A

def read_gentourng(n):
    m = n*(n-1)//2
    result = subprocess.run(['gentourng', str(n)], capture_output=True, text=True)
    text = result.stdout if result.stdout else result.stderr
    reps = []
    for line in text.splitlines():
        line = line.strip()
        if len(line) == m and all(c in '01' for c in line):
            reps.append(int(line, 2))
    PAIRS = pairs_of(n)
    return [adj_from_gentourng_bits(b, n, m, PAIRS) for b in reps]

def all_tournaments_small(n):
    """All labeled tournaments (n<=4 for class gen; n<=6 for labeled sweeps)."""
    PAIRS = pairs_of(n)
    for bits in range(2**len(PAIRS)):
        A = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(PAIRS):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def canon(A, n):
    """Lex-min canonical form (upper-triangle bitstring) over all permutations."""
    best = None
    for p in itertools.permutations(range(n)):
        bits = 0
        k = 0
        for i in range(n):
            for j in range(i+1, n):
                bits = (bits << 1) | A[p[i]][p[j]]
        if best is None or bits < best:
            best = bits
    return best

def aut_count(A, n):
    cnt = 0
    for p in itertools.permutations(range(n)):
        if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n) if i != j):
            cnt += 1
    return cnt

def int_det(Mrows):
    """Exact integer determinant, fraction-free Bareiss. Mrows: list of lists of ints."""
    M = [row[:] for row in Mrows]
    n = len(M)
    sign = 1
    prev = 1
    for k in range(n-1):
        if M[k][k] == 0:
            for r in range(k+1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n-1][n-1]

def S_of(A, n):
    return [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

def detIS(A, n):
    S = S_of(A, n)
    M = [[(1 if i == j else 0) + S[i][j] for j in range(n)] for i in range(n)]
    return int_det(M)

def isqrt_exact(v):
    r = int(round(v ** 0.5))
    for c in (r-1, r, r+1):
        if c >= 0 and c*c == v: return c
    raise ValueError(f'not a perfect square: {v}')

def H_count(A, n):
    """Number of Hamiltonian paths, bitmask DP."""
    out = [0]*n
    for i in range(n):
        for j in range(n):
            if A[i][j]: out[i] |= (1 << j)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c:
                nxt = out[last] & ~mask
                while nxt:
                    b = nxt & (-nxt)
                    j = b.bit_length() - 1
                    dp[mask | b][j] += c
                    nxt ^= b
    full = (1 << n) - 1
    return sum(dp[full])

def three_cycles(A, n):
    """List of frozensets {a,b,c} that induce a 3-cycle."""
    out = []
    for a, b, c in itertools.combinations(range(n), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[b][a] and A[c][b] and A[a][c]):
            out.append((a, b, c))
    return out

def is_locally_transitive(A, n):
    for (a, b, c) in three_cycles(A, n):
        for v in range(n):
            if v in (a, b, c): continue
            if A[v][a] and A[v][b] and A[v][c]: return False   # out-vortex
            if A[a][v] and A[b][v] and A[c][v]: return False   # in-vortex
    return True

def is_acyclic(A, n):
    return sorted(sum(row) for row in A) == list(range(n))

def switch(A, n, mask):
    """Reverse all arcs between mask-set and complement. Returns new adjacency."""
    inset = [(mask >> v) & 1 for v in range(n)]
    B = [row[:] for row in A]
    for i in range(n):
        for j in range(n):
            if i != j and inset[i] != inset[j]:
                B[i][j] = A[j][i]
    return B

def is_switch_of_transitive(A, n):
    for mask in range(1 << (n-1)):       # wlog vertex n-1 not switched
        if is_acyclic(switch(A, n, mask), n):
            return True
    return False

def SSt(A, n):
    S = S_of(A, n)
    return [[sum(S[i][k]*S[j][k] for k in range(n)) for j in range(n)] for i in range(n)]

def is_DRT(A, n):
    """SS^T = nI - J  (doubly regular tournament characterization)."""
    P = SSt(A, n)
    return all(P[i][j] == (n-1 if i == j else -1) for i in range(n) for j in range(n))

def is_skew_conference(A, n):
    P = SSt(A, n)
    return all(P[i][j] == ((n-1) if i == j else 0) for i in range(n) for j in range(n))

def is_switch_of_DRT(A, n):
    for mask in range(1 << (n-1)):
        if is_DRT(switch(A, n, mask), n):
            return True
    return False

def complement(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def involutions(n):
    a, b = 1, 1   # I(0), I(1)
    if n == 0: return 1
    for k in range(2, n+1):
        a, b = b, b + (k-1)*a
    return b

def matching_poly_Kn(n):
    """Coefficients of E[det(xI+S)] = sum_k C(n,2k)(2k-1)!! x^(n-2k), degree n -> 0."""
    coeffs = [0]*(n+1)
    for k in range(0, n//2 + 1):
        df = 1
        for t in range(2*k-1, 0, -2): df *= t
        coeffs[2*k] = comb(n, 2*k) * df     # coefficient of x^(n-2k)
    return coeffs

def char_poly_int(M, n):
    """Exact integer char poly det(xI - M) via Faddeev-LeVerrier. Returns coeffs c[0..n], c[0]=1 leading."""
    from fractions import Fraction
    Mf = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    I = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    c = [Fraction(1)]
    Mk = [row[:] for row in Mf]
    for k in range(1, n+1):
        ck = -sum(Mk[i][i] for i in range(n)) / k
        c.append(ck)
        if k < n:
            Mk = [[Mk[i][j] + (ck if i == j else 0) for j in range(n)] for i in range(n)]
            Mk = [[sum(Mf[i][t]*Mk[t][j] for t in range(n)) for j in range(n)] for i in range(n)]
    return [int(x) for x in c]

# ---------------------------------------------------------------- main census

def census(n, do_aut=True, do_hermite=False):
    print(f'\n{"="*72}\n  n = {n}\n{"="*72}')
    t0 = time.time()
    if n <= 4:
        seen = {}
        for A in all_tournaments_small(n):
            c = canon(A, n)
            if c not in seen: seen[c] = A
        classes = list(seen.values())
    else:
        classes = read_gentourng(n)
    print(f'  iso classes: {len(classes)}   ({time.time()-t0:.1f}s)')

    rows = []
    for A in classes:
        det = detIS(A, n)
        assert det % (1 << (n-1)) == 0, f'2^(n-1) divisibility FAILS: {det}'
        d = det >> (n-1)
        H = H_count(A, n)
        lt = is_locally_transitive(A, n)
        scores = tuple(sorted(sum(row) for row in A))
        aut = aut_count(A, n) if do_aut else None
        pf = None
        if n % 2 == 0:
            dS = int_det(S_of(A, n))
            pf = isqrt_exact(dS)
        rows.append(dict(A=A, det=det, d=d, H=H, lt=lt, scores=scores, aut=aut, pf=pf))

    # --- distributions
    dist = Counter(r['d'] for r in rows)
    print(f'  d-spectrum (per class): {dict(sorted(dist.items()))}')
    if do_aut:
        lab_dist = Counter()
        for r in rows:
            lab_dist[r['d']] += factorial(n) // r['aut']
        print(f'  d-spectrum (labeled):   {dict(sorted(lab_dist.items()))}')
        total = sum(lab_dist.values())
        assert total == 2**comb(n, 2), 'labeled total mismatch'
        mean_num = sum((dv << (n-1)) * cnt for dv, cnt in lab_dist.items())
        inv = involutions(n)
        ok = (mean_num == inv * total)
        print(f'  HYP-2383 involution average: E[det] = {mean_num}/{total} '
              f'= {mean_num/total:.6f} vs I({n}) = {inv}  -> {"CONFIRMED" if ok else "REFUTED"}')

    # --- floor (HYP-2384)
    floor_ok = True
    n_lt = n_d1 = 0
    for r in rows:
        d1 = (r['d'] == 1)
        sw = is_switch_of_transitive(r['A'], n)
        if not (d1 == r['lt'] == sw): floor_ok = False
        n_lt += r['lt']; n_d1 += d1
    print(f'  HYP-2384 floor: d=1 classes = {n_d1}, locally-transitive = {n_lt}, '
          f'all-predicates-agree = {floor_ok}')
    if do_aut:
        lab_lt = sum(factorial(n)//r['aut'] for r in rows if r['lt'])
        dfac = 1
        for t in range(2*n-2, 0, -2): dfac *= t
        print(f'    labeled local orders = {lab_lt} vs (2n-2)!! = {dfac}  '
              f'{"OK" if lab_lt == dfac else "MISMATCH"}')

    # --- ceiling (HYP-2385)
    dmax = max(r['d'] for r in rows)
    if n % 2 == 1:
        bound = (n+1)**((n-1)//2)
        attained = (dmax << (n-1)) == bound
        print(f'  HYP-2385 ceiling (odd): max det = {dmax << (n-1)} vs bound (n+1)^((n-1)/2) = {bound}'
              f'   attained={attained}')
        maxers = [r for r in rows if r['d'] == dmax]
        if attained:
            alldrt = all(is_switch_of_DRT(r['A'], n) for r in maxers)
            print(f'    classes attaining max: {len(maxers)}; all switchings-of-DRT: {alldrt}')
    else:
        bound = n**(n//2)
        attained = (dmax << (n-1)) == bound
        print(f'  HYP-2385 ceiling (even): max det = {dmax << (n-1)} vs bound n^(n/2) = {bound}'
              f'   attained={attained}')
        maxers = [r for r in rows if r['d'] == dmax]
        if attained:
            allsc = all(is_skew_conference(r['A'], n) for r in maxers)
            print(f'    classes attaining max: {len(maxers)}; all skew-conference: {allsc}')
        else:
            print(f'    (no skew conference at n={n}: max-d classes = {len(maxers)}, '
                  f'd_max = {dmax}, scores of maxers: {sorted(set(r["scores"] for r in maxers))})')

    # --- alignment (HYP-2386)
    Hmax = max(r['H'] for r in rows)
    argH = [r for r in rows if r['H'] == Hmax]
    argD = [r for r in rows if r['d'] == dmax]
    setH = set(id(r) for r in argH); setD = set(id(r) for r in argD)
    ds = [r['d'] for r in rows]; Hs = [r['H'] for r in rows]
    mean_d = sum(ds)/len(ds); mean_H = sum(Hs)/len(Hs)
    cov = sum((a-mean_d)*(b-mean_H) for a, b in zip(ds, Hs))
    var_d = sum((a-mean_d)**2 for a in ds); var_H = sum((b-mean_H)**2 for b in Hs)
    pear = cov/ (var_d*var_H) ** 0.5 if var_d*var_H > 0 else float('nan')
    print(f'  HYP-2386 alignment: H_max = {Hmax} ({len(argH)} classes), d_max = {dmax} '
          f'({len(argD)} classes); argmax sets equal: {setH == setD}; '
          f'Pearson(d,H) per class = {pear:.4f}')
    print(f'    d of H-maximizers: {sorted(set(r["d"] for r in argH))}; '
          f'H of d-maximizers: {sorted(set(r["H"] for r in argD))}')
    if n % 2 == 0:
        print(f'    |Pf(S)| of H-maximizers: {sorted(set(r["pf"] for r in argH))} '
              f'(HYP-2312 says 1); |Pf(S)| of d-maximizers: {sorted(set(r["pf"] for r in argD))}')

    # --- complement invariance (sanity for merged-metagraph descent)
    bad = 0
    for r in rows[:200]:
        if detIS(complement(r['A'], n), n) != r['det']: bad += 1
    print(f'  complement-invariance d(T)=d(T^op): violations in first 200 classes: {bad}')

    # --- Hermite (HYP-2387), labeled-exhaustive
    if do_hermite and n <= 6:
        tot = [0]*(n+1)
        cnt = 0
        for A in all_tournaments_small(n):
            S = S_of(A, n)
            negS = [[-x for x in row] for row in S]      # det(xI - (-S)) = det(xI + S)
            cp = char_poly_int(negS, n)
            for i in range(n+1): tot[i] += cp[i]
            cnt += 1
        mp = matching_poly_Kn(n)
        avg = [t / cnt for t in tot]
        match = all(abs(avg[i] - mp[i]) < 1e-9 for i in range(n+1))
        print(f'  HYP-2387 Hermite: E[det(xI+S)] coeffs = {avg}')
        print(f'             matching poly K_{n}    = {mp}   -> {"CONFIRMED" if match else "REFUTED"}')

    print(f'  [{time.time()-t0:.1f}s]')
    return rows

# ---------------------------------------------------------------- QR7 detail

def qr7_detail():
    print(f'\n{"="*72}\n  QR_7 detail (THM-213 cross-check + Pfaffian minor census)\n{"="*72}')
    n, QR = 7, {1, 2, 4}
    A = [[1 if i != j and (j-i) % 7 in QR else 0 for j in range(7)] for i in range(7)]
    S = S_of(A, n)
    det = detIS(A, n)
    print(f'  det(I+S) = {det} = 8^3 = {8**3}; d = {det >> 6}')
    # Pfaffian minor census by subset size
    bysize = defaultdict(Counter)
    for k in range(0, n+1, 2):
        for K in itertools.combinations(range(n), k):
            sub = [[S[i][j] for j in K] for i in K]
            if k == 0: pf2 = 1
            else: pf2 = int_det(sub)
            bysize[k][pf2] += 1
    tot = 0
    for k in sorted(bysize):
        print(f'  |K|={k}: Pf^2 distribution {dict(sorted(bysize[k].items()))}')
        tot += sum(pf2*c for pf2, c in bysize[k].items())
    print(f'  sum of Pf^2 over even subsets = {tot}  (should equal det(I+S) = {det})')
    # THM-213: |Pf(S_00)| = 7^((7-3)/4) = 7
    sub = [[S[i][j] for j in range(1, 7)] for i in range(1, 7)]
    print(f'  |Pf(S minus vertex 0)| = {isqrt_exact(int_det(sub))}  (THM-213 predicts 7)')

# ---------------------------------------------------------------- n=8 Monte Carlo average

def mc_average_n8(samples=2_000_000, seed=20260610):
    import random
    rng = random.Random(seed)
    n = 8
    PAIRS = pairs_of(n)
    tot = 0
    for s in range(samples):
        A = [[0]*n for _ in range(n)]
        for (i, j) in PAIRS:
            if rng.getrandbits(1): A[i][j] = 1
            else: A[j][i] = 1
        tot += detIS(A, n)
    mean = tot / samples
    print(f'\n  n=8 Monte Carlo E[det(I+S)] over {samples} samples: {mean:.3f} '
          f'vs I(8) = {involutions(8)} (proof says exact)')

if __name__ == '__main__':
    print('THE DETERMINANT LENS: tournament det(I+S) census  [mac-mini-2026-06-10-S2]')
    print(f'involution numbers I(n), n=2..10: {[involutions(k) for k in range(2, 11)]}')
    for n in range(3, 8):
        census(n, do_aut=True, do_hermite=(n <= 6))
    qr7_detail()
    census(8, do_aut=False, do_hermite=False)
    mc_average_n8(samples=200_000)
    print('\nDONE.')
