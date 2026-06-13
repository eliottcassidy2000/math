#!/usr/bin/env python3
"""
split_inert_analyzer.py — opus-2026-03-16-S73

SPLIT/INERT PRIME ANALYZER FOR CIRCULANT CRYPTANALYSIS

Leverages the project's structural knowledge about circulant matrices (THM-125)
and splitting behavior of primes in cyclotomic fields to analyze security of
QC-LDPC codes (BIKE, HQC) and circulant-based cryptographic constructions.

Key insight: THM-125 says circulant matrices over F_q decompose into independent
eigenspaces when q splits completely in Q(ζ_p). This gives a factor-p speedup
for rank computation. The splitting behavior is controlled by:
    q splits in Q(ζ_p)  ⟺  q ≡ 1 (mod p)
    q inert in Q(ζ_p)   ⟺  ord_p(q) = p-1

For Paley tournaments (p ≡ 3 mod 4), there's additional QR structure:
    √(-p) ∈ F_q  ⟺  q is a QR mod p (Gauss sum factorization)

This tool:
1. Classifies primes by splitting behavior for a given circulant size p
2. Identifies optimal "attack primes" for eigenspace decomposition
3. Computes defense quality rankings (which p values resist decomposition)
4. Detects torsion in chain complexes via two-prime comparison
5. Assesses security of QC-LDPC parameters

Usage:
    python3 split_inert_analyzer.py --analyze-p 23
    python3 split_inert_analyzer.py --defense-ranking
    python3 split_inert_analyzer.py --torsion-detect 7
    python3 split_inert_analyzer.py --qc-ldpc-audit 587
    python3 split_inert_analyzer.py --full-report 23
"""

import argparse
import math
import sys
from collections import defaultdict

# ─────────────────────────────────────────────────────────────────
# Core number theory
# ─────────────────────────────────────────────────────────────────

def is_prime(n):
    """Deterministic primality test for small n."""
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def primes_up_to(N):
    """Sieve of Eratosthenes."""
    sieve = [True] * (N + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(N**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, N + 1, i):
                sieve[j] = False
    return [i for i in range(2, N + 1) if sieve[i]]


def ord_mod(a, n):
    """Multiplicative order of a modulo n. Returns None if gcd(a,n)>1."""
    if math.gcd(a, n) > 1:
        return None
    r, x = 1, a % n
    while x != 1:
        x = (x * a) % n
        r += 1
        if r > n:
            return None
    return r


def legendre(a, p):
    """Legendre symbol (a/p) via Euler's criterion."""
    if a % p == 0:
        return 0
    val = pow(a, (p - 1) // 2, p)
    return val if val == 1 else -1


def jacobi(a, n):
    """Jacobi symbol (a/n) for odd n > 0."""
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be positive odd")
    a = a % n
    result = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a = a % n
    return result if n == 1 else 0


# ─────────────────────────────────────────────────────────────────
# Splitting analysis
# ─────────────────────────────────────────────────────────────────

def classify_prime(q, p):
    """
    Classify prime q's splitting behavior in Q(ζ_p).

    Returns dict with:
        order: multiplicative order of q mod p
        type: 'split', 'inert', or 'partial'
        num_factors: number of prime ideal factors (= (p-1)/order)
        residue_degree: f = order of q mod p
        speedup: factor by which eigenspace decomposition speeds up rank computation
    """
    if q == p:
        return {
            'order': None,
            'type': 'ramified',
            'num_factors': 1,
            'residue_degree': 1,
            'speedup': 1,
            'note': 'q = p (ramified)'
        }

    f = ord_mod(q, p)
    if f is None:
        return {'order': None, 'type': 'ramified', 'num_factors': 1,
                'residue_degree': 1, 'speedup': 1}

    g = (p - 1) // f  # number of prime ideal factors

    if f == 1:
        split_type = 'split'
    elif f == p - 1:
        split_type = 'inert'
    else:
        split_type = 'partial'

    # Speedup: THM-125 gives g independent eigenspaces
    # Each eigenspace is (p/g)-dimensional instead of p-dimensional
    # For Gaussian elimination: O(n^3) → g × O((n/g)^3) = O(n^3/g^2)
    speedup = g  # linear speedup in eigenspace count

    return {
        'order': f,
        'type': split_type,
        'num_factors': g,
        'residue_degree': f,
        'speedup': speedup
    }


def splitting_table(p, q_max=200):
    """
    Build complete splitting table for Q(ζ_p) over primes q ≤ q_max.

    Returns list of (q, classification) sorted by speedup (descending).
    """
    results = []
    for q in primes_up_to(q_max):
        if q == p:
            continue
        cl = classify_prime(q, p)
        cl['q'] = q
        results.append(cl)

    results.sort(key=lambda x: (-x['speedup'], x['q']))
    return results


def attack_primes(p, q_max=1000, top_k=10):
    """
    Find the best primes q for attacking p-circulant codes.

    Best = fully splitting (q ≡ 1 mod p), giving p-1 independent eigenspaces.
    """
    best = []
    for q in primes_up_to(q_max):
        if q == p:
            continue
        if q % p == 1:
            best.append(q)
            if len(best) >= top_k:
                break
    return best


def defense_ranking(p_max=200, q_scan=1000):
    """
    Rank primes p by how well they resist eigenspace decomposition.

    Defense quality = smallest q that splits completely in Q(ζ_p).
    Higher = better defense (attacker needs larger field).
    """
    results = []
    for p in primes_up_to(p_max):
        if p < 5:
            continue
        # Find smallest splitting prime
        smallest_split = None
        for q in primes_up_to(q_scan):
            if q == p:
                continue
            if q % p == 1:
                smallest_split = q
                break
        if smallest_split is None:
            smallest_split = float('inf')

        # Frobenius density: fraction of primes that split
        split_count = 0
        total = 0
        for q in primes_up_to(q_scan):
            if q == p:
                continue
            total += 1
            if q % p == 1:
                split_count += 1
        density = split_count / total if total > 0 else 0

        results.append({
            'p': p,
            'smallest_split_q': smallest_split,
            'split_density': density,
            'theoretical_density': 1.0 / (p - 1),
            'is_paley': p % 4 == 3
        })

    results.sort(key=lambda x: -x['smallest_split_q'])
    return results


# ─────────────────────────────────────────────────────────────────
# Paley tournament construction and boundary matrices
# ─────────────────────────────────────────────────────────────────

def paley_tournament(p):
    """
    Construct Paley tournament on p vertices (p ≡ 3 mod 4 prime).

    A[i][j] = 1 iff (j-i) is a quadratic residue mod p.
    """
    if p % 4 != 3 or not is_prime(p):
        raise ValueError(f"p={p} must be prime ≡ 3 mod 4")

    QR = set()
    for k in range(1, p):
        QR.add(pow(k, 2, p))

    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in QR:
                A[i][j] = 1
    return A


def boundary_d2(A, n):
    """
    Build ∂_2 boundary matrix for path homology.

    ∂_2 maps allowed 3-paths to allowed 2-paths.
    An allowed path (v0,v1,...,vk) requires A[vi][vi+1]=1 for all i
    and all faces are allowed.
    """
    # Enumerate allowed 2-paths (edges)
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(n):
            if i != j and A[i][j]:
                edge_idx[(i, j)] = len(edges)
                edges.append((i, j))

    # Enumerate allowed 3-paths (i→j→k with both edges present)
    paths3 = []
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]:
                continue
            for k in range(n):
                if k == i or k == j or not A[j][k]:
                    continue
                # Check all faces are allowed
                # Face 0: (j,k) — already checked (A[j][k]=1)
                # Face 1: (i,k) — need A[i][k]=1 OR A[k][i]=1?
                # For tournaments, exactly one of (i,k) or (k,i) exists
                # Face 1 of (i,j,k): skip j → path (i,k)
                # This face must be in A_1 (allowed 1-paths = edges)
                # In a tournament, (i,k) is always an edge (either direction)
                # But for ∂_2, face 1 is (i,k) which must be an allowed 1-path
                # Every pair has an edge in a tournament, so face 1 is always allowed
                paths3.append((i, j, k))

    if not paths3 or not edges:
        return [], [], [[]]

    # Build boundary matrix: ∂_2(i,j,k) = (j,k) - (i,k) + (i,j)
    nrows = len(edges)
    ncols = len(paths3)
    M = [[0] * ncols for _ in range(nrows)]

    for col, (i, j, k) in enumerate(paths3):
        # Face 0: (j,k) with sign +1
        if (j, k) in edge_idx:
            M[edge_idx[(j, k)]][col] += 1
        # Face 1: (i,k) with sign -1
        if (i, k) in edge_idx:
            M[edge_idx[(i, k)]][col] -= 1
        elif (k, i) in edge_idx:
            # In tournament, if A[i][k]=0 then A[k][i]=1
            # But (i,k) as a face means the path i→k
            # If i→k is not an edge, this face doesn't map to an allowed path
            pass
        # Face 2: (i,j) with sign +1
        if (i, j) in edge_idx:
            M[edge_idx[(i, j)]][col] += 1

    return edges, paths3, M


def rank_mod_q(M, q):
    """
    Compute rank of integer matrix M over F_q via Gaussian elimination.

    M is a list of lists. Works in-place on a copy.
    """
    if not M or not M[0]:
        return 0

    nrows = len(M)
    ncols = len(M[0])
    # Copy and reduce mod q
    W = [[M[i][j] % q for j in range(ncols)] for i in range(nrows)]

    rank = 0
    for col in range(ncols):
        # Find pivot
        pivot = None
        for row in range(rank, nrows):
            if W[row][col] % q != 0:
                pivot = row
                break
        if pivot is None:
            continue
        # Swap
        W[rank], W[pivot] = W[pivot], W[rank]
        # Scale pivot row
        inv = pow(W[rank][col], q - 2, q)  # Fermat's little theorem
        for j in range(ncols):
            W[rank][j] = (W[rank][j] * inv) % q
        # Eliminate
        for row in range(nrows):
            if row == rank:
                continue
            if W[row][col] % q != 0:
                factor = W[row][col]
                for j in range(ncols):
                    W[row][j] = (W[row][j] - factor * W[rank][j]) % q
        rank += 1

    return rank


# ─────────────────────────────────────────────────────────────────
# Torsion detection via two-prime method
# ─────────────────────────────────────────────────────────────────

def boundary_d3(A, n):
    """
    Build ∂_3 boundary matrix for path homology.

    ∂_3 maps allowed 4-paths to allowed 3-paths.
    """
    # Enumerate allowed 3-paths
    paths3 = []
    path3_idx = {}
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]:
                continue
            for k in range(n):
                if k in (i, j) or not A[j][k]:
                    continue
                path3_idx[(i, j, k)] = len(paths3)
                paths3.append((i, j, k))

    # Enumerate allowed 4-paths (i→j→k→l)
    paths4 = []
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]:
                continue
            for k in range(n):
                if k in (i, j) or not A[j][k]:
                    continue
                for l in range(n):
                    if l in (i, j, k) or not A[k][l]:
                        continue
                    paths4.append((i, j, k, l))

    if not paths4 or not paths3:
        return paths3, paths4, [[]]

    # Build ∂_3: faces of (i,j,k,l) are:
    # face 0: (j,k,l) with sign +1
    # face 1: (i,k,l) with sign -1
    # face 2: (i,j,l) with sign +1
    # face 3: (i,j,k) with sign -1
    nrows = len(paths3)
    ncols = len(paths4)
    M = [[0] * ncols for _ in range(nrows)]

    for col, (i, j, k, l) in enumerate(paths4):
        for face_num, face in enumerate([(j, k, l), (i, k, l), (i, j, l), (i, j, k)]):
            sign = 1 if face_num % 2 == 0 else -1
            if face in path3_idx:
                M[path3_idx[face]][col] += sign

    return paths3, paths4, M


def detect_torsion(p_tournament, q_list=None, d=2):
    """
    Detect torsion in path homology chain complex of Paley tournament P_p.

    Method: compute rank(∂_d) over multiple F_q. If rank drops at some q,
    then the chain complex has q-torsion, revealing algebraic structure.

    Returns dict mapping q → rank, plus detected torsion primes.
    """
    if q_list is None:
        q_list = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

    # Remove p itself (ramified)
    q_list = [q for q in q_list if q != p_tournament and is_prime(q)]

    A = paley_tournament(p_tournament)
    n = p_tournament

    if d == 2:
        edges, paths3, M = boundary_d2(A, n)
        if not M or not M[0]:
            return {'ranks': {}, 'torsion_primes': [], 'matrix_size': (0, 0)}
    elif d == 3:
        paths3, paths4, M = boundary_d3(A, n)
        if not M or not M[0]:
            return {'ranks': {}, 'torsion_primes': [], 'matrix_size': (0, 0)}
    else:
        return {'error': f'∂_{d} not implemented (only d=2,3 currently)'}

    # Compute rank at each prime
    ranks = {}
    for q in q_list:
        ranks[q] = rank_mod_q(M, q)

    # Generic rank = maximum observed
    generic = max(ranks.values()) if ranks else 0

    # Torsion primes = where rank drops
    torsion = [q for q, r in ranks.items() if r < generic]

    return {
        'ranks': ranks,
        'generic_rank': generic,
        'torsion_primes': torsion,
        'matrix_size': (len(M), len(M[0]) if M else 0),
        'num_edges': len(edges) if d == 2 else None,
        'num_3paths': len(paths3) if d == 2 else None
    }


# ─────────────────────────────────────────────────────────────────
# QC-LDPC security assessment
# ─────────────────────────────────────────────────────────────────

def assess_qc_ldpc(p, security_bits=128):
    """
    Assess security of a QC-LDPC code with circulant block size p.

    Checks:
    1. How many small primes q split completely (give max decomposition)
    2. Defense quality (smallest splitting q)
    3. Eigenspace dimension reduction factor
    4. Estimated attack speedup vs generic
    """
    # Find splitting primes
    split_qs = attack_primes(p, q_max=10000, top_k=50)

    # Smallest splitting prime
    smallest_split = split_qs[0] if split_qs else None

    # Count splitting primes below various bounds
    split_below = {}
    for bound in [100, 1000, 10000]:
        split_below[bound] = len([q for q in split_qs if q <= bound])

    # Partial decomposition primes (order divides p-1 but > 1)
    partial_qs = []
    for q in primes_up_to(min(p * 10, 5000)):
        if q == p:
            continue
        f = ord_mod(q, p)
        if f is not None and 1 < f < p - 1:
            g = (p - 1) // f
            partial_qs.append({'q': q, 'order': f, 'factors': g, 'speedup': g})

    partial_qs.sort(key=lambda x: -x['speedup'])

    # Theoretical Chebotarev density
    cheb_density = 1.0 / (p - 1)

    # Eigenspace dimension without decomposition
    full_dim = p

    # Security assessment
    if smallest_split and smallest_split <= 2:
        risk = 'CRITICAL'
        detail = f'q=2 splits! Binary eigenspace decomposition available.'
    elif smallest_split and smallest_split <= 64:
        risk = 'HIGH'
        detail = f'Small splitting prime q={smallest_split} enables efficient decomposition.'
    elif smallest_split and smallest_split <= 1024:
        risk = 'MEDIUM'
        detail = f'Splitting prime q={smallest_split} available but requires moderate field.'
    elif partial_qs and partial_qs[0]['speedup'] >= p // 4:
        risk = 'MEDIUM'
        detail = f'No small full split, but partial decomposition gives {partial_qs[0]["speedup"]}× speedup at q={partial_qs[0]["q"]}.'
    else:
        risk = 'LOW'
        detail = f'Smallest splitting prime q={smallest_split}. Good defense.'

    return {
        'p': p,
        'is_prime': is_prime(p),
        'is_paley': is_prime(p) and p % 4 == 3,
        'smallest_split_q': smallest_split,
        'split_primes_below_100': split_below.get(100, 0),
        'split_primes_below_1000': split_below.get(1000, 0),
        'split_primes_below_10000': split_below.get(10000, 0),
        'chebotarev_density': cheb_density,
        'best_partial': partial_qs[:5] if partial_qs else [],
        'full_dimension': full_dim,
        'risk_level': risk,
        'risk_detail': detail,
        'attack_primes': split_qs[:10]
    }


# ─────────────────────────────────────────────────────────────────
# Reporting
# ─────────────────────────────────────────────────────────────────

def print_splitting_table(p, q_max=200):
    """Print formatted splitting table."""
    print(f"\n{'='*70}")
    print(f"SPLITTING TABLE: Q(ζ_{p}) over primes q ≤ {q_max}")
    print(f"{'='*70}")
    print(f"{'q':>5} {'ord_p(q)':>9} {'type':>10} {'#factors':>9} {'speedup':>8}")
    print(f"{'-'*5:>5} {'-'*9:>9} {'-'*10:>10} {'-'*9:>9} {'-'*8:>8}")

    table = splitting_table(p, q_max)
    for cl in table:
        q = cl['q']
        print(f"{q:>5} {str(cl['order']):>9} {cl['type']:>10} {cl['num_factors']:>9} {cl['speedup']:>8}×")

    # Summary
    split_count = sum(1 for cl in table if cl['type'] == 'split')
    inert_count = sum(1 for cl in table if cl['type'] == 'inert')
    partial_count = sum(1 for cl in table if cl['type'] == 'partial')
    total = len(table)

    print(f"\n--- Summary ---")
    print(f"Total primes: {total}")
    print(f"Split (q≡1 mod {p}): {split_count} ({100*split_count/total:.1f}%)")
    print(f"Inert (ord={p-1}): {inert_count} ({100*inert_count/total:.1f}%)")
    print(f"Partial: {partial_count} ({100*partial_count/total:.1f}%)")
    print(f"Chebotarev prediction: {100/(p-1):.1f}% split")

    return table


def print_defense_ranking(p_max=200, top_k=20):
    """Print defense quality ranking."""
    print(f"\n{'='*70}")
    print(f"DEFENSE RANKING: Primes p ≤ {p_max} by resistance to decomposition")
    print(f"{'='*70}")
    print(f"{'rank':>4} {'p':>5} {'min_split_q':>12} {'split%':>8} {'Cheb%':>8} {'Paley':>6}")
    print(f"{'-'*4:>4} {'-'*5:>5} {'-'*12:>12} {'-'*8:>8} {'-'*8:>8} {'-'*6:>6}")

    results = defense_ranking(p_max)
    for i, r in enumerate(results[:top_k]):
        sq = r['smallest_split_q']
        sq_str = str(sq) if sq != float('inf') else '∞'
        paley = '✓' if r['is_paley'] else ''
        print(f"{i+1:>4} {r['p']:>5} {sq_str:>12} "
              f"{100*r['split_density']:>7.2f}% "
              f"{100*r['theoretical_density']:>7.2f}% "
              f"{paley:>6}")

    return results


def print_torsion_report(p_tournament, q_list=None, check_d3=True):
    """Print torsion detection report."""
    print(f"\n{'='*70}")
    print(f"TORSION DETECTION: Paley tournament P_{p_tournament}")
    print(f"{'='*70}")

    all_results = {}
    for d in [2, 3] if check_d3 else [2]:
        result = detect_torsion(p_tournament, q_list, d=d)
        all_results[d] = result

        if 'error' in result:
            print(f"\n∂_{d}: {result['error']}")
            continue

        print(f"\n--- ∂_{d} ---")
        print(f"Matrix size: {result['matrix_size'][0]} × {result['matrix_size'][1]}")
        print(f"\n{'q':>5} {'rank':>6} {'drop':>6}")
        print(f"{'-'*5:>5} {'-'*6:>6} {'-'*6:>6}")

        generic = result['generic_rank']
        for q in sorted(result['ranks'].keys()):
            r = result['ranks'][q]
            drop = generic - r
            marker = ' ← TORSION' if drop > 0 else ''
            print(f"{q:>5} {r:>6} {drop:>6}{marker}")

        print(f"\nGeneric rank: {generic}")
        if result['torsion_primes']:
            print(f"TORSION DETECTED at primes: {result['torsion_primes']}")
            print(f"  → These primes reveal algebraic structure in the chain complex")
            print(f"  → Attackers can exploit rank drops for distinguishing/decoding")
        else:
            print(f"No torsion detected at ∂_{d}")

    return all_results


def print_qc_ldpc_audit(p):
    """Print QC-LDPC security assessment."""
    print(f"\n{'='*70}")
    print(f"QC-LDPC SECURITY AUDIT: circulant block size p = {p}")
    print(f"{'='*70}")

    result = assess_qc_ldpc(p)

    print(f"\n--- Basic Properties ---")
    print(f"p = {p} ({'prime' if result['is_prime'] else 'NOT PRIME'})")
    print(f"Paley tournament: {'Yes (p ≡ 3 mod 4)' if result['is_paley'] else 'No'}")
    print(f"Full eigenspace dimension: {result['full_dimension']}")

    print(f"\n--- Splitting Analysis ---")
    print(f"Smallest fully-splitting prime: q = {result['smallest_split_q']}")
    print(f"Splitting primes below 100: {result['split_primes_below_100']}")
    print(f"Splitting primes below 1000: {result['split_primes_below_1000']}")
    print(f"Splitting primes below 10000: {result['split_primes_below_10000']}")
    print(f"Chebotarev density: {result['chebotarev_density']:.6f} "
          f"(≈ 1 in {int(1/result['chebotarev_density'])} primes)")

    if result['attack_primes']:
        print(f"\n--- Attack Primes (q ≡ 1 mod {p}) ---")
        for q in result['attack_primes']:
            print(f"  q = {q}  (speedup: {p-1}×)")

    if result['best_partial']:
        print(f"\n--- Best Partial Decompositions ---")
        for pp in result['best_partial']:
            print(f"  q = {pp['q']:>5}, ord_{p}(q) = {pp['order']}, "
                  f"factors = {pp['factors']}, speedup = {pp['speedup']}×")

    print(f"\n{'='*70}")
    print(f"RISK LEVEL: {result['risk_level']}")
    print(f"  {result['risk_detail']}")
    print(f"{'='*70}")

    # Recommendations
    print(f"\n--- Recommendations ---")
    if result['risk_level'] in ('CRITICAL', 'HIGH'):
        print(f"  ⚠ AVOID using p={p} for QC-LDPC codes")
        print(f"  ⚠ Eigenspace decomposition gives up to {p-1}× speedup")
        print(f"  ⚠ Consider larger p with better defense profile")
    elif result['risk_level'] == 'MEDIUM':
        print(f"  △ p={p} has moderate resistance")
        print(f"  △ Increase security parameter to compensate for partial decomposition")
    else:
        print(f"  ✓ p={p} has good resistance to eigenspace decomposition")
        print(f"  ✓ Smallest splitting prime is sufficiently large")

    return result


def print_full_report(p):
    """Print comprehensive analysis for a given prime p."""
    print(f"\n{'#'*70}")
    print(f"# FULL SPLIT/INERT ANALYSIS: p = {p}")
    print(f"# Generated by split_inert_analyzer.py (opus-2026-03-16-S73)")
    print(f"{'#'*70}")

    # 1. Splitting table
    print_splitting_table(p, q_max=200)

    # 2. QC-LDPC audit
    print_qc_ldpc_audit(p)

    # 3. Torsion detection (only for Paley primes, limited by pure Python speed)
    # ∂_2: feasible up to p~31, ∂_3: feasible up to p~19 in pure Python
    if is_prime(p) and p % 4 == 3 and p <= 19:
        print_torsion_report(p, check_d3=True)
    elif is_prime(p) and p % 4 == 3 and p <= 31:
        print_torsion_report(p, check_d3=False)
        print(f"\n[∂_3 torsion detection skipped: P_{p} too large for pure Python]")
    elif is_prime(p) and p % 4 == 3:
        print(f"\n[Torsion detection skipped: P_{p} too large for pure Python]")
        print(f"[Use NumPy version for p > 31]")


# ─────────────────────────────────────────────────────────────────
# Demo / self-test
# ─────────────────────────────────────────────────────────────────

def run_demo():
    """Run demonstration analysis."""
    print("SPLIT/INERT PRIME ANALYZER — Demo")
    print("="*50)

    # 1. Analyze a few important primes
    for p in [7, 11, 23, 47, 587]:
        if not is_prime(p):
            continue
        result = assess_qc_ldpc(p)
        paley = "Paley" if result['is_paley'] else "non-Paley"
        print(f"\np={p:>4} ({paley}): risk={result['risk_level']:<8} "
              f"min_split_q={result['smallest_split_q']:<6} "
              f"Cheb={result['chebotarev_density']:.4f}")

    # 2. Defense ranking
    print_defense_ranking(p_max=100, top_k=10)

    # 3. Torsion detection on small Paley tournaments
    for p in [7, 11, 19, 23]:
        print_torsion_report(p)

    # 4. BIKE/HQC parameter check
    print(f"\n{'='*70}")
    print(f"POST-QUANTUM PARAMETER CHECK")
    print(f"{'='*70}")
    # BIKE uses circulant block sizes around r = 12323 (Level 1)
    # HQC uses n = 17669 (Level 1)
    # These are too large for full analysis, but we can check splitting
    for name, p_val in [("BIKE-L1 (r≈12323)", 12323),
                         ("HQC-L1 (n≈17669)", 17669)]:
        if is_prime(p_val):
            splits = attack_primes(p_val, q_max=100000, top_k=3)
            print(f"\n{name}: p={p_val}")
            if splits:
                print(f"  Smallest splitting q: {splits[0]}")
                print(f"  Chebotarev density: 1/{p_val-1} ≈ {1/(p_val-1):.2e}")
            else:
                print(f"  No splitting prime found below 100000")
        else:
            # Factor it
            factors = []
            temp = p_val
            for d in range(2, int(temp**0.5) + 1):
                while temp % d == 0:
                    factors.append(d)
                    temp //= d
            if temp > 1:
                factors.append(temp)
            print(f"\n{name}: p={p_val} = {'×'.join(map(str, factors))} (COMPOSITE)")
            print(f"  Analyze each prime factor separately for circulant structure")


# ─────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Split/Inert Prime Analyzer for Circulant Cryptanalysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --analyze-p 23         Splitting table for Q(ζ_23)
  %(prog)s --defense-ranking      Rank primes by decomposition resistance
  %(prog)s --torsion-detect 7     Detect torsion in P_7 chain complex
  %(prog)s --qc-ldpc-audit 587    Security assessment for p=587
  %(prog)s --full-report 23       Complete analysis for p=23
  %(prog)s --demo                 Run demonstration
        """)

    parser.add_argument('--analyze-p', type=int,
                        help='Print splitting table for Q(ζ_p)')
    parser.add_argument('--defense-ranking', action='store_true',
                        help='Print defense quality ranking')
    parser.add_argument('--defense-max', type=int, default=200,
                        help='Max p for defense ranking (default: 200)')
    parser.add_argument('--torsion-detect', type=int,
                        help='Detect torsion in P_p chain complex')
    parser.add_argument('--qc-ldpc-audit', type=int,
                        help='QC-LDPC security assessment for block size p')
    parser.add_argument('--full-report', type=int,
                        help='Full analysis for prime p')
    parser.add_argument('--demo', action='store_true',
                        help='Run demonstration')
    parser.add_argument('--q-max', type=int, default=200,
                        help='Max q for splitting table (default: 200)')

    args = parser.parse_args()

    if args.demo:
        run_demo()
    elif args.analyze_p:
        print_splitting_table(args.analyze_p, args.q_max)
    elif args.defense_ranking:
        print_defense_ranking(args.defense_max)
    elif args.torsion_detect:
        print_torsion_report(args.torsion_detect)
    elif args.qc_ldpc_audit:
        print_qc_ldpc_audit(args.qc_ldpc_audit)
    elif args.full_report:
        print_full_report(args.full_report)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
