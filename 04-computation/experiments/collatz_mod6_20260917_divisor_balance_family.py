#!/usr/bin/env python3
"""Divisor-balance family: density, linear cells F = a*S + b*U, k-free refinement, Omega blindness.

Lane: divisor_balance_family (session collatz-mod6-20260917).
Inherits DB1-DB3 of 05-knowledge/results/arithmetic_braids_20260917_divisors.md:
  F = prod(a_i+1) - 2,  S = 2^r - 1 - [N squarefree],  U = r - [N prime],
  F = S + U  iff  N in {p, p^3, p^2 q r}.

Run from anywhere (absolute paths only):
  python3 04-computation/experiments/collatz_mod6_20260917_divisor_balance_family.py
  python3 -O 04-computation/experiments/collatz_mod6_20260917_divisor_balance_family.py

Dependencies: python3 + numpy (sieve only). All load-bearing arithmetic is exact
integer arithmetic; every check raises RuntimeError explicitly (active under -O).
RAM < 600 MB, wall time ~1-2 minutes on one core.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from math import isqrt, prod

import numpy as np

T0 = time.time()


def require(condition, message):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + message)


def banner(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


# ----------------------------------------------------------------------------
# 0. Exact definitions (profile formulas, inherited DB1) and direct divisor audit
# ----------------------------------------------------------------------------

def spf_sieve(limit):
    """Smallest prime factor for 0..limit (pure python, used for direct audits)."""
    spf = list(range(limit + 1))
    for i in range(2, isqrt(limit) + 1):
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factor_spf(n, spf):
    out = {}
    while n > 1:
        p = spf[n]
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    return out


def divisors_from_factors(fac):
    divs = [1]
    for p, a in fac.items():
        divs = [d * p ** b for d in divs for b in range(a + 1)]
    return divs


def direct_counts(n, spf, kmax=6):
    """Direct divisor-set counts: F, S_k for k=1..kmax (S_2 = S), U.  Exact."""
    fac = factor_spf(n, spf)
    divs = divisors_from_factors(fac)
    proper = [d for d in divs if 1 < d < n]
    F = len(proper)
    U = sum(1 for d in proper if d in fac)  # primes p < n dividing n
    Sk = {}
    for k in range(1, kmax + 1):
        Sk[k] = sum(1 for d in proper if all(d % (p ** k) for p in fac))
    return F, Sk, U


def profile_F(prof):
    return prod(a + 1 for a in prof) - 2


def profile_S(prof):
    r = len(prof)
    sqfree = all(a == 1 for a in prof)
    return (1 << r) - 1 - (1 if sqfree else 0)


def profile_U(prof):
    return len(prof) - (1 if tuple(prof) == (1,) else 0)


def profile_Sk(prof, k):
    """# k-free proper nontrivial divisors: box prod min(a_i+1,k) minus 1 minus [N k-free]."""
    kfree = all(a < k for a in prof)
    return prod(min(a + 1, k) for a in prof) - 1 - (1 if kfree else 0)


def profile_of(fac):
    return tuple(sorted(fac.values(), reverse=True))


banner("0. Definitions and direct audit of the profile formulas (FINITE-EXACT)")
AUDIT_N = 200_000
spf_small = spf_sieve(AUDIT_N)
mismatch = 0
for n in range(2, AUDIT_N + 1):
    fac = factor_spf(n, spf_small)
    prof = profile_of(fac)
    F, Sk, U = direct_counts(n, spf_small, kmax=6)
    require(F == profile_F(prof), f"F formula at {n}")
    require(Sk[2] == profile_S(prof), f"S formula at {n}")
    require(U == profile_U(prof), f"U formula at {n}")
    for k in range(1, 7):
        require(Sk[k] == profile_Sk(prof, k), f"S_{k} formula at {n}")
print(f"profile formulas F, S=S_2, U, S_1..S_6 agree with direct divisor enumeration for all 2<=N<={AUDIT_N}")
print(f"  [time {time.time()-T0:.1f}s]")


# ----------------------------------------------------------------------------
# 1. Density of F = S + U solutions by shape up to 10^7, crossover (FINITE-EXACT)
# ----------------------------------------------------------------------------

banner("1. Density by shape: primes, prime cubes, p^2qr up to 10^7 (FINITE-EXACT)")
X = 10_000_000
is_prime = np.ones(X + 1, dtype=bool)
is_prime[:2] = False
for i in range(2, isqrt(X) + 1):
    if is_prime[i]:
        is_prime[i * i::i] = False
primes = np.nonzero(is_prime)[0]
print(f"pi({X}) = {len(primes)}")

omega = np.zeros(X + 1, dtype=np.int8)
Omega = np.zeros(X + 1, dtype=np.int8)
for p in primes.tolist():
    omega[p::p] += 1
    pk = p
    while pk <= X:
        Omega[pk::pk] += 1
        pk *= p
print(f"omega/Omega sieve complete [time {time.time()-T0:.1f}s]")

# hostile control: sieve vs trial factorization on 2..AUDIT_N
for n in range(2, AUDIT_N + 1):
    fac = factor_spf(n, spf_small)
    require(int(omega[n]) == len(fac), f"omega sieve at {n}")
    require(int(Omega[n]) == sum(fac.values()), f"Omega sieve at {n}")
print(f"sieve omega/Omega agree with trial factorization for 2<=N<={AUDIT_N}")

shape_prime = (Omega == 1)
shape_cube = (Omega == 3) & (omega == 1)
shape_p2qr = (Omega == 4) & (omega == 3)   # Omega=4, omega=3 forces profile (2,1,1)
cum_prime = np.cumsum(shape_prime, dtype=np.int64)
cum_cube = np.cumsum(shape_cube, dtype=np.int64)
cum_p2qr = np.cumsum(shape_p2qr, dtype=np.int64)

# direct-definition control on 2..AUDIT_N: F=S+U solutions == union of the three shapes
for n in range(2, AUDIT_N + 1):
    prof = profile_of(factor_spf(n, spf_small))
    sol = profile_F(prof) == profile_S(prof) + profile_U(prof)
    shape = bool(shape_prime[n] or shape_cube[n] or shape_p2qr[n])
    require(sol == shape, f"F=S+U shape membership at {n}")
    require(shape_p2qr[n] == (prof == (2, 1, 1)), f"p2qr indicator at {n}")
print(f"F=S+U <=> shape in {{p, p^3, p^2qr}} verified directly for 2<=N<={AUDIT_N}")


def pi_of(y):
    return int(cum_prime[y]) if y >= 0 else 0


def semiprime_sqfree_count(y):
    """#{q<r primes: q*r <= y} via prime sums (independent of the Omega sieve)."""
    total = 0
    for q in primes:
        q = int(q)
        if q * q > y:
            break
        total += pi_of(y // q) - pi_of(q)
    return total


def p2qr_count_formula(x):
    """#{p^2 q r <= x, p,q,r distinct primes} = sum_p [ pi2sf(x/p^2) - #{r != p: p r <= x/p^2} ]."""
    total = 0
    for p in primes:
        p = int(p)
        y = x // (p * p)
        if y < 6:
            break
        sub = pi_of(y // p) - (1 if p <= y // p else 0)
        total += semiprime_sqfree_count(y) - sub
    return total


P2_partial = float(np.sum(1.0 / (primes.astype(np.float64) ** 2)))
P2_TAIL_BOUND = 1.0 / X   # sum_{n > X} 1/n^2 < 1/X
print(f"prime zeta P(2) partial sum over p<=10^7 = {P2_partial:.10f} (tail < {P2_TAIL_BOUND:.1e}); literature value 0.4522474200...")

print()
print("x         | #primes  | #p^3 | #p^2qr   | formula  | p2qr/primes | P(2) x loglog x/log x | asym/exact")
print("----------+----------+------+----------+----------+-------------+-----------------------+-----------")
density_rows = []
for x in (10 ** 5, 10 ** 6, 10 ** 7):
    n1, n3, n4 = int(cum_prime[x]), int(cum_cube[x]), int(cum_p2qr[x])
    f4 = p2qr_count_formula(x)
    require(f4 == n4, f"p2qr formula vs sieve at x={x}: {f4} vs {n4}")
    require(n3 == pi_of(int(round(x ** (1 / 3))) if int(round(x ** (1 / 3))) ** 3 <= x else int(round(x ** (1 / 3))) - 1),
            f"cube count vs pi(x^(1/3)) at {x}")
    asym = 0.4522474200 * x * math.log(math.log(x)) / math.log(x)
    density_rows.append((x, n1, n3, n4))
    print(f"{x:<9} | {n1:>8} | {n3:>4} | {n4:>8} | {f4:>8} | {n4/n1:>11.4f} | {asym:>21.0f} | {asym/n4:.4f}")

# exact crossover
diff = cum_p2qr - cum_prime
first = int(np.argmax(diff > 0))
require(diff[first] > 0, "crossover not found below 10^7")
require(bool(shape_p2qr[first]), "crossover point must itself be a p^2qr number")
last_prime_lead = int(np.nonzero(diff <= 0)[0][-1])
# (a) flips of the PREDICATE diff > 0 (a drop from lead 1 to a tie counts as a flip)
pred = diff > 0
flip_pts = (np.nonzero(pred[1:] != pred[:-1])[0] + 1).tolist()
flip_kinds = "".join("Q" if shape_p2qr[n] else ("P" if shape_prime[n] else "?") for n in flip_pts)
# (b) TRUE sign changes of the difference: consecutive nonzero values of opposite sign (ties skipped)
nz = np.nonzero(diff)[0]
sg = np.sign(diff[nz])
sign_change_pts = [int(nz[i + 1]) for i in np.nonzero(sg[1:] != sg[:-1])[0]]
sign_change_dirs = [(int(sg[i]), int(sg[i + 1])) for i in np.nonzero(sg[1:] != sg[:-1])[0]]
neg_after = np.nonzero(diff[first:] < 0)[0] + first
ties_after = np.nonzero(diff[first:] == 0)[0] + first
require(flip_pts == [145119, 145121, 145132, 145133, 145138, 145139, 145148, 147709, 147725, 147727, 147908, 147919, 147925],
        f"predicate flip points: {flip_pts}")
require(flip_kinds == "QPQPQPQPQPQPQ", f"flip kinds alternate p^2qr / prime: {flip_kinds}")
require(sign_change_pts == [145119, 147739, 147908], f"true sign changes of the difference: {sign_change_pts}")
require(sign_change_dirs == [(-1, 1), (1, -1), (-1, 1)], f"sign change directions: {sign_change_dirs}")
require(len(neg_after) == 155 and int(neg_after.min()) == 147739 and int(neg_after.max()) == 147893,
        f"primes strictly ahead after the crossover: {len(neg_after)} points in [{neg_after.min()},{neg_after.max()}]")
require(bool((np.diff(neg_after) == 1).all()), "primes-ahead set after the crossover is one interval")
require(int(diff[first:].min()) == -5, f"minimum of the difference after the crossover: {int(diff[first:].min())}")
require(len(ties_after) == 73 and int(ties_after.max()) == 147924 == last_prime_lead,
        f"ties after the crossover: {len(ties_after)}, last {int(ties_after.max())}")
require(bool((diff[last_prime_lead + 1:] >= 1).all()), "p^2qr leads strictly from 147925 to 10^7")
fac_first = factor_spf(first, spf_sieve(first))
print()
print(f"CROSSOVER (FINITE-EXACT): first n with #{{p^2qr<=n}} > #{{primes<=n}} is n = {first} = {fac_first}")
print(f"  counts there: primes={int(cum_prime[first])}, p^2qr={int(cum_p2qr[first])}")
print(f"  the PREDICATE '#p^2qr > #primes' flips {len(flip_pts)} times on [1,10^7], at {flip_pts}")
print(f"    (kinds {flip_kinds}: every second flip is a prime dropping the lead to a TIE, not a sign change)")
print(f"  TRUE sign changes of the difference #p^2qr - #primes on [1,10^7]: {len(sign_change_pts)} at {sign_change_pts}, directions {sign_change_dirs}")
print(f"  primes strictly ahead after {first}: exactly the {len(neg_after)} integers of [{int(neg_after.min())},{int(neg_after.max())}], minimum difference {int(diff[first:].min())}")
print(f"  ties after {first}: {len(ties_after)}, the last at {int(ties_after.max())}; last n <= 10^7 with #primes >= #p^2qr: {last_prime_lead}")
print(f"  p^2qr leads strictly on [{last_prime_lead+1},10^7]; at 10^7 the lead is p^2qr - primes = {int(diff[X])}")
naive_scale = math.exp(math.exp(1.0 / 0.4522474200))   # solve P(2) loglog x = 1
print(f"  HEURISTIC naive scale from P(2) loglog x = 1: x = {naive_scale:.0f}; actual crossover / naive = {first / naive_scale:.2f}")
del pred, nz, sg, neg_after, ties_after
print(f"  [time {time.time()-T0:.1f}s]")


# ----------------------------------------------------------------------------
# 2. Linear family F = alpha*S + beta*U, (alpha,beta) in {0,1,2,3}^2 (PROVED + FINITE-EXACT)
# ----------------------------------------------------------------------------

banner("2. Linear cells F = alpha*S + beta*U on exponent profiles")


def factorizations(T, r, minf=2):
    """Unordered factorizations of T into exactly r integer factors >= minf (nondecreasing)."""
    if r == 1:
        return [[T]] if T >= minf else []
    out = []
    f = minf
    while f ** r <= T:
        if T % f == 0:
            for rest in factorizations(T // f, r - 1, f):
                out.append([f] + rest)
        f += 1
    return out


def target_T(alpha, beta, r):
    """Nonsquarefree equation: prod(a_i+1) = 2 + alpha(2^r-1) + beta r."""
    return 2 + alpha * ((1 << r) - 1) + beta * r


def support_bound(alpha, beta):
    """PROVED bound on support r of nonsquarefree solutions (None = infinite cell)."""
    if alpha in (0, 1):
        # need 3*2^(r-1) <= T(r); T grows at most like 2^r + O(r), so bounded; compute
        r = 1
        while 3 * (1 << (r - 1)) <= target_T(alpha, beta, r):
            r += 1
        return r - 1
    if (alpha, beta) == (2, 0):
        return None
    return 7  # Lemma A + 2-adic valuation argument in the note


cells = {}
for alpha in range(4):
    for beta in range(4):
        sols = {(1,)}  # prime always
        # squarefree composite: (1-alpha)(2^r-2) = beta r, r>=2
        sqf_desc = None
        if alpha == 1 and beta == 0:
            sqf_desc = "ALL squarefree composites (every r>=2)"
        elif alpha == 1:
            sqf_desc = "none"
        elif alpha >= 2:
            sqf_desc = "none"
        else:  # alpha == 0: 2^r - 2 = beta r
            found = [r for r in range(2, 64) if (1 << r) - 2 == beta * r]
            for r in found:
                sols.add(tuple([1] * r))
            sqf_desc = "none" if not found else ", ".join(f"squarefree with r={r}" for r in found)
        bound = support_bound(alpha, beta)
        nonsq = []
        rmax = bound if bound is not None else 10
        for r in range(1, rmax + 1):
            T = target_T(alpha, beta, r)
            for fac in factorizations(T, r):
                if max(fac) >= 3:
                    prof = tuple(sorted((b - 1 for b in fac), reverse=True))
                    nonsq.append(prof)
        for prof in nonsq:
            sols.add(prof)
        cells[(alpha, beta)] = dict(sols=sols, sqf_desc=sqf_desc, bound=bound, nonsq=nonsq)

# finite hostile search: every profile with support<=10 and exponents<=12
from itertools import combinations_with_replacement

window_hits = {cell: set() for cell in cells}
nprof = 0
for r in range(1, 11):
    for comb in combinations_with_replacement(range(12, 0, -1), r):
        prof = tuple(comb)  # nonincreasing
        nprof += 1
        F, S, U = profile_F(prof), profile_S(prof), profile_U(prof)
        for (alpha, beta) in cells:
            if F == alpha * S + beta * U:
                window_hits[(alpha, beta)].add(prof)
print(f"window search: {nprof} profiles with support<=10, exponents<=12")


def predicted_in_window(alpha, beta, prof):
    """Membership in the PROVED solution set, evaluated on a window profile."""
    if prof == (1,):
        return True
    if all(a == 1 for a in prof):
        r = len(prof)
        if alpha == 1 and beta == 0:
            return True
        if alpha == 0:
            return (1 << r) - 2 == beta * r
        return False
    if (alpha, beta) == (2, 0):
        return prof == tuple([3] + [1] * (len(prof) - 1))
    return prof in cells[(alpha, beta)]["sols"]


for (alpha, beta), cell in cells.items():
    pred = {prof for prof in window_hits[(alpha, beta)]}
    for r in range(1, 11):
        for comb in combinations_with_replacement(range(12, 0, -1), r):
            prof = tuple(comb)
            require(predicted_in_window(alpha, beta, prof) == (prof in window_hits[(alpha, beta)]),
                    f"cell {(alpha,beta)} window mismatch at profile {prof}")
print("window search agrees with the proved classification in all 16 cells")

# direct integer control for all cells, N <= AUDIT_N
direct_sets = {cell: 0 for cell in cells}
for n in range(2, AUDIT_N + 1):
    prof = profile_of(factor_spf(n, spf_small))
    F, S, U = profile_F(prof), profile_S(prof), profile_U(prof)
    for (alpha, beta) in cells:
        hit = F == alpha * S + beta * U
        require(hit == predicted_in_window(alpha, beta, prof), f"integer control cell {(alpha,beta)} at N={n}")
        direct_sets[(alpha, beta)] += hit
print(f"direct integer control passed for all cells on 2<=N<={AUDIT_N}")


def fmt_prof(prof):
    letters = "pqrstuvwxyz"
    return "".join(f"{letters[i]}^{a}" if a > 1 else letters[i] for i, a in enumerate(prof))


print()
print("cell (a,b) | finite? | r-bound | #profiles | nonsquarefree profiles (exponents) | squarefree composites | #N<=2e5")
print("-----------+---------+---------+-----------+------------------------------------+-----------------------+--------")
omega_injective = {}
cell_finite = {}
for (alpha, beta), cell in sorted(cells.items()):
    finite = cell["bound"] is not None and not (alpha == 1 and beta == 0)
    cell_finite[(alpha, beta)] = finite
    nonsq_sorted = sorted(set(cell["nonsq"]), key=lambda p: (len(p), p))
    profs = " ".join(fmt_prof(p) for p in nonsq_sorted) if nonsq_sorted else "-"
    count = ("inf" if not finite else str(len(cell["sols"])))
    # all solution profiles of the cell on the support window r<=10 (for (1,0) this adds the squarefree composites)
    sols_list = sorted(cell["sols"] | ({tuple([1] * r) for r in range(2, 11)} if (alpha, beta) == (1, 0) else set()),
                       key=lambda p: (len(p), p))
    omegas = [sum(p) for p in sols_list]
    omega_injective[(alpha, beta)] = len(set(omegas)) == len(omegas)
    print(f"({alpha},{beta})      | {'yes' if finite else 'NO ':<7} | {str(cell['bound']) if cell['bound'] is not None else 'inf':<7} | {count:<9} | {profs:<34} | {cell['sqf_desc']:<21} | {direct_sets[(alpha,beta)]}")

print()
print("Omega-injectivity (distinct Omega values across all solution profiles; infinite cells listed on the window r<=10):")
for cell, flag in sorted(omega_injective.items()):
    alpha, beta = cell
    sols_list = sorted(cells[cell]["sols"] | ({tuple([1] * r) for r in range(2, 11)} if cell == (1, 0) else set()),
                       key=lambda p: (len(p), p))
    tag = "finite" if cell_finite[cell] else "INFINITE"
    print(f"  {cell}: {'INJECTIVE' if flag else 'collides '} ({tag:<8}) Omega values {[sum(p) for p in sols_list]}")
# (2,0) is Omega-injective for ALL r: Omega(p^3 q_1...q_(r-1)) = r+2 determines r, hence the profile (PROVED in note)
require(omega_injective[(2, 0)], "(2,0) must be Omega-injective (Omega = r+2)")
require(not omega_injective[(1, 0)], "(1,0) collides (p^2 and pq share Omega=2)")
require(not omega_injective[(0, 2)] and omega_injective[(1, 1)], "prime-cube cells: (0,2) collides, (1,1) injective")
inj_cells = [c for c, f in omega_injective.items() if f and cell_finite[c]]
inj_all = [c for c, f in omega_injective.items() if f]
require(inj_cells == [(0, 0), (1, 1), (2, 1), (2, 3), (3, 0), (3, 1)], f"Omega-injective finite cells: {inj_cells}")
three_inj = [c for c in inj_cells if len(cells[c]["sols"]) == 3]
require(three_inj == [(1, 1), (2, 3)], f"three-profile Omega-injective cells: {three_inj}")
print(f"  cells with Omega-injective finite solution sets: {inj_cells}")
print(f"  all Omega-injective cells (the infinite cell (2,0) included; (2,0) has Omega = r+2, injective for every r): {inj_all}")
print(f"  among the prime-cube cells (alpha+beta=2): (0,2) collides [1,3,3,3]; (1,1) and (2,0) are injective; (1,1) is the unique FINITE injective one")
print(f"  Omega-injective finite cells with exactly three profiles: {three_inj}")

# prime-power solution p^(alpha+beta+1) in every cell with alpha+beta>=1
for (alpha, beta), cell in cells.items():
    if alpha + beta >= 1:
        require((alpha + beta + 1,) in cell["sols"], f"prime power p^(a+b+1) missing in {(alpha,beta)}")
    require(tuple(x for x in cell["sols"] if len(x) == 1 and x != (1,)) != () or alpha + beta == 0, "")
print("prime-power law: every cell with alpha+beta>=1 contains exactly the prime power p^(alpha+beta+1) (checked)")
for (alpha, beta), cell in cells.items():
    pp = sorted(x for x in cell["sols"] if len(x) == 1)
    require(pp == ([(1,)] if alpha + beta == 0 else [(1,), (alpha + beta + 1,)]), f"prime-power set in {(alpha,beta)}: {pp}")

# the (2,0) infinite family, explicit check up to r=10 from the T-factorization
for r in range(1, 11):
    fams = [tuple(sorted((b - 1 for b in f), reverse=True)) for f in factorizations(target_T(2, 0, r), r) if max(f) >= 3]
    require(fams == [tuple([3] + [1] * (r - 1))], f"(2,0) family at r={r}: {fams}")
print("(2,0): F = 2S has, for every support r<=10, exactly the profile (3,1^(r-1)) beyond the prime (PROVED for all r in note)")
# 2-adic bound check for alpha in {2,3}: no nonsquarefree solutions at r=8..10 (consistent with proved r<=7)
for alpha in (2, 3):
    for beta in range(4):
        if (alpha, beta) == (2, 0):
            continue
        for r in range(8, 13):
            fams = [f for f in factorizations(target_T(alpha, beta, r), r) if max(f) >= 3]
            require(fams == [], f"unexpected solution at alpha={alpha}, beta={beta}, r={r}")
print("alpha in {2,3}: no nonsquarefree solutions at supports 8..12 (consistent with the proved bound r<=7)")
# hand-check values quoted in the note: T(r) for (2,1) and (3,3), r=1..7; (1,3) at r=5; the bound T/2^r <= 4.5
T21 = [target_T(2, 1, r) for r in range(1, 8)]
T33 = [target_T(3, 3, r) for r in range(1, 8)]
require(T21 == [5, 10, 19, 36, 69, 134, 263] and T33 == [8, 17, 32, 59, 110, 209, 404], "T(r) hand-check lists")
require(target_T(1, 3, 5) == 48 == 3 * 2 ** 4, "(1,3) at r=5 is the tight bound T=48")
ratio_max = max(Fraction(target_T(a, b, r), 1 << r) for a in (2, 3) for b in range(4) for r in range(1, 40) if (a, b) != (2, 0))
require(ratio_max == Fraction(17, 4) == Fraction(target_T(3, 3, 2), 4) and ratio_max <= Fraction(9, 2),
        f"true max of T/2^r over alpha in {{2,3}}, (alpha,beta) != (2,0): {ratio_max} (note's bound 4.5)")
print(f"hand-check values: T(2,1;r) r=1..7 = {T21}; T(3,3;r) r=1..7 = {T33}; T(1,3;5) = {target_T(1,3,5)} = 3*2^4")
print(f"  T/2^r over alpha in {{2,3}}, (alpha,beta) != (2,0): true maximum {float(ratio_max)} at (3,3), r=2; the note's bound alpha + 3r/2^r <= 4.5 holds")
print(f"  [time {time.time()-T0:.1f}s]")


# ----------------------------------------------------------------------------
# 3. Lattice identity and the k-free refinement F = S_k + U (PROVED + FINITE-EXACT)
# ----------------------------------------------------------------------------

banner("3. Box = cube + directions + apex; k-free refinement F = S_k + U")


def box_points(prof):
    from itertools import product
    return list(product(*[range(a + 1) for a in prof]))


for prof in ((3,), (2, 1, 1)):
    box = box_points(prof)
    cube = [e for e in box if max(e) <= 1]
    Q = [e for e in box if max(e) >= 2]  # non-squarefree divisors (order filter)
    r = len(prof)
    require(len(box) == len(cube) + r + 1, f"box = cube + r + 1 fails at {prof}")
    require(len(Q) == r + 1, f"|Q| = r+1 fails at {prof}")
    print(f"profile {prof}: |box|={len(box)} = |cube|={len(cube)} + directions r={r} + apex 1;  non-squarefree divisors |Q|={len(Q)}")
# the canonical bijection for p^2qr: d -> d/p maps Q\{N} onto the atoms
N = 2 ** 2 * 3 * 5
divs = sorted(divisors_from_factors({2: 2, 3: 1, 5: 1}))
Q = [d for d in divs if any(d % (p * p) == 0 for p in (2, 3, 5))]
# atom map l -> p*lcm(p,l): 2 -> 4, 3 -> 12, 5 -> 20 (inverse: p^2 -> p, p^2 q -> q)
image = sorted((2 if d == 4 else d // 4) for d in Q if d != N)
require(image == [2, 3, 5], f"atom map l -> p*lcm(p,l) is not a bijection Q\\{{N}} <-> atoms: {image}")
print(f"N=60=2^2*3*5: non-squarefree divisors Q={Q}; Q minus N = image of the atoms under l -> p*lcm(p,l) (2->4, 3->12, 5->20); apex N=60")
upper_face = sorted(d // 2 for d in Q)
require(upper_face == [2, 6, 10, 30] and upper_face == sorted(2 * d for d in divisors_from_factors({3: 1, 5: 1})), "d -> d/p sends Q onto the upper cube face p*D(m)")
print(f"  hostile: d -> d/p sends Q onto the upper cube face p*D(m) = {upper_face}, not onto the atoms")
# general: |Q| for profile (2,1^(r-1)) is 2^(r-1); equals r+1 iff r=3
for r in range(1, 9):
    prof = tuple([2] + [1] * (r - 1))
    q = prod(a + 1 for a in prof) - (1 << r)
    require(q == 1 << (r - 1), "Q size for (2,1^(r-1))")
    require((q == r + 1) == (r == 3), "Q=r+1 iff r=3")
print("profile (2,1^(r-1)): |Q| = 2^(r-1) = (upper cube face through p) and 2^(r-1) = r+1 iff r=3 (checked r<=8)")


def sk_solutions_predicted(k):
    if k == 1:
        return {(1,), (2,), (1, 1)}
    out = {(1,), (k + 1,), (k, 1, 1)}
    if k >= 3:
        out.add((k, 2))
    return out


print()
print("k | predicted solution profiles of F = S_k + U | window check (r<=8, exp<=12) | direct N<=2e5")
print("--+---------------------------------------------+------------------------------+-------------")
for k in range(1, 7):
    pred = sk_solutions_predicted(k)
    hits = set()
    for r in range(1, 9):
        for comb in combinations_with_replacement(range(12, 0, -1), r):
            prof = tuple(comb)
            if profile_F(prof) == profile_Sk(prof, k) + profile_U(prof):
                hits.add(prof)
    require(hits == pred, f"S_{k} window classification: {sorted(hits)} vs {sorted(pred)}")
    cnt = 0
    for n in range(2, AUDIT_N + 1):
        prof = profile_of(factor_spf(n, spf_small))
        hit = profile_F(prof) == profile_Sk(prof, k) + profile_U(prof)
        require(hit == (prof in pred), f"S_{k} direct control at N={n}")
        cnt += hit
    print(f"{k} | {' '.join(fmt_prof(p) for p in sorted(pred, key=lambda p:(len(p),p))):<43} | {len(hits)} profiles, all predicted | {cnt}")
# a few explicit F, S_3, U values
for prof in ((4,), (3, 2), (3, 1, 1), (3,), (2, 1, 1)):
    print(f"  profile {fmt_prof(prof):<8} F={profile_F(prof):>3} S_3={profile_Sk(prof,3):>3} U={profile_U(prof)}  F-S_3-U={profile_F(prof)-profile_Sk(prof,3)-profile_U(prof)}")
print(f"  [time {time.time()-T0:.1f}s]")


# ----------------------------------------------------------------------------
# 4. Omega blindness: partitions per Omega value and sandwich endpoint census
# ----------------------------------------------------------------------------

banner("4. Omega classes forget the repeated prime; partitions per Omega; endpoint census")


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, first):
            yield (first,) + rest


expected_sol = {1: 1, 3: 1, 4: 1}
print("Omega k | #profiles p(k) | #with F=S+U | which | defects D over the class")
print("--------+----------------+-------------+-------+-------------------------")
for k in range(1, 9):
    parts = list(partitions(k))
    sols = [p for p in parts if profile_F(p) == profile_S(p) + profile_U(p)]
    require(len(sols) == expected_sol.get(k, 0), f"solution count at Omega={k}")
    defects = sorted(set(profile_F(p) - profile_S(p) - profile_U(p) for p in parts))
    print(f"{k:<7} | {len(parts):<14} | {len(sols):<11} | {' '.join(fmt_prof(p) for p in sols) or '-':<5} | {defects}")
# the two Omega=4, omega=2 profiles have different D: Omega and omega jointly do not determine D
d1 = profile_F((2, 2)) - profile_S((2, 2)) - profile_U((2, 2))
d2 = profile_F((3, 1)) - profile_S((3, 1)) - profile_U((3, 1))
require((d1, d2) == (2, 1), "D at p^2q^2, p^3q")
print(f"witness: p^2q^2 and p^3q share (Omega,omega)=(4,2) but D = {d1} vs {d2}; Omega=4 alone: D in {sorted(set(profile_F(p)-profile_S(p)-profile_U(p) for p in partitions(4)))}")

# sandwich endpoint census: centers 6k, 1<=k<=10^6, endpoints 6k-1 (left), 6k+1 (right)
K = 1_000_000
ks = np.arange(1, K + 1, dtype=np.int64)
left = 6 * ks - 1
right = 6 * ks + 1
require(int(right[-1]) <= X, "endpoints within sieve")
print()
print("endpoint shape census over centers 6k, 1<=k<=10^6 (FINITE-EXACT):")
print("shape   | left 6k-1 | right 6k+1")
print("--------+-----------+-----------")
for name, arr in (("prime", shape_prime), ("p^3", shape_cube), ("p^2qr", shape_p2qr)):
    print(f"{name:<7} | {int(arr[left].sum()):>9} | {int(arr[right].sum()):>10}")
# p^3 endpoints: p^3 = 6k-1 iff p = 5 mod 6 ; p^3 = 6k+1 iff p = 1 mod 6 (p>=5), since p^2 = 1 mod 6
ENDMAX = int(right[-1])   # largest endpoint 6K+1 = 6000001 (NOT the sieve bound X)
cubes_left = [int(p) ** 3 for p in primes if int(p) ** 3 <= ENDMAX and int(p) % 6 == 5]
cubes_right = [int(p) ** 3 for p in primes if int(p) ** 3 <= ENDMAX and int(p) % 6 == 1]
require(len(cubes_left) == 21 == int(shape_cube[left].sum()), f"p^3 left census: {len(cubes_left)} vs {int(shape_cube[left].sum())}")
require(len(cubes_right) == 19 == int(shape_cube[right].sum()), f"p^3 right census: {len(cubes_right)} vs {int(shape_cube[right].sum())}")
require(set(cubes_left) == set(int(n) for n in left[shape_cube[left]]), "left cube endpoints are exactly the p = 5 mod 6 cubes")
require(set(cubes_right) == set(int(n) for n in right[shape_cube[right]]), "right cube endpoints are exactly the p = 1 mod 6 cubes")
# side law p^3 = p (mod 6) for every prime p >= 5; the boundary primes 2, 3 give 8 = 2, 27 = 3 (mod 6), never endpoints
require(all(int(p) ** 3 % 6 == int(p) % 6 for p in primes if int(p) ** 3 <= ENDMAX), "p^3 = p mod 6 for all primes with p^3 <= 6K+1")
require(8 % 6 == 2 and 27 % 6 == 3 and all((c + 1) % 6 != 0 and (c - 1) % 6 != 0 for c in (8, 27)), "8 and 27 are never endpoints 6k-1 or 6k+1")
pmax_cube = max(int(p) for p in primes if int(p) ** 3 <= ENDMAX)
print(f"  p^3 endpoints (largest endpoint {ENDMAX}, largest cube prime p = {pmax_cube}): left <=> p = 5 mod 6 ({len(cubes_left)} of them, first {cubes_left[:4]}); right <=> p = 1 mod 6 ({len(cubes_right)}, first {cubes_right[:4]})")
print(f"  side law p^3 = p (mod 6) holds for every prime with p^3 <= {ENDMAX}; p = 2, 3 give 8 = 2, 27 = 3 (mod 6), never endpoints, so the iff is over all primes")
# both endpoints solutions of F=S+U (any shape): joint count
sol_ind = shape_prime | shape_cube | shape_p2qr
both = int((sol_ind[left] & sol_ind[right]).sum())
twin = int((shape_prime[left] & shape_prime[right]).sum())
print(f"  centers with both endpoints F=S+U solutions: {both}; of which twin primes: {twin}; mixed/other: {both-twin}")
print(f"  [time {time.time()-T0:.1f}s]")

import resource
peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 ** 2)   # bytes on macOS
print(f"  [time total {time.time()-T0:.1f}s, peak RSS {peak_mb:.0f} MB]")
print()
print("ALL CHECKS PASSED")
