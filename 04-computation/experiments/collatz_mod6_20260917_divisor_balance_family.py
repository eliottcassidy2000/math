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
flips = int(np.count_nonzero((diff[1:] > 0) != (diff[:-1] > 0)))
fac_first = factor_spf(first, spf_sieve(first))
print()
print(f"CROSSOVER (FINITE-EXACT): first n with #{{p^2qr<=n}} > #{{primes<=n}} is n = {first} = {fac_first}")
print(f"  counts there: primes={int(cum_prime[first])}, p^2qr={int(cum_p2qr[first])}")
print(f"  last n <= 10^7 with #primes >= #p^2qr: {last_prime_lead}; sign changes of the difference on [1,10^7]: {flips}")
print(f"  at 10^7 the lead is p^2qr - primes = {int(diff[X])}")
# lead history at the tail: strictly increasing lead from some point?
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
for (alpha, beta), cell in sorted(cells.items()):
    finite = cell["bound"] is not None and not (alpha == 1 and beta == 0)
    nonsq_sorted = sorted(set(cell["nonsq"]), key=lambda p: (len(p), p))
    profs = " ".join(fmt_prof(p) for p in nonsq_sorted) if nonsq_sorted else "-"
    count = ("inf" if not finite else str(len(cell["sols"])))
    sols_list = sorted(cell["sols"], key=lambda p: (len(p), p))
    omegas = [sum(p) for p in sols_list]
    omega_injective[(alpha, beta)] = finite and len(set(omegas)) == len(omegas)
    print(f"({alpha},{beta})      | {'yes' if finite else 'NO ':<7} | {str(cell['bound']) if cell['bound'] is not None else 'inf':<7} | {count:<9} | {profs:<34} | {cell['sqf_desc']:<21} | {direct_sets[(alpha,beta)]}")

print()
print("Omega-injectivity (distinct Omega values across all solution profiles, finite cells only):")
for cell, flag in sorted(omega_injective.items()):
    sols_list = sorted(cells[cell]["sols"], key=lambda p: (len(p), p))
    print(f"  {cell}: {'INJECTIVE' if flag else 'collides '} Omega values {[sum(p) for p in sols_list]}")
inj_cells = [c for c, f in omega_injective.items() if f]
print(f"  cells with Omega-injective finite solution sets: {inj_cells}")

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
cubes_left = [int(p) ** 3 for p in primes if int(p) ** 3 <= X and int(p) % 6 == 5]
cubes_right = [int(p) ** 3 for p in primes if int(p) ** 3 <= X and int(p) % 6 == 1]
require(len(cubes_left) == int(shape_cube[left].sum()), "p^3 left census")
require(len(cubes_right) == int(shape_cube[right].sum()), "p^3 right census")
print(f"  p^3 endpoints: left <=> p = 5 mod 6 ({len(cubes_left)} of them, first {cubes_left[:4]}); right <=> p = 1 mod 6 ({len(cubes_right)}, first {cubes_right[:4]})")
# both endpoints solutions of F=S+U (any shape): joint count
sol_ind = shape_prime | shape_cube | shape_p2qr
both = int((sol_ind[left] & sol_ind[right]).sum())
twin = int((shape_prime[left] & shape_prime[right]).sum())
print(f"  centers with both endpoints F=S+U solutions: {both}; of which twin primes: {twin}; mixed/other: {both-twin}")
print(f"  [time {time.time()-T0:.1f}s]")

print()
print("ALL CHECKS PASSED")
