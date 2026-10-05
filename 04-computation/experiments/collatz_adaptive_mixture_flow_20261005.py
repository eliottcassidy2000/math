#!/usr/bin/env python3
"""Exact checks of adaptive rooted flows; no universal convergence assumption.

Only finite orbit checks use discovery. The theorem about all rooted sources
is proved separately in the companion note. Checks survive python -O.
"""
from fractions import Fraction as F
from collections import deque
from functools import lru_cache
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json

CHECKS = 0


def require(condition, witness=None):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(witness)


def v2(n):
    if n <= 0:
        raise ValueError(n)
    return (n & -n).bit_length() - 1


def U(n):
    t = 3 * n + 1
    return t >> v2(t)


def sibling(n):
    return 4 * n + 1


def unsibling(n):
    j = 0
    while n % 8 == 5:
        n = (n - 1) // 4
        j += 1
    return n, j


def inverse_base(y):
    if y % 3 == 0:
        raise ValueError(y)
    a = 1 if y % 3 == 2 else 2
    return ((1 << a) * y - 1) // 3


@lru_cache(None)
def weight(L, K):
    return F(2 * factorial(K) * factorial(L + 1), factorial(L + K + 2))


def beta_integer(a, b):
    return F(factorial(a - 1) * factorial(b - 1), factorial(a + b - 1))


def singular_weight(L, K):
    if K < 1:
        raise ValueError("The singular mixture excludes ROOT and requires K>=1")
    return beta_integer(K, L + 2)


def cost_kernel(q):
    """Discover then verify every rooted base of total sibling cost <=q.

    No depth cutoff is assumed to be complete: exceeding a resource cap is
    an error. The forward decoder and exhaustive boundary check are separate.
    """
    entries = {1: (0, 0, None, None)}
    todo = deque([1])
    while todo:
        c = todo.popleft()
        L, K, _, _ = entries[c]
        y = c
        for k in range(q - K + 1):
            if y % 3 and not (c == 1 and k == 0):
                b = inverse_base(y)
                require(b not in entries, ("duplicate inverse address", b))
                entries[b] = (L + 1, K + k, c, k)
                todo.append(b)
            y = sibling(y)
        require(len(entries) <= 200000 and L <= 2000, "kernel search cap, no completion claim")

    verify_cost_kernel(entries, q)
    return entries


def verify_cost_kernel(entries, q):
    require(entries.get(1) == (0, 0, None, None), "bad ROOT boundary")
    for b, (L, K, parent, depth) in entries.items():
        require(type(b) is int and b > 0 and b % 2 == 1)
        require(type(L) is int and type(K) is int and L >= 0 and 0 <= K <= q)
        require(unsibling(b) == (b, 0) and K <= q)
        if b != 1:
            require(parent in entries and type(depth) is int and depth >= 0)
            require(unsibling(U(b)) == (parent, depth))
            lp, kp, _, _ = entries[parent]
            require(L == lp + 1 and K == kp + depth)
        # Independent closed-form fibre, rather than repeated S in discovery.
        for k in range(q - K + 1):
            y = 4 ** k * b + (4 ** k - 1) // 3
            if y % 3 and not (b == 1 and k == 0):
                a = 1 if y % 6 == 5 else 2
                child = (2 ** a * y - 1) // 3
                require(child in entries and entries[child] == (L + 1, K + k, b, k))


def require_rejected(entries, q):
    try:
        verify_cost_kernel(entries, q)
    except RuntimeError:
        require(True)
    else:
        require(False, "invalid kernel was accepted")


def polynomial_integral(L, K):
    # Independent integration: expand (1-r)^(L+1), integrate each monomial.
    return 2 * sum((F((-1) ** i * comb(L + 1, i), K + i + 1)
                    for i in range(L + 2)), F(0))


base_cache = {1: (0, 0)}


def labels(n):
    b, j = unsibling(n)
    c, path, seen = b, [], set()
    while c not in base_cache:
        require(c not in seen and len(path) < 10000, ("unresolved finite input", n, c))
        seen.add(c)
        parent, k = unsibling(U(c))
        path.append((c, parent, k))
        c = parent
    for c, parent, k in reversed(path):
        L, K = base_cache[parent]
        base_cache[c] = L + 1, K + k
    L, K = base_cache[b]
    return L, K + j


def direct_labels(n):
    # Separate reader: actual odd orbit, without sibling/base reduction.
    if n == 1:
        return (0, 0), 0
    tau = K = 0
    seen = set()
    while n != 1:
        require(n not in seen and tau < 10000, ("unresolved direct input", n))
        seen.add(n)
        a = v2(3 * n + 1)
        K += (a - 1) // 2
        n = (3 * n + 1) >> a
        tau += 1
    return (tau - 1, K), tau


def rational_record(x):
    return {"numerator": x.numerator, "denominator": x.denominator}


def main():
    fixed_bounds = []
    for r in (F(1, 64), F(1, 16), F(1, 4), F(1, 2), F(3, 4), F(15, 16)):
        D = 1 + r + r * r
        kappa = 1 - r * r / D
        first = r * (1 + r * r) / D
        child3 = (1 - r) * r
        row3 = (r + r * r) / D
        second_bound = child3 * row3 + (first - child3) * kappa
        base_bound = 1 + first + second_bound / (1 - kappa)
        require(base_bound == 2 + 2 * r - r * r)
        # Independently enumerate the first root generation and its row sums.
        first_finite = second_finite = F(0)
        target = 1
        for k in range(25):
            if target % 3 and k:
                child = inverse_base(target)
                require(unsibling(child)[1] == 0 and U(child) == target)
                d = (1 - r) * r ** k
                first_finite += d
                forbidden = (-child) % 3
                second_finite += d * (1 - r ** forbidden / D)
                if k == 1:
                    require(child == 3 and forbidden == 0)
            target = sibling(target)
        require(first_finite < first and second_finite <= second_bound)
        fixed_bounds.append({"r": str(r), "base_bound": str(base_bound),
                             "full_bound": str(base_bound / (1 - r))})

    # A two-counter array: exact integration, Pascal splitting, complete tails.
    for L in range(61):
        for K in range(61):
            w = weight(L, K)
            require(w > 0 and w == 2 * beta_integer(K + 1, L + 2))
            require(w == weight(L + 1, K) + weight(L, K + 1))
            require(weight(L + 1, K) / w == F(L + 2, L + K + 3))
            require(weight(L, K + 1) / w == F(K + 1, L + K + 3))
            require(sum((weight(L + 1, K + j) for j in range(8)), F(0))
                    + weight(L, K + 8) == w)
            if L <= 18 and K <= 18:
                require(w == polynomial_integral(L, K))
            if K:
                z = singular_weight(L, K)
                require(z == w * F(L + K + 2, 2 * K))
                require(z == singular_weight(L + 1, K) + singular_weight(L, K + 1))
            t = L + K
            best = (F(L, t) ** L * F(K, t) ** K) if t else F(1)
            require(F(2 * (L + 1), (t + 1) * (t + 2)) * best <= w <= best)
    require(weight(0, 0) == 1)
    require(sum((weight(0, j) for j in range(100)), F(0)) + F(2, 101) == 2)

    # General beta-prior bound: beta>1 is needed for the infinite root fibre.
    beta_bounds = []
    for alpha in range(1, 6):
        for beta in range(2, 7):
            z = beta_integer(alpha, beta)
            integral = (2 * beta_integer(alpha, beta - 1)
                        + 2 * beta_integer(alpha + 1, beta - 1)
                        - beta_integer(alpha + 2, beta - 1)) / z
            formula = F(3 * (alpha + beta - 1), beta - 1) - F(beta, alpha + beta)
            require(integral == formula)
            if (alpha, beta) in ((1, 2), (2, 2), (1, 3)):
                beta_bounds.append({"alpha": alpha, "beta": beta, "bound": str(formula)})
    require(beta_bounds[0]["bound"] == "16/3")

    # Exhaust the total-cost budget, retaining arbitrary base depth. This
    # finite completed certificate controls the entire unbounded remainder.
    largest_kernel = cost_kernel(8)
    small_kernel = {b: row for b, row in largest_kernel.items() if row[1] <= 2}
    missing_zero_cost_child = dict(small_kernel)
    del missing_zero_cost_child[9]
    require_rejected(missing_zero_cost_child, 2)
    bad_root = dict(small_kernel)
    bad_root[1] = (1, 0, None, None)
    require_rejected(bad_root, 2)
    bad_bill = dict(small_kernel)
    bad_bill[3] = (1, 0, 1, 1)
    require_rejected(bad_bill, 2)
    kernel_bounds = []
    previous_proper = previous_singular = None
    for q in range(2, 9):
        kernel = {b: row for b, row in largest_kernel.items() if row[1] <= q}
        P_integral = sum((beta_integer(K + 1, L + 1)
                          for L, K, _, _ in kernel.values()), F(0))
        proper = 2 * P_integral
        singular = 1 + sum((beta_integer(K, L + 1)
                            for b, (L, K, _, _) in kernel.items() if b != 1), F(0))
        for c, (L, K, _, _) in kernel.items():
            h = q - K + 1
            e = (-c - h) % 3
            for j in range(3):
                if j != e:
                    proper += 2 * beta_integer(q + j, L + 1)
                    singular += beta_integer(q - 1 + j, L + 1)
        for r in (F(1, 16), F(1, 4), F(1, 2), F(3, 4)):
            D = 1 + r + r * r
            P = sum(((1 - r) ** L * r ** K for L, K, _, _ in kernel.values()), F(0))
            boundary = F(0)
            boundary_numerator = F(0)
            for c, (L, K, _, _) in kernel.items():
                h = q - K + 1
                e = (-c - h) % 3
                boundary += (1 - r) ** L * r ** (K + h) * (1 - r ** e / D)
                boundary_numerator += (1 - r) ** L * (D - r ** e)
            require(boundary * D == r ** (q + 1) * boundary_numerator)
            R = P + boundary * D / (r * r)
            require(R <= 2 + 2 * r - r * r)
            if q == 2:
                polynomial = 1 + 5*r - 3*r**2 - 4*r**3 + 10*r**4 - 10*r**5 + 5*r**6 - r**7
                require(R == polynomial)
        if q == 2:
            require(set(kernel) == {1, 3, 17, 11, 7, 9})
            require(proper == F(407, 84) and singular == F(61, 14))
        if previous_proper is not None:
            require(proper < previous_proper and singular < previous_singular)
        previous_proper, previous_singular = proper, singular
        kernel_bounds.append({"budget": q, "base_count": len(kernel),
                              "max_base_depth": max(row[0] for row in kernel.values()),
                              "proper_full_mass_bound": rational_record(proper),
                              "singular_nonroot_mass_bound": rational_record(singular)})
    require(len(largest_kernel) == 3591)
    require(previous_proper < F(41, 10) and previous_singular < F(23, 8))
    require(sum((singular_weight(0, j) for j in range(1, 101)), F(0)) + F(1, 101) == 1)

    rows = []
    total_head = flux_head = singular_head = singular_flux_head = F(0)
    maxL = maxK = 0
    limit = 1 << 15
    for n in range(1, limit, 2):
        L, K = labels(n)
        independent, tau = direct_labels(n)
        require((L, K) == independent, n)
        require(n == 1 or K >= 1, n)
        maxL, maxK = max(maxL, L), max(maxK, K)
        w = weight(L, K)
        total_head += w
        if n % 3 == 0:
            flux_head += w
        if n > 1:
            singular_head += singular_weight(L, K)
            if n % 3 == 0:
                singular_flux_head += singular_weight(L, K)
        if n > 1 and U(n) != 1:
            lm, km = labels(U(n))
            require(L == lm + 1 and K == km + (v2(3 * n + 1) - 1) // 2)
        if n < 4096 and n > 1 and n % 3:
            b = inverse_base(n)
            require(labels(b) == (L + 1, K), (n, b))
            p = b
            partial = F(0)
            for j in range(6):
                require(U(p) == n and labels(p) == (L + 1, K + j))
                partial += weight(*labels(p))
                p = sibling(p)
            require(partial + weight(L, K + 6) == w)
            require(sum((singular_weight(L + 1, K + j) for j in range(6)), F(0))
                    + singular_weight(L, K + 6) == singular_weight(L, K))
        if n in (1, 3, 5, 7, 9, 13, 17, 21, 27, 53, 113, 155, 111):
            rows.append({"n": n, "L": L, "K": K, "odd_steps": tau,
                         "weight": str(w)})
    require(total_head < F(16, 3) and flux_head < 1)
    require(singular_head < F(23, 8) and singular_flux_head < 1)

    # Transport through the existing six-edge common-future paid controller.
    word = (1, 2, 1, 1, 1, 2)
    controller_rows = []
    for t in range(128):
        n, h = 155 + 2048 * t, 111 + 1458 * t
        x = n
        for a in word:
            require(v2(3 * x + 1) == a and x > 1)
            x = (3 * x + 1) >> a
        require(x == sibling(h) and h < n)
        l, k = labels(h)
        require(labels(n) == (l + 6, k + 1))
        require(weight(*labels(n)) < weight(l, k))
        if t < 3:
            controller_rows.append({"source": n, "child": h,
                                    "child_labels": [l, k],
                                    "ratio": str(weight(l + 6, k + 1) / weight(l, k))})

    # A bounded 3-divisible inverse section; every unit target is reached.
    section_exponents = {1: 6, 2: 5, 4: 4, 5: 1, 7: 2, 8: 3}
    for n in range(3, 8192, 2):
        if n % 3 == 0:
            continue
        a = section_exponents[n % 9]
        z = ((1 << a) * n - 1) // 3
        require(z % 3 == 0 and z % 2 == 1 and 0 < z < 22 * n and U(z) == n)
        L, K = labels(n)
        j = (a - 1) // 2
        require(labels(z) == (L + 1, K + j))
        num = L + 2
        den = 1
        for i in range(j):
            num *= K + 1 + i
        for i in range(j + 1):
            den *= L + K + 3 + i
        require(weight(*labels(z)) / weight(L, K) == F(num, den))

    # Hostile controls: counts do not identify a legal route or its source.
    require(labels(53) == labels(113) == (1, 3))
    require(v2(3 * 53 + 1) == 5 and v2(3 * 113 + 1) == 2)
    # A ROOT-anchored cost classification cannot be restarted mid-orbit:
    # arbitrary long legitimate G-prefixes have zero sibling cost.
    R = 64
    for i in range(R):
        n = 3 ** i * 2 ** (R + 3 - i) - 1
        m = 3 ** (i + 1) * 2 ** (R + 2 - i) - 1
        require(v2(3 * n + 1) == 1 and U(n) == m > n > 1)
        require(unsibling(n) == (n, 0) and unsibling(m) == (m, 0))
    # Sending counter mass to an unrelated finite-grid boundary is not ROOT.
    require(weight(1, 0) + weight(0, 1) == 1)
    # Exact weight 1 at ROOT is NOT the target's killed incoming inequality.
    require(U(1) == 1)
    # Prior uniform in r has infinite total mass just on the rooted sibling ray:
    # integral r^j dr = 1/(j+1); finite partial sums exceed any fixed bound.
    harmonic_head = sum((F(1, j + 1) for j in range(1000)), F(0))
    require(harmonic_head > 7)

    out = {
        "status": "PROVED scoped identities; FINITE-EXACT checks; universal positivity OPEN",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": CHECKS,
        "universe": {"odd_sources": [1, limit - 1], "counter_rectangle": [60, 60],
                     "inverse_fibre_targets_below": 4096, "tested_siblings_per_target": 6,
                     "controller_parameters": [0, 127], "section_targets_below": 8192},
        "fixed_parameter_mass_bounds": fixed_bounds,
        "beta_prior_bounds": beta_bounds,
        "adaptive_mass_upper_bound": "16/3",
        "root_sibling_mass": "2",
        "three_divisible_flux": "1",
        "singular_nonroot_mass_upper_bound": "23/8 (strict; verified kernel)",
        "improved_proper_full_mass_upper_bound": "41/10 (strict; verified kernel)",
        "cost_kernel_bounds": kernel_bounds,
        "cost_kernel_columns": ["base", "L", "K", "parent", "edge_depth"],
        "cost_kernel": [[b, *row] for b, row in sorted(largest_kernel.items())],
        "finite_head_mass": rational_record(total_head),
        "finite_head_three_divisible_flux": rational_record(flux_head),
        "finite_head_singular_mass": rational_record(singular_head),
        "finite_head_singular_flux": rational_record(singular_flux_head),
        "max_head_labels": [maxL, maxK],
        "sample_weights": rows,
        "controller_transports": controller_rows,
        "controls": {"same_counts_different_words": [53, 113],
                     "nonroot_zero_cost_prefix_length": R,
                     "uniform_prior_root_fibre_diverges": True,
                     "partial_counters_do_not_certify_root": True},
        "limits": "Finite checks do not establish support positivity at untested integers."
    }
    dest = Path(__file__).resolve().parents[2] / "05-knowledge/results/collatz_adaptive_mixture_flow_20261005.json"
    kernel_rows = out.pop("cost_kernel")
    encoded = json.dumps(out, indent=2, sort_keys=True)
    encoded = (encoded[:-2] + ',\n  "cost_kernel": [\n'
               + ',\n'.join('    ' + json.dumps(row) for row in kernel_rows)
               + '\n  ]\n}\n')
    dest.write_text(encoded)
    print(json.dumps({"checks": CHECKS, "head_mass_decimal": float(total_head),
                      "head_flux_decimal": float(flux_head),
                      "singular_head_mass_decimal": float(singular_head),
                      "singular_head_flux_decimal": float(singular_flux_head),
                      "output": str(dest)}, sort_keys=True))


if __name__ == "__main__":
    main()
