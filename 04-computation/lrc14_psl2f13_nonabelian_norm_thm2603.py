#!/usr/bin/env python3
"""Exact PSL_2(F_13) seven-edge norm and target-transition audit.

Universe and conventions
------------------------
* Matrices lie in SL_2(F_13), with A and -A identified in PSL_2.
* Matrices act on column projective coordinates x in P^1(F_13).
* U^a acts as x -> x+a and fixes infinity.
* A chronological list E_0,...,E_6 therefore composes as E_6...E_0.

The script verifies both the positive finite-group signal and the hostile
transition consequence.  In particular it prints the actual THM-2602
consequence object: the target-residue twisted-return support, not merely
the centrality of an intermediate matrix product.
"""

from collections import Counter, deque


P = 13
INF = 13
I = (1, 0, 0, 1)
MINUS_I = (P - 1, 0, 0, P - 1)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mat(a, b, c, d):
    return (a % P, b % P, c % P, d % P)


def mul(A, B):
    a, b, c, d = A
    e, f, g, h = B
    return mat(a * e + b * g, a * f + b * h,
               c * e + d * g, c * f + d * h)


def inv(A):
    a, b, c, d = A
    require((a * d - b * c) % P == 1, "matrix is not in SL2")
    return mat(d, -b, -c, a)


def power(A, n):
    require(n >= 0, "negative exponent passed to power")
    out = I
    base = A
    while n:
        if n & 1:
            out = mul(out, base)
        base = mul(base, base)
        n //= 2
    return out


def neg(A):
    return tuple((-x) % P for x in A)


def psl(A):
    """Canonical representative of {A,-A}."""
    B = neg(A)
    return min(A, B)


def pmul(A, B):
    return psl(mul(A, B))


def pinv(A):
    return psl(inv(A))


def ppow(A, n):
    return psl(power(A, n))


def U(a):
    return mat(1, a, 0, 1)


def G(t):
    return mat(0, 1, -1, t)


W = mat(0, -1, 1, 0)


def bruhat(r, s):
    return mul(mul(U(r), W), U(s))


def trace(A):
    return (A[0] + A[3]) % P


def is_central(A):
    return A == I or A == MINUS_I


def projective_order(A):
    x = psl(A)
    y = psl(I)
    for n in range(1, 200):
        y = pmul(y, x)
        if y == psl(I):
            return n
    raise RuntimeError("projective order search exceeded group order bound")


def act(A, x):
    a, b, c, d = A
    if x == INF:
        if c == 0:
            return INF
        return a * pow(c, -1, P) % P
    denominator = (c * x + d) % P
    if denominator == 0:
        return INF
    return (a * x + b) * pow(denominator, -1, P) % P


def permutation(A):
    ans = tuple(act(A, x) for x in range(P + 1))
    require(sorted(ans) == list(range(P + 1)), "non-permutation action")
    return ans


def generated_group(generators):
    gens = [psl(A) for A in generators]
    gens += [pinv(A) for A in gens]
    identity = psl(I)
    seen = {identity}
    queue = deque([identity])
    while queue:
        x = queue.popleft()
        for s in gens:
            y = pmul(x, s)
            if y not in seen:
                seen.add(y)
                queue.append(y)
    return seen


def conjugate(A, B):
    return mul(mul(A, B), inv(A))


def product(matrices):
    out = I
    for A in matrices:
        out = mul(out, A)
    return out


def cycles_of(A):
    perm = permutation(A)
    seen = set()
    cycles = []
    for x in range(P + 1):
        if x in seen:
            continue
        cycle = []
        y = x
        while y not in seen:
            seen.add(y)
            cycle.append(y)
            y = perm[y]
        cycles.append(tuple(cycle))
    return tuple(cycles)


def fixed_point(A):
    fixed = [x for x in range(P + 1) if act(A, x) == x]
    require(len(fixed) == 1, "expected a unique projective fixed point")
    return fixed[0]


def cyclotomic_sum_is_zero(exponent_counts):
    """Test sum_e c_e zeta_13^e=0 over Q using Phi_13."""
    require(len(exponent_counts) == P, "wrong cyclotomic coefficient vector")
    return len(set(exponent_counts)) == 1


def central_norm_data(t, a, reverse=False):
    g = G(t)
    gi = inv(g)
    factors = [conjugate(power(g, k), U(a)) for k in range(7)]
    algebraic = list(reversed(factors)) if reverse else factors
    norm = product(algebraic)
    # Column-vector chronological order that realizes this algebraic product.
    chronological = list(reversed(algebraic))
    return factors, algebraic, chronological, norm


def path_leakage(chronological):
    leaked = []
    avoided = []
    for q in range(P):
        x = q
        hit_inf = False
        for A in chronological:
            x = act(A, x)
            hit_inf = hit_inf or x == INF
        require(x == q, "central norm did not return projective state")
        (leaked if hit_inf else avoided).append(q)
    return tuple(leaked), tuple(avoided)


def main():
    traces7_expected = {3, 5, 6, 7, 8, 10}  # +/-{3,5,6}
    success_expected = {
        (3, "+"): (6, 8, 9, 10, 11),
        (3, "-"): (2, 3, 4, 5, 7),
        (5, "+"): (2, 8, 10, 11, 12),
        (5, "-"): (1, 2, 3, 5, 11),
        (6, "+"): (1, 3, 9, 11, 12),
        (6, "-"): (1, 2, 4, 10, 12),
    }

    groups = {}
    for t in (3, 5, 6):
        g = G(t)
        require(power(g, 7) == MINUS_I, f"g_{t}^7 != -I")
        require(projective_order(g) == 7, f"g_{t} projective order")
        groups[t] = generated_group([U(1), g])
        require(len(groups[t]) == 1092, f"<U,g_{t}> is not PSL2(13)")
    group = groups[3]
    require(groups[3] == groups[5] == groups[6], "generator pairs differ")

    # The natural action is faithful: it has 1,092 distinct permutations.
    action = {permutation(A) for A in group}
    require(len(action) == 1092, "projective-line action is not faithful")

    trace_order7 = set()
    order_histogram = Counter()
    for A in group:
        o = projective_order(A)
        order_histogram[o] += 1
        if o == 7:
            # A canonical PSL representative chooses one of the two SL traces;
            # include both signs because trace is defined only up to sign in PSL.
            trace_order7.add(trace(A))
            trace_order7.add((-trace(A)) % P)
    require(trace_order7 == traces7_expected, "wrong order-seven trace locus")
    require(order_histogram[7] == 468, "wrong number of order-seven elements")

    # Normalizer of the target translation deck.
    deck = {psl(U(a)) for a in range(P)}
    normalizer = set()
    for h in group:
        image = {pmul(pmul(h, x), pinv(h)) for x in deck}
        if image == deck:
            normalizer.add(h)
    require(len(normalizer) == 78, "wrong Sylow-13 normalizer size")
    normalizer_orders = sorted({projective_order(h) for h in normalizer})
    require(7 not in normalizer_orders, "order-seven target-deck normalizer found")

    success_rows = {}
    leakage_rows = {}
    successful_factor_generated_orders = set()
    all_successful_a = set()
    for t in (3, 5, 6):
        g = G(t)
        require(tuple(sorted(len(c) for c in cycles_of(g))) == (7, 7),
                "g does not have two seven-cycles on P1")
        fixed_orbit = tuple(act(power(g, k), INF) for k in range(7))
        require(len(set(fixed_orbit)) == 7, "conjugate decks repeat")

        conjugate_decks = []
        for k in range(7):
            Vk = conjugate(power(g, k), U(1))
            Pk = {psl(power(Vk, a)) for a in range(P)}
            conjugate_decks.append(Pk)
            require(fixed_point(Vk) == fixed_orbit[k], "wrong deck fixed point")
        require(len({frozenset(x) for x in conjugate_decks}) == 7,
                "fewer than seven target decks")
        for i in range(7):
            for j in range(i):
                require(conjugate_decks[i] & conjugate_decks[j] == {psl(I)},
                        "distinct prime decks intersect nontrivially")

        for reverse, sign in ((False, "+"), (True, "-")):
            successes = []
            for a in range(1, P):
                factors, algebraic, chronological, norm = central_norm_data(
                    t, a, reverse=reverse)
                moving = mul(U(a), g) if not reverse else mul(inv(g), U(a))
                moving_trace = (t - a) % P if not reverse else (t + a) % P
                require(trace(moving) == moving_trace, "moving trace sign error")

                # Telescope identities in the stated algebraic convention.
                if not reverse:
                    rhs = mul(power(mul(U(a), g), 7), power(inv(g), 7))
                else:
                    rhs = mul(power(g, 7), power(mul(inv(g), U(a)), 7))
                require(norm == rhs, "nonabelian norm telescope failed")

                central_by_trace = moving_trace in traces7_expected
                require(is_central(norm) == central_by_trace,
                        "trace criterion and norm centrality disagree")
                if not is_central(norm):
                    continue
                successes.append(a)
                all_successful_a.add(a)
                factor_generated_order = len(generated_group(factors))
                successful_factor_generated_orders.add(factor_generated_order)
                require(factor_generated_order == 1092,
                        "successful factor bank generates a proper subgroup")
                leaked, avoided = path_leakage(chronological)
                require(3 <= len(leaked) <= 5, "unexpected infinity leakage")

                # Full P1 norm is identity in PSL.  Twisting by the base target
                # holonomy U^(7a) fixes only infinity; the affine target trace is 0.
                twist = U(7 * a)
                full_returns = tuple(
                    x for x in range(P + 1)
                    if act(twist, act(norm, x)) == x
                )
                target_returns = tuple(x for x in full_returns if x != INF)
                require(full_returns == (INF,), "twisted full return is not cemetery-only")
                require(target_returns == (), "spurious target-residue twisted return")

                # Restrict every chronological factor to affine target states.
                # Its product is identity on precisely the paths avoiding infinity;
                # after the nonzero target twist it has no diagonal support.
                partial_product = {}
                for q in range(P):
                    x = q
                    alive = True
                    for A in chronological:
                        x = act(A, x)
                        if x == INF:
                            alive = False
                            break
                    if alive:
                        partial_product[q] = x
                require(set(partial_product) == set(avoided),
                        "affine partial-product support mismatch")
                require(all(q == x for q, x in partial_product.items()),
                        "affine restricted norm is not diagonal")
                restricted_twisted = [
                    q for q, x in partial_product.items()
                    if act(twist, x) == q
                ]
                require(restricted_twisted == [],
                        "affine restricted norm passes twisted-return test")
                leakage_rows[(t, sign, a)] = (len(leaked), len(avoided))

            success_rows[(t, sign)] = tuple(successes)
            require(tuple(successes) == success_expected[(t, sign)],
                    "central norm success atlas mismatch")
    require(all_successful_a == set(range(1, P)), "trace atlas misses a target root")
    require(successful_factor_generated_orders == {1092},
            "successful factor-bank generated orders vary")

    # The 13 x 13 rank-one Bruhat square is the exact algebraic analogue of
    # an independent source/future action carrier.  The group right coordinate
    # s_B has the opposite sign from THM-2615's future coordinate s=-s_B.
    # Thus using exponent -lambda*r-nu*s_B below is exactly the typed transform
    # -lambda*r+nu*s.  The order-seven wall depends only on r+s_B=r-s, so its
    # support lies on lambda=nu and contributes only q=lambda-nu=0.
    bruhat_square = {(r, s): psl(bruhat(r, s))
                     for r in range(P) for s in range(P)}
    require(len(set(bruhat_square.values())) == P * P,
            "Bruhat coordinate square is not injective")
    order7_wall = set()
    for r in range(P):
        for s in range(P):
            B = bruhat(r, s)
            require(B == mat(r, r * s - 1, 1, s), "Bruhat formula failed")
            require(mul(mul(U(4), B), U(9)) == bruhat(r + 4, s + 9),
                    "independent left/right translation failed")
            on_wall = projective_order(B) == 7
            require(on_wall == ((r + s) % P in traces7_expected),
                    "Bruhat trace wall mismatch")
            if on_wall:
                order7_wall.add((r, s))
    require(len(order7_wall) == 78, "wrong Bruhat order-seven wall size")
    for t in (3, 5, 6):
        require(psl(G(t)) == psl(bruhat(0, -t)), "g Bruhat coordinate mismatch")
        for a in range(P):
            require(psl(mul(U(a), G(t))) == psl(bruhat(a, -t)),
                    "forward moving matrix Bruhat mismatch")
            require(psl(mul(inv(G(t)), U(a))) == psl(bruhat(t, a)),
                    "reverse moving matrix Bruhat mismatch")

    bruhat_fourier_support = set()
    for lam in range(P):
        for nu in range(P):
            counts = [0] * P
            for r, s in order7_wall:
                counts[(-lam * r - nu * s) % P] += 1
            if not cyclotomic_sum_is_zero(counts):
                bruhat_fourier_support.add((lam, nu))
    expected_fourier_support = {(c, c) for c in range(P)}
    require(bruhat_fourier_support == expected_fourier_support,
            "Bruhat wall Fourier support is not the diagonal")
    require({(lam - nu) % P for lam, nu in bruhat_fourier_support} == {0},
            "Bruhat wall creates a primitive target-difference mode")

    # Lawful adjacent chart transport.  Q_k=P1\{g^k infinity}; V_k acts as
    # the target C13 translation in chart k.  The complete set of equivariant
    # chart maps is T_k(c)=V_(k+1)^c g=g V_k^c.
    torsor_sum_histogram = Counter({0: 1})
    for _ in range(7):
        next_histogram = Counter()
        for old_sum, count in torsor_sum_histogram.items():
            for c in range(P):
                next_histogram[(old_sum + c) % P] += count
        torsor_sum_histogram = next_histogram
    require(set(torsor_sum_histogram.values()) == {P ** 6},
            "equivariant connection sums are not uniform")

    for t in (3, 5, 6):
        g = G(t)
        V = [conjugate(power(g, k), U(1)) for k in range(8)]
        fixed = [fixed_point(Vk) for Vk in V]
        require(fixed[7] == fixed[0], "seven-chart target bundle does not close")
        for k in range(7):
            Qk = set(range(P + 1)) - {fixed[k]}
            Qnext = set(range(P + 1)) - {fixed[k + 1]}
            for c in range(P):
                Tk = mul(power(V[k + 1], c), g)
                require(Tk == mul(g, power(V[k], c)), "intertwiner identity failed")
                require({act(Tk, x) for x in Qk} == Qnext,
                        "intertwiner misses adjacent target fibre")
                require(psl(mul(Tk, V[k])) == psl(mul(V[k + 1], Tk)),
                        "target-deck equivariance failed")

        for a in range(1, P):
            # Natural connection c_k=0.  Edge E_k=g V_k^a is locally U^a.
            edges = [mul(g, power(V[k], a)) for k in range(7)]
            chronological_product = product(reversed(edges))
            require(psl(chronological_product) == psl(U(7 * a)),
                    "lawful natural transport did not restore U^(7a)")

            # One explicit connection invoice; the general law depends only on
            # the sum of c_k, and each sum has 13^6 witnesses.
            cs = [(-7 * a) % P] + [0] * 6
            corrected_edges = [
                mul(mul(power(V[k + 1], cs[k]), g), power(V[k], a))
                for k in range(7)
            ]
            corrected_product = product(reversed(corrected_edges))
            require(psl(corrected_product) == psl(I),
                    "external C13 connection invoice failed")

            # A second nontrivial tuple checks the full sum formula.
            cs2 = [(a + 2 * k + t) % P for k in range(7)]
            edges2 = [
                mul(mul(power(V[k + 1], cs2[k]), g), power(V[k], a))
                for k in range(7)
            ]
            product2 = product(reversed(edges2))
            require(psl(product2) == psl(U(7 * a + sum(cs2))),
                    "C13 connection holonomy formula failed")

    print("PSL2(F13) seven-edge nonabelian norm / transition audit")
    print(f"group_order={len(group)} faithful_P1_degree=14")
    print(f"projective_order_histogram={sorted(order_histogram.items())}")
    print(f"order7_signed_trace_set={sorted(trace_order7)}")
    print(f"target_deck_normalizer_size={len(normalizer)} orders={normalizer_orders}")
    for key in sorted(success_rows):
        print(f"central_norm[{key[0]},{key[1]}]={list(success_rows[key])}")
    leak_hist = Counter(leak for leak, avoid in leakage_rows.values())
    avoid_hist = Counter(avoid for leak, avoid in leakage_rows.values())
    print(f"successful_norms={len(leakage_rows)} leak_count_hist={sorted(leak_hist.items())}")
    print(f"successful_factor_bank_generated_orders={sorted(successful_factor_generated_orders)}")
    print(f"affine_avoid_count_hist={sorted(avoid_hist.items())}")
    print("twisted_return_full_P1={infinity}; twisted_return_target_F13=empty")
    print("bruhat_square=169 order7_wall=78 Fourier_support=lambda=nu only")
    print("bruhat_target_difference_support={0}")
    print("lawful_natural_chart_holonomy=U^(7a), not identity")
    print(f"equivariant_connection_choices_per_total_holonomy={P ** 6}")
    print("cancellation_condition=sum(c_k)=-7a mod 13")
    print("PASS")


if __name__ == "__main__":
    main()
