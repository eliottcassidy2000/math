#!/usr/bin/env python3
"""HYP-2982 analytic sieve packet weights for LRC14.

This is a finite atlas, not an LRC proof.  It makes explicit the arithmetic
weights mentioned in the prompt and asks which ones are safe proof carriers.

The main computed channels are:

* M(x) = sum_{n<=x} mu(n)
* A(x) = sum_{n<=x} mu(n)/n
* G(x) = sum_{n<=x} mu(n)^2 / phi(n)
* Phi(x) = sum_{q<=x} phi(q)
* prime sums theta(x) and sum_{p<=x} 1/p

The transfer lesson is then written as a Tournament Analysis on proof carriers:
vertices are analytic packets and LRC proof obligations, not runners.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from math import gcd, log, sqrt


N = 200_000
CHECKPOINTS = [
    14,
    25,
    27,
    36,
    41,
    63,
    64,
    84,
    98,
    128,
    168,
    256,
    280,
    512,
    1024,
    2048,
    4096,
    4312,
    8192,
    16384,
    32768,
    65536,
    131072,
    200000,
]
LRC_QS = [14, 25, 27, 36, 41, 63, 84, 98, 168, 280, 4312]


def linear_mu_phi(limit: int):
    primes: list[int] = []
    is_comp = [False] * (limit + 1)
    mu = [0] * (limit + 1)
    phi = [0] * (limit + 1)
    mu[1] = 1
    phi[1] = 1
    for i in range(2, limit + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
            phi[i] = i - 1
        for p in primes:
            v = i * p
            if v > limit:
                break
            is_comp[v] = True
            if i % p == 0:
                mu[v] = 0
                phi[v] = phi[i] * p
                break
            mu[v] = -mu[i]
            phi[v] = phi[i] * (p - 1)
    return primes, mu, phi


def prime_factorization(n: int, primes: list[int]) -> dict[int, int]:
    out: dict[int, int] = {}
    m = n
    for p in primes:
        if p * p > m:
            break
        while m % p == 0:
            out[p] = out.get(p, 0) + 1
            m //= p
    if m > 1:
        out[m] = out.get(m, 0) + 1
    return out


def prefix_tables(primes: list[int], mu: list[int], phi: list[int]):
    is_prime = [False] * (N + 1)
    for p in primes:
        is_prime[p] = True

    M = [0] * (N + 1)
    A = [0.0] * (N + 1)
    G = [0.0] * (N + 1)
    Phi = [0] * (N + 1)
    theta = [0.0] * (N + 1)
    recip_p = [0.0] * (N + 1)

    for n in range(1, N + 1):
        M[n] = M[n - 1] + mu[n]
        A[n] = A[n - 1] + mu[n] / n
        G[n] = G[n - 1]
        if mu[n] != 0:
            G[n] += 1.0 / phi[n]
        Phi[n] = Phi[n - 1] + phi[n]
        theta[n] = theta[n - 1]
        recip_p[n] = recip_p[n - 1]
        if is_prime[n]:
            theta[n] += log(n)
            recip_p[n] += 1.0 / n
    return M, A, G, Phi, theta, recip_p


def beats(edges: dict[tuple[int, int], int], i: int, j: int) -> bool:
    if i < j:
        return edges[(i, j)] == 1
    return edges[(j, i)] == 0


def score_hist(edges: dict[tuple[int, int], int], n: int) -> dict[int, int]:
    scores = Counter({i: 0 for i in range(n)})
    for i, j in combinations(range(n), 2):
        scores[i if beats(edges, i, j) else j] += 1
    return dict(sorted(Counter(scores.values()).items()))


def directed_3cycles(edges: dict[tuple[int, int], int], n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        out = Counter()
        for i, j in ((a, b), (a, c), (b, c)):
            out[i if beats(edges, i, j) else j] += 1
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def hamiltonian_count(edges: dict[tuple[int, int], int], n: int) -> int:
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for size in range(2, n + 1):
        nxt: dict[tuple[int, int], int] = {}
        for mask in range(1 << n):
            if mask.bit_count() != size:
                continue
            for last in range(n):
                if not ((mask >> last) & 1):
                    continue
                prev_mask = mask ^ (1 << last)
                total = 0
                for prev in range(n):
                    if ((prev_mask >> prev) & 1) and beats(edges, prev, last):
                        total += dp.get((prev_mask, prev), 0)
                if total:
                    nxt[(mask, last)] = total
        dp.update(nxt)
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def carrier_tournament():
    vertices = [
        (
            "raw_prime_count",
            (1, 0, 0, 0, 0, 0, 0, 0, 0),
            "counts primes but forgets residues, smoothing, and packet labels",
        ),
        (
            "mertens_mobius_cancellation",
            (1, 1, 0, 0, 0, 0, 0, 0, 0),
            "keeps signed squarefree cancellation but not unit capacity",
        ),
        (
            "mu_over_n_explicit_tail",
            (1, 1, 1, 0, 0, 0, 0, 0, 0),
            "normalizes Mobius cancellation by scale; useful tail guard",
        ),
        (
            "squarefree_unit_normalizer_G",
            (1, 1, 1, 1, 0, 0, 0, 0, 0),
            "G(x)=sum mu^2/phi keeps inverse primitive-unit packet cost",
        ),
        (
            "selberg_quadratic_upper_packet",
            (1, 1, 1, 1, 1, 0, 0, 0, 0),
            "quadratic upper-bound sieve declares lambda_d and kernel",
        ),
        (
            "large_sieve_minor_arc_packet",
            (1, 1, 1, 1, 1, 1, 0, 0, 0),
            "large sieve controls families of phases/exceptions",
        ),
        (
            "circle_method_major_minor_split",
            (1, 1, 1, 1, 1, 1, 1, 0, 0),
            "major/minor arcs retain local main term plus oscillatory error",
        ),
        (
            "explicit_formula_smoothing_packet",
            (1, 1, 1, 1, 1, 1, 1, 1, 0),
            "smoothing transform plus zero/exponential-sum terms are declared",
        ),
        (
            "kaczynski_boundary_approach_packet",
            (1, 1, 1, 1, 1, 1, 1, 1, 1),
            "exceptional/boundary approach class is retained as proof data",
        ),
        (
            "labelled_lrc_interval_packet",
            (2, 2, 2, 2, 2, 2, 2, 2, 2),
            "HYP-2981-style packet: endpoint, route, interval, and handoff labels",
        ),
    ]

    edges: dict[tuple[int, int], int] = {}
    for i, j in combinations(range(len(vertices)), 2):
        vi = vertices[i][1]
        vj = vertices[j][1]
        if vi > vj:
            edges[(i, j)] = 1
        elif vj > vi:
            edges[(i, j)] = 0
        else:
            edges[(i, j)] = 1
    return vertices, edges


def main() -> None:
    print("HYP-2982: analytic sieve packet weights for LRC14")
    print("=" * 72)
    print(f"sieve limit N={N}")
    primes, mu, phi = linear_mu_phi(N)
    M, A, G, Phi, theta, recip_p = prefix_tables(primes, mu, phi)

    print("\nArithmetic sums at LRC/Goldbach packet checkpoints")
    print(
        "Q        pi(Q)   theta/Q     M(Q)    M/sqrtQ      sum_mu/n"
        "      G=sum_mu2/phi   G-logQ     Phi(Q)"
    )
    for q in CHECKPOINTS:
        piq = sum(1 for p in primes if p <= q)
        print(
            f"{q:6d} {piq:8d} {theta[q]/q:10.6f} {M[q]:8d}"
            f" {M[q]/sqrt(q):11.6f} {A[q]:13.8f}"
            f" {G[q]:16.8f} {G[q]-log(q):10.6f} {Phi[q]:10d}"
        )

    print("\nLRC exact-period denominator packet readings")
    print("q     factorization        phi(q)  mu(q)  mu2/phi(q)  Phi(q)  G(q)  Selberg 1/G")
    for q in LRC_QS:
        fac = "*".join(
            f"{p}^{e}" if e > 1 else str(p)
            for p, e in prime_factorization(q, primes).items()
        )
        mu2_over_phi = (mu[q] * mu[q]) / phi[q] if phi[q] else 0.0
        print(
            f"{q:3d}  {fac:<18s} {phi[q]:7d} {mu[q]:6d}"
            f" {mu2_over_phi:12.8f} {Phi[q]:7d} {G[q]:7.4f} {1/G[q]:11.8f}"
        )
    print(
        "Squarefree blindness warning: mu^2/phi keeps q=14 and prime q=41, "
        "but zeroes prime-power or repeated-prime exact-period packets such as "
        "25,27,36,63,84,98,168,280,4312.  Those denominators need a "
        "Ramanujan/endpoint/Fejer/Kaczynski side channel before the quotient is theorem-safe."
    )

    print("\nLarge-sieve / Selberg normalizer readout")
    print(
        "G(z)=sum_{d<=z} mu(d)^2/phi(d) is the dimension-one inverse-unit"
        " normalizer.  Numerically G(z)-log z is already slowly stabilizing,"
        " so inverse primitive-unit cost is logarithmic after squarefree"
        " normalization, while Phi(z)=sum phi(q) is quadratic packet capacity."
    )
    for z in [14, 27, 41, 128, 1024, 8192, 65536, 200000]:
        print(
            f"z={z:6d}: G={G[z]:10.6f}, 1/G={1/G[z]:10.7f}, "
            f"log(z)={log(z):10.6f}, Phi/(z^2)={Phi[z]/(z*z):.6f}"
        )

    print("\nPrompt bridge: Goldbach explicit machinery as packet guardrails")
    bridge_rows = [
        (
            "sums over primes",
            "theta(x), prime races, local residues",
            "raw prime count is too coarse; residues and phase history matter",
        ),
        (
            "sum mu(n)",
            "squarefree signed cancellation",
            "Mertens cancellation is a tail guard, not a local certificate",
        ),
        (
            "sum mu(n)/n",
            "scale-normalized Mobius tail",
            "candidate error budget for quotient kernels",
        ),
        (
            "sum mu(n)^2/phi(n)",
            "inverse primitive-unit capacity",
            "Selberg/large-sieve normalizer for denominator packets",
        ),
        (
            "large sieve",
            "family bound over phases/residues",
            "preconditions late-q packets but must keep exceptional labels",
        ),
        (
            "circle method",
            "major arcs plus minor arcs",
            "same split as LRC primal twist witnesses vs blocker hypergraph",
        ),
        (
            "exponential sums",
            "oscillatory core",
            "LRC analogue is Fejer/Toeplitz negative quadratic form",
        ),
        (
            "explicit formula / saddle point",
            "smoothing transform plus zero terms",
            "smoothing is a proof object; changing it changes the packet",
        ),
        (
            "Kaczynski/Kaczorowski",
            "boundary approach / Goldbach exceptional sets",
            "exceptional-set boundary must be labelled, not averaged away",
        ),
    ]
    for name, carrier, lesson in bridge_rows:
        print(f"- {name}: carrier={carrier}; lesson={lesson}.")

    print("\nTournament Analysis")
    print("vertices are analytic proof carriers, not runners or arcs")
    print(
        "pairwise observable: retention vector = local data, Mobius sign,"
        " scale tail, unit capacity, quadratic kernel, family phase control,"
        " major/minor split, smoothing/zero packet, boundary approach class"
    )
    print("tie Hamiltonian path: listed vertex order")
    vertices, edges = carrier_tournament()
    n = len(vertices)
    path = " > ".join(name for name, _score, _why in reversed(vertices))
    print(f"score_hist={score_hist(edges, n)}")
    print(f"directed_3cycles={directed_3cycles(edges, n)}")
    print(f"hamiltonian_paths={hamiltonian_count(edges, n)}")
    print(f"Hamiltonian path: {path}")
    print("\nCarrier notes")
    for name, score, why in vertices:
        print(f"{name:<38s} score={score} :: {why}")

    print("\nLRC14 theorem-facing readout")
    print(
        "A Selberg/large-sieve packet is useful for LRC14 only as a middle"
        " certificate: it may bound a family of late denominators or"
        " exceptional phases, but HYP-2978/HYP-2981 still require endpoint,"
        " route, safe-measure, Fejer interval, Ramanujan, or state-lift labels."
    )
    print(
        "New subtarget: prove a labelled analytic packet lemma for the"
        " HYP-2901 Node-3 wall: after small denominator packets are killed,"
        " G(z)-type inverse-unit normalization plus a declared smoothing"
        " transform forces a twist/Fejer witness or routes to HYP-2908."
    )


if __name__ == "__main__":
    main()
