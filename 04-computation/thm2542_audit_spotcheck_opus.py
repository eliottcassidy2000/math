#!/usr/bin/env python3
"""Independent hostile-audit spot-check for THM-2542 (opus, 2026-07-27).

Written from scratch; deliberately shares no code with the companion
referee lrc14_seven_chart_cech_holonomy_thm2542.py.

Target: the load-bearing overlap-cocycle composition of THM-2542, attacked
from the rawest available definitions:

  A. Rederive THM-2542 eq. (3) [= THM-2535 (34e)] directly from THM-2535's
     raw cut operator (11): R_{tau,a,b}(v) = sum_r d(v - tau*rep(a r + b), r),
     with symbolic positive rational masses.
  B. Build the transported chart coordinate maps with their actual domains
     (teeth 0 and 1, i.e. clock rows {k, k+eta}), verify every overlap map
     is a genuinely CONSTANT translation by s-t on the whole 13-point deck,
     that non-adjacent charts do not overlap, and that every triple overlap
     is empty (so the Cech cocycle condition holds vacuously - checked, not
     assumed).  Also verify the sharpness hostile: with UNRESTRICTED chart
     domains (all seven teeth) the coordinate change is NOT constant, so the
     stated domain restriction is load-bearing for the seven-cycle nerve.
  C. Compose the seven overlap maps around the oriented cycle as explicit
     permutation composition; assert the holonomy is translation by 7(s-t)
     and fixed-point-free, for all 936 triples (eta, s, t).
  D. EXHAUSTIVE coboundary search: enumerate all 13^6 gauge vectors with
     h_0 = 0 (these hit every coboundary of the 7-cycle) and assert the only
     constant coboundary is zero.  Hence the constant cochain (a,...,a),
     a != 0, is not a coboundary: the class is nonzero by exhaustion, with
     no use of the sum shortcut.
  E. Affine-gauge / carry-covariance hostile: gauge single charts by every
     element of AGL_1(F_13) (156 per chart) and full random gauge vectors;
     assert the holonomy permutation never becomes the identity and matches
     the conjugation prediction x -> x + A_0*7(s-t).
  F. Cover hostiles by exact linear chase (no sum shortcut): every one-edge
     refinement (8-cycle) and a random full refinement (91-cycle) remain
     non-coboundaries; the n-fold cover pullback is non-trivializable for
     n = 1..12 and trivializable at n = 13 with an explicitly verified gauge.
  G. Characteristic-zero certificate: nonvanishing of the (11e) mode factors
     in Z[x]/Phi_7 and Z[y]/Phi_13 (exact polynomial arithmetic), so the
     F_547 referee's nonvanishing claims also hold over the cyclotomics.

All arithmetic is exact (int / Fraction / polynomial over Z).  Exit 0 with
final PASS line iff every assertion holds.
"""

from fractions import Fraction
from itertools import combinations, product
import random

P = 13   # root deck
Q = 7    # clock cycle

random.seed(25420727)  # deterministic run


# ----------------------------------------------------------------------
# Part A: THM-2542 (3) from THM-2535's raw cut operator (11).
# ----------------------------------------------------------------------
def rep(x):
    """Canonical representative of x mod 7 in {0,...,6}, lifted to Z."""
    return x % Q


def cut_operator(d, tau, a, b):
    """THM-2535 (11): R_{tau,a,b}(v) = sum_r d(v - tau*rep(a r + b), r)."""
    return {
        v: sum(d.get(((v - tau * rep(a * r + b)) % P, r), Fraction(0))
               for r in range(Q))
        for v in range(P)
    }


def part_a():
    checks = 0
    for am in range(1, P):                      # marker root a != 0
        xs = [Fraction(random.randint(1, 97), random.randint(1, 89))
              for _ in range(Q)]                # positive rational masses
        for kappa in range(Q):
            # THM-2535 (34c): two-clock signed table, rows {kappa, kappa+1}.
            d = {(am, kappa): xs[kappa],
                 (am, (kappa + 1) % Q): -xs[(kappa + 1) % Q]}
            # THM-2535 (34d): tau_edge=-a, cut scale 1, cut intercept -kappa.
            r_chart = cut_operator(d, (-am) % P, 1, (-kappa) % Q)
            for v in range(P):
                expected = Fraction(0)
                if v == am:
                    expected += xs[kappa]
                if v == 0:
                    expected -= xs[(kappa + 1) % Q]
                assert r_chart[v] == expected, (am, kappa, v)
                checks += 1
    return checks  # 12 * 7 * 13 = 1092


# ----------------------------------------------------------------------
# Part B: charts with true domains; overlap maps; nerve structure.
# ----------------------------------------------------------------------
def chart_coord(k, h, r, eta, s, t):
    """Chart k reads deck point h at clock row r (row in chart k's domain).

    THM-2542 (11b) with the type repair made explicit: the tooth index
    eta^{-1}(r-k) is computed in F_7 and lifted by rep() to {0,...,6}
    before multiplying the F_13 element (t-s).
    """
    eta_inv = pow(eta, -1, Q)
    tooth = rep(eta_inv * (r - k))          # integer in 0..6
    return (h + (t - s) * tooth) % P


def part_b():
    overlap_checks = 0
    nerve_checks = 0
    seam_checks = 0
    for eta in range(1, Q):
        domains = {k: {k % Q, (k + eta) % Q} for k in range(Q)}
        # nerve: chart k-eta and chart k overlap exactly at row k; all other
        # pairs are disjoint; every triple overlap is empty.
        for j, k in combinations(range(Q), 2):
            inter = domains[j] & domains[k]
            if (j - k) % Q in (eta, Q - eta):
                assert len(inter) == 1
            else:
                assert inter == set()
            nerve_checks += 1
        for j, k, l in combinations(range(Q), 3):
            assert domains[j] & domains[k] & domains[l] == set()
            nerve_checks += 1
        for s in range(P):
            for t in range(P):
                if s == t:
                    continue
                for k in range(Q):
                    km = (k - eta) % Q
                    # overlap row is k; transition chart km -> chart k.
                    pairs = [(chart_coord(km, h, k, eta, s, t),
                              chart_coord(k, h, k, eta, s, t))
                             for h in range(P)]
                    deltas = {(v2 - v1) % P for v1, v2 in pairs}
                    assert deltas == {(s - t) % P}, (eta, s, t, k)
                    assert len({v1 for v1, _ in pairs}) == P  # bijection
                    overlap_checks += 1
                # event bookkeeping (THM-2542 sect. 1): the clock-k event
                # sits at root t in the incoming chart k-eta and at root s
                # in the outgoing chart k.
                for k in range(Q):
                    assert chart_coord(k, s, k, eta, s, t) == s
                    assert chart_coord((k - eta) % Q, s, k, eta, s, t) == t
        # Sharpness hostile: WITHOUT the domain restriction the coordinate
        # change chart(k-eta) -> chart(k) on all seven rows is NOT constant:
        # translation by s-t on six rows and by 7(s-t) on the wrap row.
        s, t = 1, 0
        for k in range(Q):
            km = (k - eta) % Q
            wrap = 0
            plain = 0
            for r in range(Q):
                d1 = (chart_coord(k, 0, r, eta, s, t)
                      - chart_coord(km, 0, r, eta, s, t)) % P
                if d1 == (s - t) % P:
                    plain += 1
                elif d1 == (7 * (s - t)) % P:
                    wrap += 1
            assert plain == 6 and wrap == 1, (eta, k)
            assert (7 * (s - t)) % P != (s - t) % P
            seam_checks += 1
    return overlap_checks, nerve_checks, seam_checks


# ----------------------------------------------------------------------
# Part C: holonomy by explicit permutation composition.
# ----------------------------------------------------------------------
def transition_perm(k, eta, s, t, gauges=None):
    """Overlap permutation chart(k-eta)-coords -> chart(k)-coords at row k.

    gauges: optional dict chart -> (A, H) in AGL_1(F_13) acting on that
    chart's coordinates by v -> A v + H.
    """
    km = (k - eta) % Q
    perm = {}
    for h in range(P):
        v1 = chart_coord(km, h, k, eta, s, t)
        v2 = chart_coord(k, h, k, eta, s, t)
        if gauges:
            a1, h1 = gauges[km]
            a2, h2 = gauges[k]
            v1 = (a1 * v1 + h1) % P
            v2 = (a2 * v2 + h2) % P
        perm[v1] = v2
    assert sorted(perm) == list(range(P))
    return perm


def holonomy_perm(eta, s, t, gauges=None):
    """Compose transitions along chart path 0 -> eta -> 2eta -> ... -> 0."""
    hol = {v: v for v in range(P)}
    for step in range(1, Q + 1):
        k = (step * eta) % Q            # arriving chart; final step k = 0
        tr = transition_perm(k, eta, s, t, gauges)
        hol = {v: tr[hol[v]] for v in range(P)}
    return hol


def part_c():
    checks = 0
    for eta in range(1, Q):
        for s in range(P):
            for t in range(P):
                if s == t:
                    continue
                hol = holonomy_perm(eta, s, t)
                expected = {v: (v + 7 * (s - t)) % P for v in range(P)}
                assert hol == expected, (eta, s, t)
                assert all(hol[v] != v for v in range(P))
                checks += 1
    return checks  # 936


# ----------------------------------------------------------------------
# Part D: exhaustive coboundary enumeration on the 7-cycle.
# ----------------------------------------------------------------------
def part_d():
    # Coverage lemma: delta(h + c*1) = delta(h); verified on 1000 samples.
    for _ in range(1000):
        h = [random.randrange(P) for _ in range(Q)]
        c = random.randrange(P)
        d1 = [(h[k] - h[k - 1]) % P for k in range(Q)]
        hc = [(x + c) % P for x in h]
        d2 = [(hc[k] - hc[k - 1]) % P for k in range(Q)]
        assert d1 == d2
    # Exhaust all coboundaries via h_0 = 0.
    constant_coboundaries = 0
    total = 0
    for tail in product(range(P), repeat=Q - 1):
        total += 1
        h = (0,) + tail
        c = tail[0]                     # h_1 - h_0
        ok = True
        for k in range(2, Q):
            if (h[k] - h[k - 1]) % P != c:
                ok = False
                break
        if ok and (h[0] - h[Q - 1]) % P != c:
            ok = False
        if ok:
            constant_coboundaries += 1
            assert c == 0 and h == (0,) * Q
    assert total == P ** (Q - 1)
    assert constant_coboundaries == 1
    return total


# ----------------------------------------------------------------------
# Part E: affine-gauge (carry-covariance) hostile.
# ----------------------------------------------------------------------
def part_e():
    identity = {v: v for v in range(P)}
    single_checks = 0
    # Exhaustive single-chart AGL_1(F_13) gauges on the 12 THM-2535 packets.
    for s in range(1, P):               # (eta, s, t) = (1, a, 0)
        for k in range(Q):
            for A in range(1, P):
                for H in range(P):
                    gauges = {j: (1, 0) for j in range(Q)}
                    gauges[k] = (A, H)
                    hol = holonomy_perm(1, s, 0, gauges)
                    assert hol != identity, (s, k, A, H)
                    single_checks += 1
    # Random full gauge vectors on random general triples, with the
    # conjugation prediction Hol' = f_0 o Hol o f_0^{-1}.
    full_checks = 0
    for _ in range(500):
        eta = random.randrange(1, Q)
        s = random.randrange(P)
        t = random.randrange(P)
        if s == t:
            continue
        gauges = {j: (random.randrange(1, P), random.randrange(P))
                  for j in range(Q)}
        hol = holonomy_perm(eta, s, t, gauges)
        a0, h0 = gauges[0]
        predicted = {(a0 * v + h0) % P:
                     (a0 * ((v + 7 * (s - t)) % P) + h0) % P
                     for v in range(P)}
        assert hol == predicted
        assert hol != identity
        full_checks += 1
    return single_checks, full_checks


# ----------------------------------------------------------------------
# Part F: refinement and cover hostiles by exact linear chase.
# ----------------------------------------------------------------------
def cycle_coboundary_solvable(g):
    """Decide (exactly) whether the cochain g on an n-cycle is a coboundary.

    Edges are e_i: vertex i -> vertex i+1 with value g[i] = h_{i+1} - h_i.
    Propagate from h_0 = 0; solvable iff the chase closes.  The general
    solution differs by a constant, so this decision is complete.
    """
    n = len(g)
    h = [0] * n
    for i in range(n - 1):
        h[i + 1] = (h[i] + g[i]) % P
    closes = (h[0] - h[n - 1]) % P == g[n - 1] % P
    return closes, h


def part_f():
    refine_checks = 0
    for a in range(1, P):
        base = [a] * Q
        for e in range(Q):
            for u in range(P):
                refined = base[:e] + [u, (a - u) % P] + base[e + 1:]
                assert len(refined) == Q + 1
                solvable, _ = cycle_coboundary_solvable(refined)
                assert not solvable, (a, e, u)
                refine_checks += 1
        # random full refinement into a 91-cycle
        refined = []
        for _ in range(Q):
            parts = [random.randrange(P) for _ in range(P - 1)]
            parts.append((a - sum(parts)) % P)
            refined.extend(parts)
        assert len(refined) == Q * P
        solvable, _ = cycle_coboundary_solvable(refined)
        assert not solvable, a
        refine_checks += 1
    cover_checks = 0
    for a in range(1, P):
        for n in range(1, P + 1):
            g = [a] * (Q * n)           # pullback to the n-fold cover
            solvable, h = cycle_coboundary_solvable(g)
            if n < P:
                assert not solvable, (a, n)
            else:
                assert solvable, a
                # verify the explicit trivializing gauge entrywise
                m = Q * P
                for i in range(m):
                    assert (h[(i + 1) % m] - h[i]) % P == a
            cover_checks += 1
    return refine_checks, cover_checks


# ----------------------------------------------------------------------
# Part G: characteristic-zero nonvanishing certificates.
# ----------------------------------------------------------------------
def poly_mod(coeffs, mod_coeffs):
    """Reduce polynomial (coeff list, low->high) mod monic mod_coeffs."""
    coeffs = list(coeffs)
    n = len(mod_coeffs) - 1
    while len(coeffs) > n:
        lead = coeffs.pop()
        if lead:
            for i in range(n):
                coeffs[-n + i] -= lead * mod_coeffs[i]
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def part_g():
    phi7 = [1, 1, 1, 1, 1, 1, 1]    # Phi_7(x) = 1 + x + ... + x^6
    checks = 0
    # 1 - x^m != 0 in Z[x]/Phi_7 for m = 1..6  (x = primitive 7th root)
    for m in range(1, 7):
        poly = [1] + [0] * (m - 1) + [-1]   # 1 - x^m, low -> high
        red = poly_mod(poly, phi7)
        assert red, m
        checks += 1
    # y^i - y^j != 0 in Z[y]/Phi_13 for 0 <= i < j <= 12
    phi13 = [1] * 13
    for i in range(13):
        for j in range(i + 1, 13):
            poly = [0] * 13
            poly[i] += 1
            poly[j] -= 1
            red = poly_mod(poly, phi13)
            assert red, (i, j)
            checks += 1
    return checks


def main():
    a_checks = part_a()
    print(f"A_cut_operator_rederivation_checks={a_checks}")
    ov, nv, sm = part_b()
    print(f"B_constant_overlap_translation_checks={ov}")
    print(f"B_nerve_pair_and_triple_overlap_checks={nv}")
    print(f"B_unrestricted_domain_seam_hostile_checks={sm}")
    c_checks = part_c()
    print(f"C_holonomy_permutation_composition_checks={c_checks}")
    d_total = part_d()
    print(f"D_exhaustive_gauges_enumerated={d_total} "
          "constant_coboundaries_found=1 (zero cochain only)")
    e_single, e_full = part_e()
    print(f"E_single_chart_affine_gauge_hostiles={e_single} "
          f"random_full_affine_gauges={e_full} holonomy_never_identity=1")
    f_refine, f_cover = part_f()
    print(f"F_refinement_noncoboundary_checks={f_refine} "
          f"cover_degree_checks={f_cover} minimal_trivializing_degree=13")
    g_checks = part_g()
    print(f"G_char_zero_cyclotomic_nonvanishing_checks={g_checks}")
    print("SPOTCHECK_ALL_EXACT=PASS")


if __name__ == "__main__":
    main()
