#!/usr/bin/env python3
"""
lrc14_support6_minkowski_count_codex_20260619.py

Codex 2026-06-19: first execution of the missing support-six
Minkowski/successive-minima count for the LRC(14) seven-sector route.

This does NOT claim to finish LRC(14).  It turns the open tail into data:

  * Work on THM-538's nonzero-offset lattice
        Lambda(E) = { n in Z^(k-1) : sum n_i e_i = 0 }.
    The zero speed is not a coordinate here; this avoids the free 0-speed
    Fourier coordinate noted in the S10 assembler.

  * Surviving THM-538 terms are exactly the coefficient vectors with at
    least six nonzero coordinates and no nonzero coordinate divisible by 7.
    Geometrically these are deleted anti-cosets: a fiber of the homomorphism
    n -> sum n_i e_i after deleting coordinate hyperplanes and 7-cosets.

  * For selected wide shapes, enumerate surviving relations in l_infty
    coefficient boxes, record shell counts, absolute Fourier-Kernel mass,
    a cheap c1-product envelope, and small successive-minima rank profiles.

  * Include Tournament Analysis, but with proof-obligations/routes as
    vertices rather than runners.  This explicitly challenges the local
    default that tournament vertices should be speeds or arcs.

Run:
    python3 04-computation/lrc14_support6_minkowski_count_codex_20260619.py
"""
from __future__ import annotations

import cmath
import itertools
import math
from collections import defaultdict
from fractions import Fraction as F
from functools import reduce
from math import comb, gcd

try:
    import sys
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass


TWO_PI_I = 2j * math.pi
C1 = 0.697303
PRIMES = (1000003, 1000033)
CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}
BK = {8: 16, 9: 13, 10: 12}


def banner(title: str) -> None:
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b) if a and b else (a or b)


def gcd_all(E: tuple[int, ...]) -> int:
    return reduce(gcd, [e for e in E if e], 0)


def measS7_fast(E) -> F:
    """Exact measure of hitting all seven fixed sectors."""
    E = tuple(sorted(set(int(e) for e in E)))
    nonzero = [e for e in E if e]
    if not nonzero:
        return F(0)
    D = 7 * reduce(lcm, nonzero, 1)
    cuts = {0, D}
    for e in nonzero:
        step = D // (7 * e)
        for x in range(0, D + 1, step):
            cuts.add(x)
    cuts = sorted(cuts)
    total = F(0)
    for a, b in zip(cuts, cuts[1:]):
        if b <= a:
            continue
        num = a + b
        den = 2 * D
        sectors = {0}
        for e in nonzero:
            sectors.add((7 * e * num) // den % 7)
        if len(sectors) == 7:
            total += F(b - a, D)
    return total


def M7(k: int) -> F:
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1)
               for t in range(7))


SUBSETS = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SIGNS = {T: (-1) ** len(T) for T in SUBSETS}
_CHAT: dict[tuple[int, tuple[int, ...]], complex] = {}


def shat(n: int, j: int) -> complex:
    if n == 0:
        return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) -
            cmath.exp(-TWO_PI_I * n * (a + 1.0 / 7.0))) / (TWO_PI_I * n)


def chat(n: int, T: tuple[int, ...]) -> complex:
    key = (n, T)
    if key in _CHAT:
        return _CHAT[key]
    if n == 0:
        val = complex(1 - len(T) / 7.0, 0.0)
    elif n % 7 == 0:
        val = 0j
    else:
        val = -sum(shat(n, j) for j in T)
    _CHAT[key] = val
    return val


def K_value(n: tuple[int, ...]) -> complex:
    total = 0j
    for T in SUBSETS:
        prod = 1 + 0j
        for ni in n:
            prod *= chat(ni, T)
            if prod == 0:
                break
        total += SIGNS[T] * prod
    return total


def support_info(v: tuple[int, ...]) -> tuple[bool, int, int]:
    """Return (survives_7_cosets, support_count, denominator_product)."""
    s = 0
    prod = 1
    for x in v:
        if x == 0:
            continue
        if x % 7 == 0:
            return False, 0, 1
        s += 1
        prod *= abs(x)
    return True, s, prod


def rank_mod(rows: list[tuple[int, ...]], p: int) -> int:
    if not rows:
        return 0
    mat = [[x % p for x in row] for row in rows]
    m = len(mat)
    n = len(mat[0])
    r = 0
    for c in range(n):
        pivot = None
        for i in range(r, m):
            if mat[i][c] % p:
                pivot = i
                break
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        inv = pow(mat[r][c], p - 2, p)
        mat[r] = [(x * inv) % p for x in mat[r]]
        for i in range(m):
            if i != r and mat[i][c] % p:
                factor = mat[i][c] % p
                mat[i] = [(mat[i][j] - factor * mat[r][j]) % p for j in range(n)]
        r += 1
        if r == m:
            break
    return r


def independent(rows: list[tuple[int, ...]], v: tuple[int, ...]) -> bool:
    trial = rows + [v]
    return all(rank_mod(trial, p) > len(rows) for p in PRIMES)


def add_to_basis(rows: list[tuple[int, ...]], v: tuple[int, ...]) -> None:
    # Every row lies in one linear relation kernel, so d-1 independent rows
    # already span the possible rational row space for this E.
    if len(rows) >= max(0, len(v) - 1):
        return
    if independent(rows, v):
        rows.append(v)


def enumerate_support6(E: tuple[int, ...], H: int, cutoff: int | None = None,
                       exact_height: int = 3) -> dict:
    """Meet-in-the-middle enumeration of support-six relations in |n|_inf <= H."""
    nonzero = tuple(e for e in E if e)
    d = len(nonzero)
    split = d // 2
    left_e = nonzero[:split]
    right_e = nonzero[split:]
    left = defaultdict(list)
    rng = range(-H, H + 1)

    for lt in itertools.product(rng, repeat=len(left_e)):
        s = sum(a * b for a, b in zip(lt, left_e))
        left[s].append(lt)

    shell_count = defaultdict(int)
    shell_bound = defaultdict(float)
    shell_exact_abs = defaultdict(float)
    shell_exact_signed = defaultdict(float)
    support_hist = defaultdict(int)
    type2_count = defaultdict(int)
    type2_bound = defaultdict(float)
    type2_exact_abs = defaultdict(float)
    type2_exact_signed = defaultdict(float)
    basis_by_height: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    total_matches = 0
    survivors = 0
    killed_7coset = 0
    min_height = None
    min_type2_height = None
    examples_by_height: dict[int, tuple[int, ...]] = {}
    type2_examples_by_height: dict[int, tuple[int, ...]] = {}

    large_idx = set()
    if cutoff is not None:
        large_idx = {i for i, e in enumerate(nonzero) if e > cutoff}

    for rt in itertools.product(rng, repeat=len(right_e)):
        s = sum(a * b for a, b in zip(rt, right_e))
        for lt in left.get(-s, ()):
            v = lt + rt
            if not any(v):
                continue
            total_matches += 1
            ok, supp, denom = support_info(v)
            if not ok:
                killed_7coset += 1
                continue
            if supp < 6:
                continue
            survivors += 1
            h = max(abs(x) for x in v)
            min_height = h if min_height is None else min(min_height, h)
            shell_count[h] += 1
            support_hist[(supp, h)] += 1
            # A deliberately simple THM-538 envelope:
            # |K(n)| <= sum_T prod_i |chat(n_i,T)| <= 64*c1^supp/prod |n_i|.
            bound = 64.0 * (C1 ** supp) / denom
            shell_bound[h] += bound
            examples_by_height.setdefault(h, v)
            add_to_basis(basis_by_height[h], v)

            is_type2 = bool(large_idx and any(v[i] for i in large_idx))
            if is_type2:
                type2_count[h] += 1
                type2_bound[h] += bound
                min_type2_height = h if min_type2_height is None else min(min_type2_height, h)
                type2_examples_by_height.setdefault(h, v)

            if h <= exact_height:
                kval = K_value(v)
                aval = abs(kval)
                shell_exact_abs[h] += aval
                shell_exact_signed[h] += kval.real
                if is_type2:
                    type2_exact_abs[h] += aval
                    type2_exact_signed[h] += kval.real

    rank_profile = []
    running_basis: list[tuple[int, ...]] = []
    for h in range(1, H + 1):
        for row in basis_by_height.get(h, ()):
            add_to_basis(running_basis, row)
        rank_profile.append((h, len(running_basis)))

    return {
        "E": E,
        "nonzero": nonzero,
        "H": H,
        "cutoff": cutoff,
        "d": d,
        "rank": max(0, d - 1),
        "covolume": math.sqrt(sum(e * e for e in nonzero)) / max(1, gcd_all(E)),
        "total_matches": total_matches,
        "survivors": survivors,
        "killed_7coset": killed_7coset,
        "min_height": min_height,
        "min_type2_height": min_type2_height,
        "shell_count": dict(shell_count),
        "shell_bound": dict(shell_bound),
        "shell_exact_abs": dict(shell_exact_abs),
        "shell_exact_signed": dict(shell_exact_signed),
        "support_hist": dict(support_hist),
        "type2_count": dict(type2_count),
        "type2_bound": dict(type2_bound),
        "type2_exact_abs": dict(type2_exact_abs),
        "type2_exact_signed": dict(type2_exact_signed),
        "rank_profile": rank_profile,
        "examples": sorted((h, v) for h, v in examples_by_height.items()),
        "type2_examples": sorted((h, v) for h, v in type2_examples_by_height.items()),
    }


def fmt_shells(stats: dict, key: str, digits: int = 6) -> str:
    data = stats[key]
    if not data:
        return "{}"
    parts = []
    for h in sorted(data):
        val = data[h]
        if isinstance(val, int):
            parts.append(f"{h}:{val}")
        else:
            parts.append(f"{h}:{val:.{digits}g}")
    return "{" + ", ".join(parts) + "}"


def report_case(name: str, E: tuple[int, ...], H: int, exact_height: int = 3) -> dict:
    k = len(E)
    cutoff = BK.get(k)
    stats = enumerate_support6(E, H=H, cutoff=cutoff, exact_height=exact_height)
    m = measS7_fast(E)
    consec = measS7_fast(tuple(range(k)))
    delta = m - M7(k)
    print(f"\nCASE {name}")
    print(f"  E={E}")
    print(f"  k={k}, span={max(E)}, gcd={gcd_all(E)}, B(k)={cutoff}")
    print(f"  meas(S7)={m}={float(m):.6f}; consec={float(consec):.6f}; "
          f"cap={float(CAP.get(k, F(1))):.6f}; corr=meas-M7={float(delta):.6f}")
    print(f"  nonzero dimension d={stats['d']}, relation rank={stats['rank']}, "
          f"covolume ||e||={stats['covolume']:.3f}")
    print(f"  box H={H}: raw relations={stats['total_matches']}, "
          f"surviving support>=6={stats['survivors']}, killed nonzero 7-coset={stats['killed_7coset']}")
    print(f"  first surviving height={stats['min_height']}; "
          f"first large-involving height={stats['min_type2_height']}")
    print(f"  shell counts       {fmt_shells(stats, 'shell_count')}")
    print(f"  shell abs envelope {fmt_shells(stats, 'shell_bound')}")
    print(f"  exact |K| shells h<={exact_height} {fmt_shells(stats, 'shell_exact_abs')}")
    print(f"  exact K shells  h<={exact_height} {fmt_shells(stats, 'shell_exact_signed')}")
    print(f"  type-II counts       {fmt_shells(stats, 'type2_count')}")
    print(f"  type-II abs envelope {fmt_shells(stats, 'type2_bound')}")
    print(f"  type-II exact |K| h<={exact_height} {fmt_shells(stats, 'type2_exact_abs')}")
    print(f"  type-II exact K  h<={exact_height} {fmt_shells(stats, 'type2_exact_signed')}")
    print(f"  rank profile from survivor shells: {stats['rank_profile']}")
    if stats["type2_examples"]:
        print("  first type-II anti-coset witnesses:")
        for h, v in stats["type2_examples"][:3]:
            s = sum(a * b for a, b in zip(v, stats["nonzero"]))
            print(f"    h={h}, n={v}, dot={s}")
    return stats


def tournament_analysis(rows: list[dict]) -> None:
    banner("TOURNAMENT ANALYSIS: proof-obligation vertices, not runner vertices")
    print("  Pairwise observable: route A -> B if A satisfies more of")
    print("  {uniform, computable, signed, multilarge, finite-check-compatible}.")
    print("  Switch/gauge: flip the edge if the session prioritizes experimental")
    print("  discovery over proof closure.  Tie Hamiltonian path: score order, then name.")

    traits = {
        "bounded_finite_check":      (1, 1, 0, 1, 1),
        "anti_coset_shell_count":    (0, 1, 1, 1, 1),
        "successive_minima_theta":   (1, 0, 1, 1, 1),
        "single_stranger_weyl":      (1, 1, 0, 0, 1),
        "cluster_collapse":         (0, 1, 0, 1, 1),
        "free_coordinate_envelope":  (0, 1, 0, 1, 0),
    }
    names = sorted(traits)
    score = {name: sum(traits[name]) for name in names}
    edges = {}
    outdeg = {name: 0 for name in names}
    for a, b in itertools.combinations(names, 2):
        if (score[a], a) >= (score[b], b):
            edges[(a, b)] = a
            outdeg[a] += 1
        else:
            edges[(a, b)] = b
            outdeg[b] += 1

    hist = defaultdict(int)
    for v in outdeg.values():
        hist[v] += 1
    ham = sorted(names, key=lambda x: (-score[x], x))
    cycles = []
    for a, b, c in itertools.combinations(names, 3):
        def beats(x, y):
            return edges[tuple(sorted((x, y)))] == x
        if beats(a, b) and beats(b, c) and beats(c, a):
            cycles.append((a, b, c))
        if beats(a, c) and beats(c, b) and beats(b, a):
            cycles.append((a, c, b))
    print(f"  score histogram (outdegree -> count): {dict(sorted(hist.items()))}")
    print(f"  directed 3-cycles: {len(cycles)}")
    print(f"  tie Hamiltonian path: {' -> '.join(ham)}")
    print("  SCCs: singleton/transitive under this proof-closure gauge.")
    print()
    print("  Assumption challenge:")
    print("    Tried vertices = speeds, relations, deleted anti-cosets, proof obligations.")
    print("    Kept quotient = anti-coset shell count because it preserves the signed")
    print("    THM-538 correction predicate K(n)!=0.  It destroys x-space component")
    print("    geometry, which is exactly why single-stranger Weyl remains a parallel")
    print("    route rather than a consequence of the count.")
    print()
    print("  Family difficulty ranking by first type-II survivor height:")
    fam = []
    for st in rows:
        fam.append((st["min_type2_height"] if st["min_type2_height"] is not None else 999,
                    st["covolume"], st["E"]))
    for h, covol, E in sorted(fam):
        print(f"    h2={h}, covol={covol:.2f}, E={E}")


def main() -> None:
    banner("LRC(14) support-six Minkowski count: deleted anti-cosets")
    print("  THM-538 means the Fourier tail is not a free product sum.")
    print("  It is a theta sum over Lambda(E), after deleting coordinate hyperplanes")
    print("  and all nonzero coefficient 7-cosets.  The count below is the first")
    print("  finite execution of that missing geometry-of-numbers object.")

    banner("Exact target constants")
    for k in (8, 9, 10):
        consec = measS7_fast(tuple(range(k)))
        print(f"  k={k}: M7={float(M7(k)):.6f}, consec={consec}={float(consec):.6f}, "
              f"cap={CAP[k]}={float(CAP[k]):.6f}, cap-consec={float(CAP[k]-consec):.6f}")

    cases = [
        ("k8 worst wide verifier: 21 = 1+2+3+4+5+6",
         (0, 1, 2, 3, 4, 5, 6, 21), 8, 4),
        ("k8 dissociated one-stranger",
         (0, 1, 2, 3, 4, 5, 6, 97), 6, 3),
        ("k8 no-scale shifted cluster",
         (0, 50, 51, 52, 53, 54, 55, 56), 4, 3),
        ("k9 worst verifier one-stranger",
         (0, 1, 2, 3, 4, 5, 6, 7, 68), 5, 3),
        ("k10 worst verifier: 22 = 1+2+3+4+5+7",
         (0, 1, 2, 3, 4, 5, 6, 7, 8, 22), 4, 3),
    ]
    rows = [report_case(name, E, H, exact_h) for name, E, H, exact_h in cases]

    banner("Interpretation")
    print("  1. The support-six floor is real: the first surviving anti-cosets are")
    print("     sparse and visible, not a diffuse harmonic product.")
    print("  2. Span alone is not a lower bound for the first survivor.  The k=8")
    print("     verifier worst case has a height-1 large-involving relation because")
    print("     21 = 1+2+3+4+5+6.  This is the anti-coset obstruction to the naive")
    print("     'large speed forces large coefficient' proof.")
    print("  3. Dissociated strangers move the first type-II shell outward quickly.")
    print("     That supports HYP-2610's single-stranger Weyl/contraction route.")
    print("  4. No-scale shifted clusters have small relations everywhere, but their")
    print("     exact sector measure is tiny.  They should be handled by collapse or")
    print("     by quotienting common translation, not by span-only Minkowski decay.")
    print("  5. The next rigorous lemma should count deleted anti-cosets by additive")
    print("     energy level: finite-check the low-height identities like")
    print("     N=sum(six core offsets), then apply a true successive-minima theta")
    print("     bound to the remaining dissociated classes.")

    tournament_analysis(rows)

    banner("Session verdict")
    print("  LRC(14) is still not proved in this run.")
    print("  The missing Minkowski count has been made concrete enough to split:")
    print("    (A) low-height anti-coset resonances, finite/additive-energy ledger;")
    print("    (B) dissociated deleted-coset theta tail, genuine successive minima;")
    print("    (C) no-scale clusters, collapse quotient rather than span decay.")


if __name__ == "__main__":
    main()
