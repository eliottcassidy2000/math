#!/usr/bin/env python3
"""
LRC(14) support-6 relative count scout.

This is a targeted computation for the post-THM-538 gap:

    meas(S7(E)) = M7(k) + sum_{0 != n in Lambda^o(E)} K(n),
    K(n) = 0 unless at least six coordinates are nonzero and not multiples of 7.

The failed bound summed the per-coordinate envelope freely and saw a harmonic
divergence.  This script measures the coupled object instead: for each fixed
six-coordinate support J it sums over the hyperplane

    sum_{j in J} n_j e_j = 0,

with all n_j nonzero and 7-coprime.  The one linear equation should make height
shells decay like H^{-2}; that is the "Minkowski count" this scout tries to make
concrete.

The output is evidence and route-finding, not a proof.  In particular, exact
support-6 is only the leading layer; support 7+ low-height resonances are also
reported separately because they explain the familiar wide but safe spikes.
"""
from __future__ import annotations

import cmath
import itertools
import math
from collections import defaultdict
from fractions import Fraction as F
from math import comb, gcd, prod
from functools import lru_cache, reduce

TWO_PI_I = 2j * math.pi
C1 = 0.697303  # measured sharp value from THM-538 / S10 audit
ENV64 = 64.0


def gcd_all(xs):
    return reduce(gcd, [abs(x) for x in xs if x], 0)


def primitive(E):
    nz = [e for e in E if e]
    g = gcd_all(nz)
    return g == 1


def measS7(E):
    """Exact rational measure that the speeds in E hit all seven sectors."""
    E = sorted(set(E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            bps.add(F(m, 7 * e))
    bps = sorted(bps)
    total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7:
            total += x1 - x0
    return total


def M7(k):
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))


def shat(n, j):
    if n == 0:
        return 1.0 / 7.0
    a = j / 7.0
    return (
        cmath.exp(-TWO_PI_I * n * a)
        - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))
    ) / (TWO_PI_I * n)


SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1) ** len(T) for T in SUB}
SUB_INDEX = tuple(range(len(SUB)))
SGN_LIST = tuple(SGN[T] for T in SUB)


@lru_cache(None)
def chat(n, T):
    if n == 0:
        return complex(1 - len(T) / 7.0, 0.0)
    if n % 7 == 0:
        return 0j
    return -sum(shat(n, j) for j in T)


def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0:
                break
        s += SGN[T] * p
    return s


@lru_cache(None)
def coeffs(H):
    return tuple(x for x in range(-H, H + 1) if x != 0 and x % 7 != 0)


def triple_aggregate(vals, H):
    """
    Aggregate triples by their weighted sum.

    Return sum -> (count, sum 1/prod |n_i|).  This lets a six-support
    hyperplane sum be computed by one meet-in-the-middle convolution.
    """
    vals = tuple(vals)
    out = defaultdict(lambda: [0, 0.0])
    cs = coeffs(H)
    for ns in itertools.product(cs, repeat=len(vals)):
        s = sum(a * n for a, n in zip(vals, ns))
        inv = 1.0 / prod(abs(n) for n in ns)
        rec = out[s]
        rec[0] += 1
        rec[1] += inv
    return out


def support6_subset_sum(vals, H):
    """
    Exact-support-6 envelope sum for one six-speed subset.

    All six coefficients are nonzero and non-7.  The returned env is the THM-538
    absolute envelope 64*c1^6*sum prod 1/|n_i|.
    """
    vals = tuple(vals)
    left = triple_aggregate(vals[:3], H)
    right = triple_aggregate(vals[3:], H)
    count = 0
    invsum = 0.0
    for s, (c1, r1) in left.items():
        c2, r2 = right.get(-s, (0, 0.0))
        if c2:
            count += c1 * c2
            invsum += r1 * r2
    return count, ENV64 * (C1 ** 6) * invsum


def signed_tuple_aggregate(vals, H):
    """
    Aggregate tuple products by weighted sum for the exact signed K contribution.

    For each sum s this stores a 64-vector:

        A_T(s) = sum_{n: vals.n=s} prod_i chat(n_i,T).

    Pairing two halves on opposite sums then evaluates the signed support layer
    exactly, including the T-dependent cancellations.
    """
    vals = tuple(vals)
    out = {}
    cs = coeffs(H)
    for ns in itertools.product(cs, repeat=len(vals)):
        s = sum(a * n for a, n in zip(vals, ns))
        if s not in out:
            out[s] = [0j] * len(SUB)
        vec = out[s]
        for ti, T in enumerate(SUB):
            p = 1.0 + 0j
            for ni in ns:
                p *= chat(ni, T)
                if p == 0:
                    break
            vec[ti] += p
    return out


def signed_support6_subset_sum(vals, H, total_dim):
    """Exact signed K sum over one exact six-coordinate support."""
    vals = tuple(vals)
    left = signed_tuple_aggregate(vals[:3], H)
    right = signed_tuple_aggregate(vals[3:], H)
    zero_power = total_dim - 6
    zero_weight = tuple((chat(0, T) ** zero_power) for T in SUB)
    signed = 0j
    for s, lv in left.items():
        rv = right.get(-s)
        if rv is None:
            continue
        for ti in SUB_INDEX:
            signed += SGN_LIST[ti] * zero_weight[ti] * lv[ti] * rv[ti]
    return signed


def support6_profile(E, H, large_cutoff=None):
    """
    Sum exact-support-6 envelopes over all six-coordinate supports.

    If large_cutoff is supplied, also separately sum the supports that touch an
    offset larger than the cutoff.  The zero offset is omitted because it cannot
    appear in a nontrivial relation coordinate in the reduced offset lattice.
    """
    offsets = tuple(e for e in sorted(E) if e != 0)
    rows = []
    total_count = 0
    total_env = 0.0
    type_count = 0
    type_env = 0.0
    for idxs in itertools.combinations(range(len(offsets)), 6):
        vals = tuple(offsets[i] for i in idxs)
        count, env = support6_subset_sum(vals, H)
        total_count += count
        total_env += env
        is_type = large_cutoff is not None and any(v > large_cutoff for v in vals)
        if is_type:
            type_count += count
            type_env += env
        rows.append((vals, count, env, is_type))
    rows.sort(key=lambda row: row[2], reverse=True)
    return {
        "offsets": offsets,
        "rows": rows,
        "total_count": total_count,
        "total_env": total_env,
        "type_count": type_count,
        "type_env": type_env,
    }


def signed_support6_profile(E, H, large_cutoff=None):
    """Exact signed support-6 contribution over all six-coordinate supports."""
    offsets = tuple(e for e in sorted(E) if e != 0)
    total = 0j
    type_total = 0j
    rows = []
    for idxs in itertools.combinations(range(len(offsets)), 6):
        vals = tuple(offsets[i] for i in idxs)
        signed = signed_support6_subset_sum(vals, H, len(offsets))
        is_type = large_cutoff is not None and any(v > large_cutoff for v in vals)
        total += signed
        if is_type:
            type_total += signed
        rows.append((vals, signed, is_type))
    rows.sort(key=lambda row: abs(row[1]), reverse=True)
    return {"total": total, "type_total": type_total, "rows": rows}


def low_height_relations(E, Hmax=4, max_support=None):
    """
    Find low-height genuine relations of support >=6.

    We enumerate supports and solve for the final coefficient.  This catches
    the height-1 and height-2 resonances that exact support-6 alone can miss.
    """
    offsets = tuple(e for e in sorted(E) if e != 0)
    d = len(offsets)
    max_support = min(max_support or d, d)
    out = []
    for r in range(6, max_support + 1):
        for idxs in itertools.combinations(range(d), r):
            vals = tuple(offsets[i] for i in idxs)
            last = vals[-1]
            for H in range(1, Hmax + 1):
                cs = coeffs(H)
                cs_set = set(cs)
                count = 0
                examples = []
                absK = 0.0
                signedK = 0j
                for ns0 in itertools.product(cs, repeat=r - 1):
                    s = sum(a * n for a, n in zip(vals[:-1], ns0))
                    if s % last:
                        continue
                    nlast = -s // last
                    if nlast not in cs_set:
                        continue
                    if max(abs(x) for x in (*ns0, nlast)) != H:
                        continue
                    full = [0] * d
                    for local_i, global_i in enumerate(idxs[:-1]):
                        full[global_i] = ns0[local_i]
                    full[idxs[-1]] = nlast
                    kv = Kk(tuple(full))
                    count += 1
                    absK += abs(kv)
                    signedK += kv
                    if len(examples) < 3:
                        examples.append(tuple(full))
                if count:
                    out.append(
                        {
                            "support": vals,
                            "support_size": r,
                            "height": H,
                            "count_on_shell": count,
                            "absK_shell": absK,
                            "signedK_shell": signedK,
                            "examples": examples,
                        }
                    )
                    break
    out.sort(key=lambda row: (row["height"], row["support_size"], -row["absK_shell"]))
    return out


def participation_tournament(profile, topn=8):
    """
    Tournament Analysis diagnostic.

    Vertices are nonzero offsets.  The pairwise observable is relative
    exact-support-6 envelope participation:

        p_i - p_j,

    where p_i sums the support-6 envelope over supports containing offset i.
    The binary relation points toward the larger participant.  Ties use the
    Hamiltonian path induced by increasing offset.  Since this gauge collapses
    to a one-dimensional score, cycles are not expected; that is itself a
    useful sanity check for this diagnostic.
    """
    offsets = profile["offsets"]
    p = {v: 0.0 for v in offsets}
    for vals, _count, env, _is_type in profile["rows"]:
        for v in vals:
            p[v] += env
    ordered = sorted(offsets, key=lambda v: (-p[v], v))
    score = {v: 0 for v in offsets}
    flips = 0
    for a, b in itertools.combinations(offsets, 2):
        if p[a] > p[b] or (p[a] == p[b] and a < b):
            score[a] += 1
        else:
            score[b] += 1
            if a < b:
                flips += 1
    hist = defaultdict(int)
    for v in offsets:
        hist[score[v]] += 1
    return {
        "ordered": ordered[:topn],
        "participation": p,
        "score_hist": dict(sorted(hist.items())),
        "flips_against_offset_order": flips,
    }


def fmt(x):
    return f"{x:.6g}"


def section(title):
    print("\n" + "=" * 88)
    print(title)
    print("=" * 88)


def permanent_constant_table():
    """
    Exact-support-6 normalized K constants by ambient dimension.

    For six genuine nonzero coordinates and the remaining coordinates zero,
    |K(n)|*prod|n_i| depends only on the nonzero residues mod 7 and the number
    of zero coordinates.  This is the permanent/sector-cancellation constant
    hidden by the blunt 64*c1^6 envelope.
    """
    section("SUPPORT-6 RESIDUE PERMANENT CONSTANT")
    envelope = ENV64 * (C1 ** 6)
    print(f"blunt envelope constant 64*c1^6 = {envelope:.9g}")
    print(f"{'ambient d':>9} {'max residue constant':>22} {'best residues':>24} {'ratio to envelope':>18}")
    for d in [6, 7, 8, 9]:
        best_const = -1.0
        best_res = None
        for rs in itertools.product(range(1, 7), repeat=6):
            vec = list(rs) + [0] * (d - 6)
            const = abs(Kk(tuple(vec))) * prod(abs(x) for x in rs)
            if const > best_const:
                best_const = const
                best_res = rs
        print(f"{d:>9} {best_const:>22.9g} {str(best_res):>24} {best_const / envelope:>18.6g}")
    print(
        "Reading: the exact support-6 permanent improves the per-relation constant,\n"
        "but only by about 11x to 80x here.  The saved wide-spread proof still needs\n"
        "signed summation over the relation hyperplanes, not just a smaller absolute\n"
        "constant."
    )


def shape_report(name, E, cutoff, Hs, signed_Hs=None):
    section(f"SHAPE: {name}  E={E}")
    k = len(E)
    m = measS7(E)
    print(f"primitive={primitive(E)}  k={k}  span={max(E)}")
    print(f"meas(S7)={m} = {float(m):.6f}   M7={float(M7(k)):.6f}   corr={float(m - M7(k)):.6f}")
    print(f"large cutoff for relative/type-II supports: > {cutoff}")
    prev_all = 0.0
    prev_type = 0.0
    last_profile = None
    print("\nExact support-6 cumulative envelope (all supports vs supports touching a large offset)")
    print(f"{'H':>4} {'R6_all':>12} {'env6_all':>12} {'shell*H^2':>12} {'R6_type':>12} {'env6_type':>12} {'type_shell*H^2':>15}")
    for H in Hs:
        prof = support6_profile(E, H, cutoff)
        shell_all = prof["total_env"] - prev_all
        shell_type = prof["type_env"] - prev_type
        print(
            f"{H:>4} {prof['total_count']:>12} {fmt(prof['total_env']):>12} "
            f"{fmt(shell_all * H * H):>12} {prof['type_count']:>12} "
            f"{fmt(prof['type_env']):>12} {fmt(shell_type * H * H):>15}"
        )
        prev_all = prof["total_env"]
        prev_type = prof["type_env"]
        last_profile = prof
    if last_profile is None:
        return
    if signed_Hs:
        print("\nExact SIGNED support-6 layer (same cutoff, cumulative |n_i|<=H)")
        print(f"{'H':>4} {'Re signed6_all':>16} {'Abs signed6_all':>16} {'Re signed6_type':>17} {'Abs signed6_type':>17}")
        for H in signed_Hs:
            signed = signed_support6_profile(E, H, cutoff)
            print(
                f"{H:>4} {signed['total'].real:>16.8g} {abs(signed['total']):>16.8g} "
                f"{signed['type_total'].real:>17.8g} {abs(signed['type_total']):>17.8g}"
            )
        print("  top signed support-6 pieces at largest signed H:")
        for vals, val, is_type in signed["rows"][:5]:
            tag = "type-II" if is_type else "core"
            print(f"    {vals}: signed=({val.real:.6g}{val.imag:+.2g}j), abs={abs(val):.6g}, {tag}")
    print("\nTop exact-support-6 hyperplanes at largest H:")
    for vals, count, env, is_type in last_profile["rows"][:5]:
        tag = "type-II" if is_type else "core"
        print(f"  {vals}: count={count}, env={fmt(env)}, {tag}")
    tour = participation_tournament(last_profile)
    print("\nTournament Analysis fingerprint (vertices=offsets, gauge=support-6 participation)")
    print(f"  Hamiltonian tie path / mass order: {tour['ordered']}")
    print(f"  score histogram: {tour['score_hist']}")
    print(f"  edge flips against natural offset order: {tour['flips_against_offset_order']}")
    lows = low_height_relations(E, Hmax=4, max_support=min(len(E) - 1, 8))
    print("\nLow-height genuine support>=6 resonances (height <= 4):")
    if not lows:
        print("  none")
    else:
        for row in lows[:8]:
            print(
                "  support=%s size=%d H=%d shell_count=%d absK_shell=%.6g signedK_shell=(%.3g%+.3gj) example=%s"
                % (
                    row["support"],
                    row["support_size"],
                    row["height"],
                    row["count_on_shell"],
                    row["absK_shell"],
                    row["signedK_shell"].real,
                    row["signedK_shell"].imag,
                    row["examples"][0],
                )
            )


def relative_stranger_table(H=28, cutoff=16):
    section("RELATIVE STRANGER TABLE: core {0..6} plus one large M")
    print(
        "This isolates the type-II support-6 envelope.  It does not sum over possible M;\n"
        "for each fixed E it measures the relations that actually touch the large coordinate."
    )
    print(f"{'M':>6} {'meas':>10} {'type_R6':>12} {'type_env6':>12} {'top_type_support':>30} {'top_env':>10}")
    for M in [14, 21, 28, 35, 49, 68, 97, 140, 211, 421]:
        E = list(range(7)) + [M]
        prof = support6_profile(E, H, cutoff)
        top_type = next((row for row in prof["rows"] if row[3]), None)
        if top_type:
            vals, _count, top_env, _ = top_type
            support = str(vals)
        else:
            top_env = 0.0
            support = "-"
        print(
            f"{M:>6} {float(measS7(E)):>10.6f} {prof['type_count']:>12} "
            f"{fmt(prof['type_env']):>12} {support:>30} {fmt(top_env):>10}"
        )


def main():
    section("LRC(14) SUPPORT-6 RELATIVE COUNT SCOUT")
    print("Goal: execute the coupled support-6 count that the free product bound missed.")
    print(f"Envelope used: |K(n)| <= 64*c1^r/prod|n_j| with c1={C1}.")
    print("Nonzero coefficients divisible by 7 are excluded because chat(7m,T)=0.")
    permanent_constant_table()

    Hs = [2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24, 28]
    shape_report("bounded AP k=8", list(range(8)), cutoff=16, Hs=Hs, signed_Hs=[2, 3, 4, 6, 8, 10, 12])
    shape_report("one-stranger resonant k=8", list(range(7)) + [21], cutoff=16, Hs=Hs, signed_Hs=[2, 3, 4, 6, 8, 10, 12])
    shape_report("one-stranger high k=8", list(range(7)) + [211], cutoff=16, Hs=Hs, signed_Hs=[2, 3, 4, 6, 8, 10, 12])
    shape_report("k=9 sampled worst wide", list(range(8)) + [68], cutoff=15, Hs=[2, 3, 4, 5, 6, 8, 10, 12, 16, 20], signed_Hs=[2, 3, 4, 6, 8, 10])
    shape_report("k=10 sampled worst wide", list(range(9)) + [22], cutoff=13, Hs=[2, 3, 4, 5, 6, 8, 10, 12, 16], signed_Hs=[2, 3, 4, 6, 8])
    relative_stranger_table(H=28, cutoff=16)

    section("READING")
    print(
        "1. The support-6 hyperplane envelope is the right coupled object, but it is\n"
        "   still numerically far too large at proof-relevant heights.  The missing\n"
        "   gain is not merely 'Minkowski makes the sum finite'; it is the exact\n"
        "   signed T-sum inside K(n), visible in the signed support-6 layer."
    )
    print(
        "2. Full-lattice successive minima are too blunt for wide sets with dense cores:\n"
        "   core-only relations stay short forever.  The useful object is relative:\n"
        "   supports that touch the large coordinates, or equivalently the quotient by\n"
        "   the bounded-core relation lattice."
    )
    print(
        "3. Low-height support-7 resonances, e.g. 1+2+3+4+5+6=21 and\n"
        "   1+2+3+4+5+7=22, explain the safe wide spikes.  They are finite subset-sum\n"
        "   walls, not a harmonic tail.  This suggests a proof split: enumerate bounded\n"
        "   subset-sum resonance walls, then apply the relative hyperplane-shell count\n"
        "   outside them."
    )
    print(
        "4. THM-410 analogy: interval reversals count triangle creation by interior\n"
        "   points of a long edge.  Here a large offset creates support>=6 resonances\n"
        "   only when it lands in a bounded core subset-sum interval.  The long object\n"
        "   is not dangerous by length alone; it is dangerous only through how many\n"
        "   bounded witnesses lie inside it."
    )
    print(
        "5. Assumption challenged: tournament vertices need not be runners/arcs.  In\n"
        "   this scout the natural vertices are proof obligations (six-support\n"
        "   hyperplanes, then offsets for a participation fingerprint).  This preserves\n"
        "   the LRC predicate 'support-6 tail below the cap margin' and destroys phase\n"
        "   location data; it is a counting quotient, not a witness-time quotient."
    )


if __name__ == "__main__":
    main()
