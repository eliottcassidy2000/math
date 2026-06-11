#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
circulant_odd_det_macmini_s2.py — HYP-2388 (odd-function dictionary)
Instance: macmini-2026-06-10-S2 (subagent)

Circulant tournaments on Z_n (odd n) = odd ±1 functions f: Z_n -> {0,±1},
f(0)=0, f(-x)=-f(x); arc x->y iff f(y-x)=+1.  Connection set
S_set = {x : f(x)=+1}; tournament condition <=> oddness (S ∩ -S = ∅, per
MISTAKE-017 — asserted exhaustively below).

For each odd n in 3..23, enumerate all 2^((n-1)/2) circulant tournaments and:
 1. Verify det(I+S) = prod_{k=1}^{(n-1)/2} (1 + t_k^2), where i*t_k are the
    eigenvalues sum_x f(x) omega^{kx} (omega = e^{2 pi i / n}).  Integer det
    by fraction-free Bareiss on exact Python ints; spectral product at
    mpmath dps=50, relative-error check < 1e-30.
 2. Find ALL spectrally-flat circulants (all |t_k| equal).  Exact test:
    flat <=> det == (n+1)^((n-1)/2)  [AM-GM: sum(1+t_k^2) = (n-1)(n+1)/2,
    so prod <= ((n+1))^((n-1)/2) with equality iff all t_k^2 = n].
    Cross-checked exactly against the DRT condition SS^T = nI - J via the
    integer autocorrelation c_d = sum_u f(u) f(u+d) == -1 for all d != 0.
 3. Record the max-det circulant ("circulant Barba point") per n: max d,
    achieving connection sets up to rotation/multiplier equivalence, spectrum.
 4. Deficit vs the global odd-n ceiling (n+1)^((n-1)/2).
 5. Paley n = 7, 11, 19, 23: verify det(I+S) = (n+1)^((n-1)/2) exactly.
Also: H(T) (number of directed Hamiltonian paths, Redei) for every multiplier
orbit at n <= 15 via bitmask DP (vertex-transitivity: H = n * #paths starting
at vertex 0 — validated against the all-starts DP at n = 5, 7).

Run from repo root:
  python3 04-computation/circulant_odd_det_macmini_s2.py 2>&1 | \
      tee 05-knowledge/results/circulant_odd_det_macmini_s2.out
"""

import sys
import time
from math import gcd

import mpmath
mpmath.mp.dps = 50
REL_TOL = mpmath.mpf("1e-30")

T0 = time.time()


def log(msg=""):
    print(msg)
    sys.stdout.flush()


# ----------------------------------------------------------------------
# Exact integer determinant: fraction-free Bareiss with row pivoting.
# ----------------------------------------------------------------------
def bareiss_det(M):
    A = [row[:] for row in M]
    nn = len(A)
    sign = 1
    prev = 1
    for k in range(nn - 1):
        if A[k][k] == 0:
            piv = None
            for r in range(k + 1, nn):
                if A[r][k] != 0:
                    piv = r
                    break
            if piv is None:
                return 0
            A[k], A[piv] = A[piv], A[k]
            sign = -sign
        akk = A[k][k]
        Ak = A[k]
        for i in range(k + 1, nn):
            Ai = A[i]
            aik = Ai[k]
            for j in range(k + 1, nn):
                Ai[j] = (Ai[j] * akk - aik * Ak[j]) // prev
            Ai[k] = 0
        prev = akk
    return sign * A[nn - 1][nn - 1]


def build_I_plus_S(n, fvec):
    # (I+S)[i][j] = 1 on diagonal, f(j-i) off diagonal  (entries ±1)
    return [[1 if i == j else fvec[(j - i) % n] for j in range(n)]
            for i in range(n)]


# ----------------------------------------------------------------------
# Hamiltonian path counts (Redei H).  Circulants are vertex-transitive:
# H = n * (#Ham paths starting at vertex 0).  Validated below vs full DP.
# ----------------------------------------------------------------------
def in_masks(S, n):
    inm = [0] * n
    for v in range(n):
        mk = 0
        for u in range(n):
            if u != v and (v - u) % n in S:
                mk |= 1 << u
        inm[v] = mk
    return inm


def ham_paths_circ(S, n):
    """H(T) for circulant tournament with connection set S (start-0 trick)."""
    inm = in_masks(S, n)
    size = 1 << n
    dp = [None] * size
    base = [0] * n
    base[0] = 1
    dp[1] = base
    for mask in range(3, size, 2):
        cur = [0] * n
        vi = mask & (mask - 1)          # drop bit 0: end vertex v != 0
        while vi:
            vb = vi & -vi
            vi ^= vb
            v = vb.bit_length() - 1
            pm = mask ^ vb
            row = dp[pm]
            s = 0
            ui = pm & inm[v]
            while ui:
                ub = ui & -ui
                ui ^= ub
                s += row[ub.bit_length() - 1]
            cur[v] = s
        dp[mask] = cur
    return n * sum(dp[size - 1])


def ham_paths_full(S, n):
    """Reference all-starts DP (validation only, small n)."""
    inm = in_masks(S, n)
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        cur = dp[mask]
        vi = mask
        while vi:
            vb = vi & -vi
            vi ^= vb
            v = vb.bit_length() - 1
            pm = mask ^ vb
            row = dp[pm]
            s = 0
            ui = pm & inm[v]
            while ui:
                ub = ui & -ui
                ui ^= ub
                s += row[ub.bit_length() - 1]
            cur[v] = s
    return sum(dp[size - 1])


# ----------------------------------------------------------------------
# Multiplier (Adam) equivalence: S ~ m*S for m in U(n).  Rotations fix the
# connection set; m = -1 gives the reversal (x -> -x is an isomorphism
# T(S) ≅ T(-S) = T^op), so reversal is included automatically.
# ----------------------------------------------------------------------
def orbit_of(S, n, units):
    return {frozenset((u * x) % n for x in S) for u in units}


def canon(S, n, units):
    return min(tuple(sorted((u * x) % n for x in S)) for u in units)


def spectrum_charpoly_int(t2_list):
    """Integer coefficients (descending powers) of P(y) = prod_k (y - t_k^2).
    Integer because charpoly(S) = x * prod(x^2 + t_k^2) in Z[x]."""
    coeffs = [mpmath.mpf(1)]
    for v in t2_list:
        new = [mpmath.mpf(0)] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            new[i] += c
            new[i + 1] -= c * v
        coeffs = new
    ints = [int(mpmath.nint(c)) for c in coeffs]
    resid = max(abs(c - i) for c, i in zip(coeffs, ints))
    return ints, float(resid)


def qr_set(p):
    return frozenset((x * x) % p for x in range(1, p))


def fmt_t2(t2_list):
    return "[" + ", ".join(mpmath.nstr(v, 10) for v in t2_list) + "]"


# ----------------------------------------------------------------------
# Validation of the vertex-transitivity H trick.
# ----------------------------------------------------------------------
def validate_H_trick():
    log("[validation] H trick (H = n * #paths from vertex 0) vs full DP:")
    for n in (5, 7):
        m = (n - 1) // 2
        for bits in range(1 << m):
            S = frozenset((x if (bits >> (x - 1)) & 1 else n - x)
                          for x in range(1, m + 1))
            a = ham_paths_circ(S, n)
            b = ham_paths_full(S, n)
            assert a == b, (n, sorted(S), a, b)
        log(f"  n={n}: all {1 << m} circulants agree.")
    log("")


# ----------------------------------------------------------------------
# Main per-n processing.
# ----------------------------------------------------------------------
def process_n(n, do_H):
    m = (n - 1) // 2
    units = [u for u in range(1, n) if gcd(u, n) == 1]
    ceiling = (n + 1) ** m
    pow2 = 1 << (n - 1)
    two_pi = 2 * mpmath.pi
    sintab = [mpmath.sin(two_pi * j / n) for j in range(n)]

    all_sets = []
    for bits in range(1 << m):
        S = frozenset((x if (bits >> (x - 1)) & 1 else n - x)
                      for x in range(1, m + 1))
        all_sets.append(S)

    info = {}          # S -> (det, t2_list)
    flat_sets = []
    max_rel = mpmath.mpf(0)

    for S in all_sets:
        # MISTAKE-017 guard: oddness <=> tournament
        negS = {(n - x) % n for x in S}
        assert not (S & negS), f"S ∩ -S != ∅ at n={n}, S={sorted(S)}"
        assert (S | negS) == set(range(1, n)), f"S ∪ -S incomplete at n={n}"

        fvec = [0] * n
        for x in range(1, n):
            fvec[x] = 1 if x in S else -1

        det = bareiss_det(build_I_plus_S(n, fvec))
        assert det > 0, (n, sorted(S), det)
        assert det % pow2 == 0, f"det not divisible by 2^(n-1): n={n} S={sorted(S)}"

        # spectral product check (task 1)
        t2 = []
        for k in range(1, m + 1):
            s = mpmath.mpf(0)
            for x in S:
                s += sintab[(k * x) % n]
            t2.append((2 * s) ** 2)
        prod = mpmath.mpf(1)
        for v in t2:
            prod *= (1 + v)
        rel = abs(prod - det) / det
        if rel > max_rel:
            max_rel = rel
        assert rel < REL_TOL, f"spectral formula mismatch n={n} S={sorted(S)}"

        # exact DRT test via autocorrelation: SS^T = nI - J <=> c_d = -1 ∀d≠0
        is_drt = all(
            sum(fvec[u] * fvec[(u + d) % n] for u in range(n)) == -1
            for d in range(1, n))
        flat = (det == ceiling)
        assert is_drt == flat, f"flat<=>DRT broken at n={n} S={sorted(S)}"
        if flat:
            # numeric flatness double-check: all t_k^2 = n
            assert all(abs(v - n) < mpmath.mpf("1e-30") for v in t2)
            flat_sets.append(S)
        info[S] = (det, t2)

    # multiplier orbits
    orbits = {}       # canon tuple -> list of member sets
    for S in all_sets:
        orbits.setdefault(canon(S, n, units), []).append(S)
    # iso-invariance consistency: det constant on each orbit
    orb_rows = []
    for c, members in sorted(orbits.items()):
        dets = {info[S][0] for S in members}
        assert len(dets) == 1, f"det not orbit-constant at n={n} orbit {c}"
        det = dets.pop()
        rep = frozenset(c)
        orb_rows.append({
            "canon": c, "size": len(members), "det": det,
            "d": det // pow2, "t2": info[rep][1],
            "flat": det == ceiling,
        })
    orb_rows.sort(key=lambda r: (-r["det"], r["canon"]))

    # H per orbit (n <= 15)
    if do_H:
        for r in orb_rows:
            H = ham_paths_circ(frozenset(r["canon"]), n)
            assert H % 2 == 1, f"Redei violated?! n={n} orbit {r['canon']} H={H}"
            r["H"] = H

    max_det = orb_rows[0]["det"]
    min_det = min(r["det"] for r in orb_rows)
    max_orbits = [r for r in orb_rows if r["det"] == max_det]
    flat_orbits = sorted({canon(S, n, units) for S in flat_sets})

    # ---------------- output ----------------
    log(f"==== n={n}  (m={m}, circulants=2^{m}={1 << m}, "
        f"multiplier-orbits={len(orb_rows)}, |U(n)|={len(units)}) ====")
    log(f"  task1 spectral formula det(I+S) = prod(1+t_k^2): "
        f"all {1 << m} PASS exactly (max rel err {mpmath.nstr(max_rel, 3)}, dps=50)")
    log(f"  oddness/tournament check (S ∩ -S = ∅): all PASS")
    log(f"  ceiling (n+1)^m = {ceiling}")
    log(f"  flat circulants (all |t_k| equal): {len(flat_sets)} sets, "
        f"{len(flat_orbits)} orbit(s); verified flat <=> DRT (SS^T=nI-J) <=> "
        f"det==ceiling, exhaustively")
    if n % 4 == 3:
        if len(flat_sets) == 0:
            log(f"    NOTE: n={n} ≡ 3 (mod 4) but NO flat circulant exists.")
        for c in flat_orbits:
            tag = ""
            if all(gcd(x, n) == 1 for x in range(1, n)) or True:
                pass
            log(f"    flat orbit rep {c} "
                f"(orbit size {len(orbits[c])})")
        # QR comparison at primes p ≡ 3 mod 4
        if n in (3, 7, 11, 19, 23):
            QR = qr_set(n)
            cQR = canon(QR, n, units)
            only_qr = (flat_orbits == [cQR])
            log(f"    QR set = {tuple(sorted(QR))}, canon {cQR}; "
                f"NQR = -QR is in the SAME multiplier orbit "
                f"(-1 is a nonresidue for p≡3 mod 4)")
            log(f"    QR/NQR the ONLY flat circulants? {only_qr}")
    else:
        log(f"    (n={n} ≡ 1 mod 4: flat impossible — DRT needs n ≡ 3 mod 4; "
            f"confirmed: {len(flat_sets)} flat)")

    log(f"  MAX det = {max_det}   d = {max_det // pow2}   "
        f"deficit = ceiling - max = {ceiling - max_det}   "
        f"ratio max/ceiling = {mpmath.nstr(mpmath.mpf(max_det) / ceiling, 8)}")
    for r in max_orbits:
        cp, resid = spectrum_charpoly_int(r["t2"])
        log(f"    max orbit rep {r['canon']} (orbit size {r['size']}); "
            f"t_k^2 = {fmt_t2(r['t2'])}")
        log(f"      prod(y - t_k^2) integer coeffs (desc) = {cp}  "
            f"(rounding resid {resid:.1e})")
    log(f"  MIN det = {min_det}   d = {min_det // pow2}"
        + ("   [floor d=1 present]" if min_det == pow2 else ""))

    det_hist = {}
    for r in orb_rows:
        det_hist[r["det"]] = det_hist.get(r["det"], 0) + r["size"]
    log(f"  det histogram over all {1 << m} circulants "
        f"{{det: count}}: { {k: det_hist[k] for k in sorted(det_hist, reverse=True)} }")

    log(f"  orbit table (rep | size | det | d{' | H' if do_H else ''} | flat):")
    for r in orb_rows:
        h = f" | H={r['H']}" if do_H else ""
        log(f"    {r['canon']} | size {r['size']} | det {r['det']} | "
            f"d {r['d']}{h} | {'FLAT' if r['flat'] else '-'}")

    if do_H:
        maxH = max(r["H"] for r in orb_rows)
        maxH_orbits = [r["canon"] for r in orb_rows if r["H"] == maxH]
        maxdet_H = {r["canon"]: r["H"] for r in max_orbits}
        aligned = all(r["H"] == maxH for r in max_orbits)
        log(f"  H summary: max-H = {maxH} at orbit(s) {maxH_orbits}; "
            f"max-det orbit H values: {maxdet_H}; "
            f"max-det == max-H aligned? {aligned}")
    log(f"  [t = {time.time() - T0:.1f}s]")
    log("")

    out = {
        "n": n, "m": m, "ceiling": ceiling,
        "n_orbits": len(orb_rows),
        "max_det": max_det, "max_d": max_det // pow2,
        "min_det": min_det, "min_d": min_det // pow2,
        "deficit": ceiling - max_det,
        "flat_count": len(flat_sets), "flat_orbits": flat_orbits,
        "max_orbits": [r["canon"] for r in max_orbits],
        "max_rel_err": float(max_rel),
    }
    if do_H:
        out["maxH"] = maxH
        out["maxH_orbits"] = maxH_orbits
        out["maxdet_H"] = maxdet_H
        out["aligned"] = aligned
    return out


def main():
    log("=" * 78)
    log("HYP-2388 odd-function dictionary: circulant tournaments on Z_n, odd n=3..23")
    log("macmini-2026-06-10-S2 subagent | exact Bareiss dets + mpmath dps=50 spectra")
    log("=" * 78)
    log("")
    validate_H_trick()

    summary = []
    for n in range(3, 24, 2):
        summary.append(process_n(n, do_H=(n <= 15)))

    log("=" * 78)
    log("SUMMARY")
    log("=" * 78)
    log("")
    log("Task 1 — spectral determinant formula det(I+S) = prod_{k=1}^{(n-1)/2}(1+t_k^2):")
    log("  verified EXACTLY for all 4094 circulant tournaments, odd n = 3..23.")
    worst = max(s["max_rel_err"] for s in summary)
    log(f"  worst relative error across all n: {worst:.3e} (dps=50)")
    log("")
    log("Task 2 — spectrally flat circulants (= circulant DRTs, SS^T = nI - J):")
    log("  n | n mod 4 | #flat sets | #flat orbits | flat orbit reps")
    for s in summary:
        log(f"  {s['n']:2d} |    {s['n'] % 4}    | {s['flat_count']:4d}       | "
            f"{len(s['flat_orbits']):3d}          | {s['flat_orbits']}")
    log("")
    log("Task 3/4 — circulant Barba point (max det) vs global ceiling (n+1)^((n-1)/2):")
    log("  n | max det(I+S) | max d | ceiling | deficit | max/ceiling")
    for s in summary:
        ratio = mpmath.nstr(mpmath.mpf(s["max_det"]) / s["ceiling"], 8)
        log(f"  {s['n']:2d} | {s['max_det']} | {s['max_d']} | {s['ceiling']} | "
            f"{s['deficit']} | {ratio}")
    log("")
    log("  sequence circ-max det(I+S), n=3,5,...,23: "
        + str([s["max_det"] for s in summary]))
    log("  sequence circ-max d = det/2^(n-1),  same n: "
        + str([s["max_d"] for s in summary]))
    log("  sequence #multiplier-orbits of circulant tournaments: "
        + str([s["n_orbits"] for s in summary]))
    log("")
    log("Task 5 — Paley verification, n in {7, 11, 19, 23}:")
    for s in summary:
        n = s["n"]
        if n in (7, 11, 19, 23):
            ok = (s["max_det"] == s["ceiling"]
                  and canon(qr_set(n), n,
                            [u for u in range(1, n) if gcd(u, n) == 1])
                  in s["max_orbits"])
            log(f"  n={n}: Paley det(I+S) = {s['ceiling']} = "
                f"{n + 1}^{s['m']} exactly: {'VERIFIED' if ok else 'FAILED'} "
                f"(d = {s['ceiling'] // (1 << (n - 1))} = {(n + 1) // 4}^{s['m']})")
    log("")
    log("H alignment (n <= 15): circulant max-det vs circulant max-H")
    log("  n | max-det orbit(s): H | circ max-H | max-H orbit(s) | aligned?")
    for s in summary:
        if "maxH" in s:
            log(f"  {s['n']:2d} | {s['maxdet_H']} | {s['maxH']} | "
                f"{s['maxH_orbits']} | {s['aligned']}")
    log("")
    log(f"Total time: {time.time() - T0:.1f}s")


if __name__ == "__main__":
    main()
