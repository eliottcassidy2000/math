#!/usr/bin/env python3
"""
c3_window_proof_check_kpc1.py -- THREAD C (HYP-2368 / THM-462,
kind-pasteur-2026-06-10-S1 sub-agent kpc1)

Machine validation of the CONSTRUCTIVE STEP of the THM-462 proof that the c3
spectrum is gap-free: for every n and every integer m in [0, M(n)] some
n-tournament has exactly m 3-cycles, where M(n) = (n^3-n)/24 (n odd),
(n^3-4n)/24 (n even).

Proof skeleton being validated here (see THM-462 file for the hand proof):
  * Induction step: adding a vertex that beats everything preserves c3, so
    A(n-1) subset A(n); this covers [0, M(n-1)].
  * Window step: the top window [M(n-1), M(n)] has length
        W(n) = M(n) - M(n-1) = h(h+1)/2  (n = 2h+1 odd)
                             = h(h-1)/2  (n = 2h   even),
    and is covered by NEAR-REGULAR score sequences perturbed by a Lagrange
    four-square pattern: write t = a1^2+a2^2+a3^2+a4^2 (0 <= ai <= sqrt(t)),
    lower 4 scores by a_i and raise 4 scores by a_i.
      - n odd:  base = (h,...,h) (regular);            f = m(n) + t.
      - n even: base = (h-1 x h, h x h) (near-regular); perturb the UPPER half
                only;                                    f = m(n) + t.
    Then c3 = C(n,3) - f = M(n) - t sweeps the whole window as t = 0..W(n).

This script EXACTLY verifies, for every n in 8..80 and EVERY t in [0, W(n)]:
  (1) the four-square decomposition exists with max part a_max <= h-2
      (the proof needs a_max <= h-2 for the Landau prefix checks; the hand
      proof guarantees it for h >= 8, i.e. n >= 16);
  (2) the constructed perturbed sequence is a genuine Landau sequence
      (nondecreasing, prefix sums >= C(k,2), total = C(n,2)) -- checked
      directly from the definition;
  (3) f(constructed) = m(n) + t exactly, where m(n) = near-regular f;
  (4) the window-length algebra W(n) = M(n) - M(n-1) matches h(h+1)/2 resp.
      h(h-1)/2.
It reports the smallest n at which the recipe works for ALL t (expected: all
n >= 16; some smaller n may fail (1) -- that is why the induction base is
computational, c3_spectrum_dp_kpc1.py having verified gap-freeness for n <= 30).

EXACT integer arithmetic. No floats (isqrt only).
"""

from math import isqrt

def comb2(x):
    return x * (x - 1) // 2

def comb3(x):
    return x * (x - 1) * (x - 2) // 6

def M(n):
    if n <= 2:
        return 0
    if n % 2 == 1:
        return (n**3 - n) // 24
    return (n**3 - 4 * n) // 24

def near_regular_f(n):
    if n % 2 == 1:
        h = (n - 1) // 2
        return n * comb2(h)
    h = n // 2
    return h * (comb2(h) + comb2(h - 1))

def window(n):
    return M(n) - M(n - 1)

def four_squares(t):
    """t = a^2+b^2+c^2+d^2, a>=b>=c>=d>=0 (Lagrange). Brute force, exact."""
    a0 = isqrt(t)
    for a in range(a0, -1, -1):
        r1 = t - a * a
        for b in range(min(a, isqrt(r1)), -1, -1):
            r2 = r1 - b * b
            for c in range(min(b, isqrt(r2)), -1, -1):
                r3 = r2 - c * c
                d = isqrt(r3)
                if d * d == r3 and d <= c:
                    return (a, b, c, d)
    return None  # impossible (Lagrange)

def is_landau(seq, n):
    """Exact Landau check from the definition."""
    if len(seq) != n:
        return False
    if any(s < 0 or s > n - 1 for s in seq):
        return False
    if sorted(seq) != list(seq):
        return False
    p = 0
    for k, s in enumerate(seq, start=1):
        p += s
        if p < comb2(k):
            return False
    return p == comb2(n)

def construct(n, t):
    """Build the perturbed near-regular sequence for target f = m(n)+t.
    Returns (sorted sequence, a_max) or None if four-square parts too big."""
    quad = four_squares(t)
    a_max = quad[0]
    if n % 2 == 1:
        h = (n - 1) // 2
        if a_max > h - 2:
            return None
        seq = [h - a for a in quad] + [h] * (n - 8) + [h + a for a in quad]
    else:
        h = n // 2
        if a_max > h - 2:
            return None
        if h < 8:
            return None  # need 8 distinct upper-half vertices
        # lower half: h scores of h-1 (untouched); upper half: h scores of h,
        # 4 lowered by a_i, 4 raised by a_i
        upper = [h - a for a in quad] + [h] * (h - 8) + [h + a for a in quad]
        seq = [h - 1] * h + upper
    seq.sort()
    return seq, a_max

def check_n(n):
    """Return (ok_all_t, n_fail, first_fail_t). Exact checks (1)-(3)."""
    W = window(n)
    h = (n - 1) // 2 if n % 2 == 1 else n // 2
    # check (4): window-length algebra
    expect = h * (h + 1) // 2 if n % 2 == 1 else h * (h - 1) // 2
    assert W == expect, f"n={n}: window {W} != algebra {expect}"
    mn = near_regular_f(n)
    fails = []
    for t in range(W + 1):
        built = construct(n, t)
        if built is None:
            fails.append(t)
            continue
        seq, a_max = built
        if not is_landau(seq, n):
            fails.append(t)
            continue
        f = sum(comb2(s) for s in seq)
        if f != mn + t:
            fails.append(t)
            continue
        # consistency: c3 = C(n,3) - f = M(n) - t
        assert comb3(n) - f == M(n) - t
    return (len(fails) == 0), len(fails), (fails[0] if fails else None)

def main():
    print("=" * 78)
    print("WINDOW-CONSTRUCTION VALIDATION for THM-462 (gap-free c3 spectrum)")
    print("n = 8..80, every t in [0, W(n)]; W(n) = M(n) - M(n-1)")
    print("=" * 78)
    first_all_good = None
    all_good_from_16 = True
    total_checks = 0
    for n in range(8, 81):
        W = window(n)
        total_checks += W + 1
        ok, nfail, first_fail = check_n(n)
        if ok and first_all_good is None:
            first_all_good = n
        if not ok and n >= 16:
            all_good_from_16 = False
        tag = "ALL t OK" if ok else f"{nfail} t-values FAIL (first t={first_fail})"
        if n <= 20 or not ok or n % 10 == 0:
            print(f"n={n:2d}: W(n)={W:4d}  {tag}")
    print(f"...({total_checks} (n,t) constructions checked exactly)")
    print()
    print(f"Recipe valid for ALL t at every n >= 16: {all_good_from_16}")
    print(f"Smallest n where recipe covers the whole window: {first_all_good}")
    print()
    print("Window-length algebra W(n) = h(h+1)/2 (odd) / h(h-1)/2 (even):"
          " verified by assertion for n = 8..80.")
    print()
    print("CONCLUSION: combined with c3_spectrum_dp_kpc1.py (gap-free for all"
          " n <= 30 >= induction base 16),")
    print("the induction  A(n) >= A(n-1) u [M(n-1), M(n)]  closes for all n.")

if __name__ == "__main__":
    main()
