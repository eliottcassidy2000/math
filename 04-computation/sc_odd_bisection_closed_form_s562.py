"""
sc_odd_bisection_closed_form_s562.py — monad-researcher-2026-06-02-S562

GOAL: Close the open target from my own S560 handoff (item c):
  "SC at ODD n is NOT base-4 Burnside of m; the odd-n SC values come from the
   one-fixed-point partition filter, an open target for its own clean closed form."

RESULT (derived + verified here): a clean closed form for the ODD-n bisection of
the self-converse (self-complementary) tournament count SC(n), parallel to the
PROVEN even-n identity SC(2m) = A(m,4) (MISTAKE-049 corrected framing, THM-283).

------------------------------------------------------------------------------
DERIVATION
------------------------------------------------------------------------------
SC(n) is a Burnside fixed-point sum over the cycle types of an anti-automorphism
(THM-283): only cycle types with all parts ≡ 2 mod 4, plus EXACTLY ONE fixed
point (a part of size 1) iff n is odd, fix a self-converse tournament. The weight
of a cycle type μ with cycle lengths (L_1,...,L_k) is 2^{c(μ)} / z(μ) where
  c(μ) = Σ_i floor(L_i/2) + Σ_{i<j} gcd(L_i, L_j)   (orbits of unordered pairs)
  z(μ) = Π_L  L^{m_L} m_L!                          (centralizer order).

Parts ≡ 2 mod 4 are exactly the DOUBLED ODD numbers: {2, 6, 10, 14, ...} = {2a : a odd}.
So a contributing cycle type biject with an odd partition λ of m (the halved
parts), via parts {2λ_i}; for odd n=2m+1 we additionally append one fixed point.

Write, for an odd partition λ ⊢ m:
  ℓ(λ)  = number of parts
  G(λ)  = Σ_{i<j} gcd(λ_i, λ_j)                       (between-cycle pair orbits)
  c2(λ) = Σ_i (λ_i-1)/2 + G(λ) = (m-ℓ)/2 + G(λ)       (A000568 / A(m,q) exponent)
  z(λ)  = Π L^{m_L} m_L!

Doubling λ -> {2λ_i}:  z(2λ) = 2^{ℓ(λ)} z(λ);  inside-orbits Σ floor(2λ_i/2)=m;
between-orbits Σ gcd(2λ_i,2λ_j) = 2 G(λ).

EVEN n = 2m (no fixed point), μ = {2λ_i}:
  c(μ) = m + 2G,  z(μ) = 2^ℓ z(λ)
  weight = 2^{m+2G} / (2^ℓ z) = 2^{m-ℓ} 4^{G} / z = 4^{(m-ℓ)/2 + G} / z = 4^{c2(λ)}/z(λ)
  => SC(2m) = Σ_{odd λ⊢m} 4^{c2(λ)}/z(λ) = A(m,4).            [PROVEN, recovered here]

ODD n = 2m+1 (one fixed point of size 1), μ = {1} ∪ {2λ_i}:
  fixed point adds floor(1/2)=0 inside, and gcd(1, 2λ_i)=1 for each of the ℓ parts,
  so c gains exactly ℓ; z gains a factor 1!·1^1 = 1.
  c(μ) = m + ℓ + 2G,  z(μ) = 2^ℓ z(λ)
  weight = 2^{m+ℓ+2G}/(2^ℓ z) = 2^{m+2G}/z = 2^m 4^{G}/z
         = 2^{ℓ} · 4^{(m-ℓ)/2 + G} / z = 2^{ℓ(λ)} · 4^{c2(λ)} / z(λ).

  =>  SC(2m+1) = Σ_{odd λ⊢m} 2^{ℓ(λ)} · 4^{c2(λ)} / z(λ)        [NEW closed form]
             =  2^m · Σ_{odd λ⊢m} 4^{G(λ)} / z(λ)               (equivalent form)

INTERPRETATION: the odd-n bisection is the SAME sum as A(m,4) = SC(2m), but each
odd-partition term carries an extra factor 2^{(number of parts)} — the price of
the single fixed vertex, which pairs (gcd 1) with each of the ℓ cycles.

This script:
  (1) recomputes SC(n) directly from THM-283 (the validated engine),
  (2) verifies BOTH new odd-n formulas equal SC(2m+1) exactly (high m),
  (3) re-verifies SC(2m)=A(m,4) (extended),
  (4) prints the full SC sequence and its two bisections for OEIS lookup,
  (5) extends everything to n=120 exactly (float-free, Fraction-asserted-integer).
"""
from fractions import Fraction
from math import gcd, factorial
from collections import Counter

# ---- known cross-check values (repo canon / OEIS A002785) --------------------
KNOWN_SC = {
    2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176, 9: 2752, 10: 8784,
    11: 279968, 12: 1492288, 13: 95458560, 14: 872687552, 15: 111698291584,
    16: 1787154671104, 17: 457509297625088, 18: 13013584213369088,
    19: 6662951988432581120,
}

# ---- odd-partition machinery -------------------------------------------------
def odd_partitions(m):
    """All partitions of m into odd parts, as sorted tuples (ascending)."""
    allowed = [p for p in range(1, m + 1) if p % 2 == 1]
    out = []
    def rec(rem, idx, cur):
        if rem == 0:
            out.append(tuple(cur))
            return
        for k in range(idx, len(allowed)):
            p = allowed[k]
            if p > rem:
                break
            cur.append(p)
            rec(rem - p, k, cur)
            cur.pop()
    rec(m, 0, [])
    return out

def ell(lam):
    return len(lam)

def G_between(lam):
    g = 0
    for i in range(len(lam)):
        for j in range(i + 1, len(lam)):
            g += gcd(lam[i], lam[j])
    return g

def c2(lam):
    """A000568 / A(m,q) pair-orbit exponent of an odd partition."""
    c = 0
    for L in lam:
        c += L // 2                       # = (L-1)/2 for odd L
    c += G_between(lam)
    return c

def z_inv(lam):
    """Fraction(1, z(lambda))."""
    denom = 1
    for L, mult in Counter(lam).items():
        denom *= (L ** mult) * factorial(mult)
    return Fraction(1, denom)

# ---- the sequences -----------------------------------------------------------
def A_q(m, q):
    """q-deformed tournament count A(m,q) = Σ_{odd λ⊢m} q^{c2(λ)}/z(λ).
       A(m,2)=A000568(m); A(m,4)=SC(2m)."""
    tot = Fraction(0)
    for lam in odd_partitions(m):
        tot += z_inv(lam) * (q ** c2(lam))
    assert tot.denominator == 1, f"A({m},{q}) not integer"
    return tot.numerator

def SC_direct(n):
    """SC(n) directly from THM-283: parts ≡2 mod4, +one fixed point iff n odd."""
    if n < 2:
        return 1
    allowed = [p for p in range(2, n + 1) if p % 4 == 2]   # {2,6,10,14,...}
    extra = 1 if n % 2 == 1 else 0
    target = n - extra
    tot = Fraction(0)
    parts_lists = []
    def rec(rem, idx, cur):
        if rem == 0:
            parts_lists.append(tuple(cur) + (1,) * extra)
            return
        for k in range(idx, len(allowed)):
            p = allowed[k]
            if p > rem:
                break
            cur.append(p)
            rec(rem - p, k, cur)
            cur.pop()
    rec(target, 0, [])
    for lengths in parts_lists:
        # pair-orbit exponent with floor(L/2) inside (handles even parts)
        c = 0
        for L in lengths:
            c += L // 2
        for i in range(len(lengths)):
            for j in range(i + 1, len(lengths)):
                c += gcd(lengths[i], lengths[j])
        denom = 1
        for L, mult in Counter(lengths).items():
            denom *= (L ** mult) * factorial(mult)
        tot += Fraction(2 ** c, denom)
    assert tot.denominator == 1, f"SC({n}) not integer"
    return tot.numerator

def SC_even_formula(m):
    """SC(2m) = A(m,4)."""
    return A_q(m, 4)

def SC_odd_formula_A(m):
    """SC(2m+1) = Σ_{odd λ⊢m} 2^{ℓ(λ)} 4^{c2(λ)} / z(λ)   [NEW form A]."""
    tot = Fraction(0)
    for lam in odd_partitions(m):
        tot += z_inv(lam) * (2 ** ell(lam)) * (4 ** c2(lam))
    assert tot.denominator == 1, f"SC_odd_A({m}) not integer"
    return tot.numerator

def SC_odd_formula_B(m):
    """SC(2m+1) = 2^m Σ_{odd λ⊢m} 4^{G(λ)} / z(λ)         [NEW form B, equivalent]."""
    tot = Fraction(0)
    for lam in odd_partitions(m):
        tot += z_inv(lam) * (4 ** G_between(lam))
    val = (2 ** m) * tot
    assert val.denominator == 1, f"SC_odd_B({m}) not integer"
    return val.numerator

# ---- main --------------------------------------------------------------------
if __name__ == "__main__":
    import time
    print("=" * 74)
    print("SC ODD-n BISECTION CLOSED FORM  —  monad-researcher-2026-06-02-S562")
    print("=" * 74)

    # [1] reproduce known SC values
    print("\n[1] Direct THM-283 SC(n) vs known canon / OEIS A002785 (n=2..19):")
    ok = True
    for n in range(2, 20):
        v = SC_direct(n)
        good = (n not in KNOWN_SC) or (v == KNOWN_SC[n])
        ok &= good
        print(f"  SC({n:2d}) = {v}{'' if good else '   <-- MISMATCH!'}")
    print(f"  All known SC reproduced: {ok}")

    # [2] verify even-n identity SC(2m)=A(m,4), extended
    print("\n[2] EVEN-n identity  SC(2m) = A(m,4)  [proven, extended to m=40]:")
    even_ok = True
    for m in range(1, 41):
        lhs, rhs = SC_direct(2 * m), SC_even_formula(m)
        match = (lhs == rhs)
        even_ok &= match
        if m <= 6 or not match:
            print(f"  m={m:2d}: SC({2*m}) = {lhs}   A({m},4) = {rhs}   "
                  f"[{'OK' if match else 'MISMATCH'}]")
    print(f"  SC(2m) == A(m,4) for m=1..40: {even_ok}")

    # [3] verify NEW odd-n closed form (two equivalent expressions)
    print("\n[3] ODD-n closed form  SC(2m+1) = Σ 2^ℓ 4^c2 / z  =  2^m Σ 4^G / z")
    print("    [NEW — verified against direct THM-283 to m=40]:")
    odd_ok = True
    for m in range(0, 41):
        direct = SC_direct(2 * m + 1)
        fa = SC_odd_formula_A(m)
        fb = SC_odd_formula_B(m)
        match = (direct == fa == fb)
        odd_ok &= match
        if m <= 6 or not match:
            print(f"  m={m:2d}: SC({2*m+1}) = {direct}   formA = {fa}   formB = {fb}"
                  f"   [{'OK' if match else 'MISMATCH'}]")
    print(f"  SC(2m+1) == formA == formB for m=0..40: {odd_ok}")

    print("\n  ==> BOTH bisections now have matching base-4 Burnside closed forms:")
    print("        SC(2m)   = Σ_{odd λ⊢m}            4^{c2(λ)} / z(λ)")
    print("        SC(2m+1) = Σ_{odd λ⊢m} 2^{ℓ(λ)} · 4^{c2(λ)} / z(λ)")
    print("      (odd = even-sum with an extra 2^{#parts} fixed-point factor)")

    # [4] full SC sequence + bisections for OEIS lookup
    print("\n[4] Full SC(n) sequence, n=1..30 (for OEIS A002785 cross-check):")
    full = [SC_direct(n) for n in range(1, 31)]
    print("  " + ", ".join(str(x) for x in full))

    print("\n  EVEN bisection  SC(2),SC(4),...,SC(40)  (= A(m,4), m=1..20):")
    print("  " + ", ".join(str(SC_direct(2 * m)) for m in range(1, 21)))
    print("\n  ODD bisection   SC(1),SC(3),...,SC(41)  (m=0..20):")
    print("  " + ", ".join(str(SC_direct(2 * m + 1)) for m in range(0, 21)))

    # [5] exact extension to n=120
    print("\n[5] EXACT extension of SC(n) to n=120 (Fraction-asserted-integer):")
    t0 = time.time()
    for n in range(20, 121):
        v = SC_direct(n)
        if n % 2 == 1:
            m = (n - 1) // 2
            assert v == SC_odd_formula_A(m) == SC_odd_formula_B(m), f"n={n} formula break"
        else:
            m = n // 2
            assert v == SC_even_formula(m), f"n={n} even formula break"
        if n <= 24 or n % 20 == 0 or n == 120:
            print(f"  SC({n}) = {v}")
    print(f"  ... (all n=20..120 verified against the bisection closed forms; "
          f"{time.time()-t0:.2f}s)")

    print("\nDONE. Odd-n bisection closed form derived and verified — S560 handoff (c) closed.")
