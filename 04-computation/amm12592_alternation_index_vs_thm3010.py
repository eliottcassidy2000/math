#!/usr/bin/env python3
"""Test: is the alternation index the shared invariant between THM-3009's
extremal atom and THM-3010's maximal-alternation extremals?

THM-3009 sec 12.2 proposed the index

    I(P) = max_{j>=1} |Delta^j P(0)| / (2^j sup|P|)   in [0,1],

(j >= 1: at j = 0 the ratio is 1 for every P, so j = 0 carries no
information).  Our capacity bound (ARCH) is saturated exactly at
P(y) = +-(-1)^y, which has I = 1.

THM-3010 (klein) works with h_0=1,h_1,...; Newton ratios R_k =
h_k^2/(h_{k-1}h_{k+1}); circuit c_k = log(R_k/R_{k-1}); and proves that the
METALLIC recurrences a_k = n a_{k-1} + a_{k-2} attain MAXIMAL circuit
alternation, via Simson a_{k-1}a_{k+1} - a_k^2 = (-1)^k.

Two readings must be separated:
  (a) the SIGN sequence of the circuit;
  (b) the VALUE sequence of the circuit.
The conjecture is interesting only if it has content for (b).
"""
from math import comb, log
from fractions import Fraction


def newton_ratios(h):
    """R_k = h_k^2/(h_{k-1} h_{k+1}) as floats, for k = 1..len(h)-2."""
    return [h[k] * h[k] / (h[k - 1] * h[k + 1]) for k in range(1, len(h) - 1)]


def circuit(h):
    """c_k = log(R_k / R_{k-1})."""
    R = newton_ratios(h)
    return [log(R[i] / R[i - 1]) for i in range(1, len(R))]


def alternation_index(P):
    """max_{j>=1} |Delta^j P(0)| / (2^j sup|P|)."""
    P = list(P)
    s = max(abs(x) for x in P)
    if s == 0:
        return 0.0
    best = 0.0
    for j in range(1, len(P)):
        d = sum((-1) ** (j - i) * comb(j, i) * P[i] for i in range(j + 1))
        best = max(best, abs(d) / (2 ** j * s))
    return best


def sign_seq(v, eps=1e-14):
    return [(1 if x > eps else (-1 if x < -eps else 0)) for x in v]


def sign_changes(v):
    s = [x for x in sign_seq(v) if x != 0]
    return sum(1 for i in range(1, len(s)) if s[i] != s[i - 1])


def metallic(n, N):
    a = [0, 1]
    while len(a) < N:
        a.append(n * a[-1] + a[-2])
    return a[1:]                      # a_1, a_2, ... (positive)


if __name__ == "__main__":
    N = 14
    print("reference points")
    alt = [(-1) ** y for y in range(N)]
    con = [1] * N
    dec = [(-1) ** y / 2.0 ** y for y in range(N)]
    print(f"  P(y) = (-1)^y            I = {alternation_index(alt):.6f}   <- ARCH extremal")
    print(f"  P(y) = 1  (constant)     I = {alternation_index(con):.6f}")
    print(f"  P(y) = (-1)^y 2^-y       I = {alternation_index(dec):.6f}   <- alternating but decaying")
    print()

    print("THM-3010 ballot columns: circuit, its sign pattern, and I")
    cols = {
        "binom(2k,k)":        [comb(2 * k, k) for k in range(0, N)],
        "Catalan":            [comb(2 * k, k) // (k + 1) for k in range(0, N)],
        "binom(2k,k-1) [B]":  [comb(2 * k, k - 1) if k >= 1 else 0 for k in range(0, N + 1)],
        "binom(2k+1,k-1)":    [comb(2 * k + 1, k - 1) if k >= 1 else 0 for k in range(0, N + 1)],
    }
    for name, h in cols.items():
        hh = [x for x in h if x > 0]
        c = circuit(hh)
        print(f"  {name:20s} circuit {''.join('+' if x>0 else '-' for x in sign_seq(c))}"
              f"  changes={sign_changes(c)}  I(values)={alternation_index(c):.6f}")
    print()

    print("THM-3010 metallic recurrences (the maximal-alternation extremals)")
    for n in (1, 2, 3, 4):
        a = metallic(n, N + 3)
        c = circuit(a)
        sgn = sign_seq(c)
        print(f"  n={n} ({'golden silver bronze  n=4'.split()[n-1] if n<=3 else 'n=4'})"
              f"  circuit {''.join('+' if x>0 else '-' if x<0 else '0' for x in sgn)}"
              f"  changes={sign_changes(c)}")
        print(f"        I(sign seq) = {alternation_index([x for x in sgn]):.6f}"
              f"    I(values) = {alternation_index(c):.6f}")
    print()
    print("""VERDICT
  On SIGN sequences the conjecture is true but empty: "maximal circuit
  alternation" MEANS the sign sequence is +-(-1)^k, and any such sequence has
  I = 1 by definition.  No content.
  On VALUE sequences it is FALSE.  Simson gives R_k = 1/(1+(-1)^k/a_k^2), so
  the circuit alternates in sign but its magnitude decays like a_k^-2, i.e.
  geometrically; a geometrically decaying alternating sequence has I -> 1/2
  (the j=1 difference dominates and later differences cannot keep up).  Our
  ARCH extremal has CONSTANT magnitude, which is what pins I = 1.
  So the shared feature is sign-alternation only; the index does NOT transfer
  as a nontrivial common invariant.""")
