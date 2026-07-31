#!/usr/bin/env python3
"""The shell-imbalance module: can two adjacent shells compensate?

Shell m is S_m = {b^m y : y in {0,1}^m, y != b^m}, words of length 2m.  A
colouring eps: S_m -> {+-1} has imbalance

    D_m(p) = sum_{w in S_m} eps_w p^{z(w)} q^{o(w)}
           = sum_k c_k p^{2m-k} q^k,    c_k = sum_{o(w)=k} eps_w.

Since {p^{2m-k}q^k} is the BERNSTEIN basis, D_m = 0 iff c = 0, and the
realisable c are exactly those with

    |c_k| <= N_k,   c_k = N_k (mod 2),   N_k = #{w in S_m : o(w) = k}.

N has generating function  (1+u^m)(1+u)^m - 1 - u^{2m}  (THM-3007).

Global fairness is sum_m D_m = 0; shell balance is every D_m = 0.  The gap
C*_1 - C* lives in the nonzero solutions.  ADJACENT COMPENSATION asks for
D_m + D_{2m} = 0 with both nonzero.  Re-expanding the degree-2m Bernstein
basis in the degree-4m one (multiply by (p+q)^{2m} = 1) turns this into the
exact convolution

    C^{(2m)}(u) = - C^{(m)}(u) (1+u)^{2m},                          (*)

so c^{(m)} DETERMINES c^{(2m)}, and the only question is whether the result
lands in its own box.  Two immediate consequences:

  * u = 1 in (*) with |c^{(2m)}| <= N^{(2m)} summed gives
    |C^{(m)}(1)| 2^{2m} <= 2(2^{2m}-1) < 2 * 2^{2m}, so |C^{(m)}(1)| < 2;
    as all c_k are even (dyadic m), C^{(m)}(1) = 0.
  * the MIDDLE class of shell 2m has only N_{2m} = 2 words, so
    |sum_j c_j^{(m)} binom(2m, j)| <= 2 -- a very tight constraint against a
    coefficient vector whose entries are even.
"""
from itertools import product
from math import comb


def N_vector(m):
    """N_k = #{w in S_m : o(w) = k}, k = 0..2m."""
    N = [0] * (2 * m + 1)
    for k in range(1, m + 1):
        N[k] += comb(m, k)                 # branch b = 0
    for k in range(m, 2 * m):
        N[k] += comb(m, k - m)             # branch b = 1
    return N


def convolve_box(c, m):
    """c^{(2m)} = -c^{(m)} * (1+u)^{2m}."""
    out = [0] * (4 * m + 1)
    for j, cj in enumerate(c):
        if cj:
            for t in range(2 * m + 1):
                out[j + t] -= cj * comb(2 * m, t)
    return out


def fits(c2, N2):
    if len(c2) > len(N2):
        if any(c2[i] for i in range(len(N2), len(c2))):
            return False
    for i in range(min(len(c2), len(N2))):
        if abs(c2[i]) > N2[i] or (c2[i] - N2[i]) % 2:
            return False
    return True


def search(m, verbose=True):
    """Exhaustive over c^{(m)} in its box with the right parity."""
    N1, N2 = N_vector(m), N_vector(2 * m)
    ranges = []
    for k in range(2 * m + 1):
        lo, hi, par = -N1[k], N1[k], N1[k] % 2
        ranges.append([v for v in range(lo, hi + 1) if (v - par) % 2 == 0])
    total, found = 0, []
    for c in product(*ranges):
        total += 1
        if not any(c):
            continue
        if sum(c) != 0:                    # forced by the u = 1 argument
            continue
        c2 = convolve_box(list(c), m)
        if fits(c2, N2):
            found.append((list(c), c2))
    if verbose:
        print(f"  m={m:3d}: searched {total} vectors, "
              f"{'NO nonzero compensating pair' if not found else f'FOUND {len(found)}'}")
        if found:
            for c, c2 in found[:3]:
                print(f"      c^(m)={c}  ->  c^(2m)={c2}")
    return found


def middle_obstruction(m):
    """The middle class of shell 2m has N_{2m} = 2; report the minimum of
    |sum_j c_j binom(2m,j)| over admissible nonzero c with sum c_j = 0."""
    N1 = N_vector(m)
    ranges = []
    for k in range(2 * m + 1):
        par = N1[k] % 2
        ranges.append([v for v in range(-N1[k], N1[k] + 1) if (v - par) % 2 == 0])
    best = None
    for c in product(*ranges):
        if not any(c) or sum(c) != 0:
            continue
        v = abs(sum(c[j] * comb(2 * m, j) for j in range(len(c))))
        if best is None or v < best[0]:
            best = (v, list(c))
    return best


if __name__ == "__main__":
    print("Can two ADJACENT shells compensate?  (shells m and 2m)")
    print()
    for m in (1, 2, 4):
        search(m)
    print()
    print("the middle-class bottleneck: |sum_j c_j binom(2m,j)| must be <= 2")
    for m in (1, 2, 4):
        b = middle_obstruction(m)
        if b is None:
            print(f"  m={m}: no admissible nonzero c with sum c_j = 0")
        else:
            print(f"  m={m}: minimum over admissible nonzero c is {b[0]}"
                  f"   (at c = {b[1]})   -> {'OK' if b[0] <= 2 else 'EXCEEDS 2'}")
    print()
    print("""READING.  The middle composition class of a shell always has exactly
TWO words (0^M 1^M and 1^M 0^M, THM-3007's cancelling corner pair), so the
compensating vector may put at most +-2 there.  But c^{(2m)} at the middle is
sum_j c_j binom(2m,j), a binomially WEIGHTED sum of an even, zero-sum vector,
and the central binomial weight is the largest in the row.  That is the
obstruction, and it is the same middle pair that makes ratio-2 blocks work in
the first place.""")
