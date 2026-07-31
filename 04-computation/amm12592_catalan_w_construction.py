#!/usr/bin/env python3
"""Middle-weighted construction via the Catalan substitution w = u/(1+u).

REDUCTION.  Write E_k(u) = (1+u)^{a_k} Lam_k(w),  w = u/(1+u),  so that
1+u = 1/(1-w) and u = w/(1-w).  Then

    E_k(u)(1+u)^{L_k} = (1+u)^{a_k+L_k} Lam_k(w) = Lam_k(w) (1-w)^{-(m-1-k)}

because a_k + L_k = m-1-k, and u^{m-1} = w^{m-1}(1-w)^{-(m-1)}.  Multiplying
the whole system by (1-w)^{m-1}:

        sum_{k=0}^{m-1} Lam_k(w) (1-w)^k  =  eps w^{m-1},   deg Lam_k <= a_k.

THE L_k DISAPPEAR ENTIRELY.  (This is the Catalan/Riordan substitution; note
the Catalan GF C(w) = (1-sqrt(1-4w))/(2w) has C(-1) = 1/phi = delta* and
sqrt(1-4w)|_{w=-1} = sqrt5 = 1/p*, the two threshold constants.)

LEGAL ATOMS.  E_k is dominated by (1+u)^{a_k} iff [u^i]E_k = C(a_k,i)P_k(a_k-i)
with |P_k| <= 1.  Writing Lam_k = sum_c lam_c w^c gives
E_k = sum_c lam_c u^c (1+u)^{a_k-c}, and since 0 <= C(a_k-c, i-c) <= C(a_k,i),
a SUFFICIENT condition is that the positive coefficients of Lam_k sum to <= 1
and the negative ones to >= -1.  Over the integers that means exactly

        Lam_k  in  { 0,  +-w^c,  w^c - w^{c'} }      (c, c' <= a_k),

the middle-weighted unit atoms.

THE LADDER.  Lam_k(1) = R_k(1) is forced, and
R_{k+1}(1) = Lam'_k(1) - R'_k(1), so the algorithm is a carry ladder whose
"digit" is the exponent c (Lam' (1) = sigma c for a single atom, c - c' for a
difference).  Feasible iff every carry stays in {0,+-1} and every required
exponent fits in [0, a_k].
"""
import sys
from fractions import Fraction
from math import comb

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from amm12592_exact_block_profile_solver import profile


def trim(p):
    p = p[:]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def val1(p):
    return sum(p)


def der1(p):
    return sum(i * c for i, c in enumerate(p))


def div_1mw(p):
    """Divide by (1-w); requires p(1)=0.  (1-w)*q = p."""
    assert val1(p) == 0, "not divisible by (1-w)"
    n = len(p)
    q = [0] * n
    acc = 0
    for i in range(n - 1, 0, -1):
        acc -= p[i]
        q[i - 1] = acc
    # verify
    chk = [0] * (n + 1)
    for i, c in enumerate(q):
        chk[i] += c
        chk[i + 1] -= c
    assert trim(chk) == trim(p), (trim(chk), trim(p))
    return trim(q)


def run(m, C, eps=1, verbose=False):
    a = profile(m, C)
    if a is None:
        return "no profile", None
    R = [0] * m
    R[m - 1] = eps
    R = trim(R)
    Lams = []
    for k in range(m):
        s = val1(R)
        if abs(s) > 1:
            return f"carry blew up at k={k}: R_k(1)={s}", None
        d = der1(R)
        if s != 0:
            # single atom sigma w^c, sigma = s, Lam'(1) = s*c
            c = max(0, min(a[k], round(d / s)))
            if abs(s * c - d) > 1:
                return (f"exponent out of range at k={k}: need c~{d/s:.1f}, "
                        f"a_k={a[k]}"), None
            L = [0] * (c + 1)
            L[c] = s
        else:
            # difference atom w^c - w^{c'},  Lam'(1) = c - c'
            if abs(d) > a[k] + 1:
                return f"difference atom too small at k={k}: need {d}, a_k={a[k]}", None
            if d >= 0:
                c, cp = min(a[k], d), 0
            else:
                c, cp = 0, min(a[k], -d)
            if abs((c - cp) - d) > 1:
                return f"difference atom out of range at k={k}", None
            L = [0] * (max(c, cp) + 1)
            if c != cp:
                L[c] += 1
                L[cp] -= 1
        L = trim(L)
        if len(L) - 1 > a[k]:
            return f"degree overflow at k={k}", None
        Lams.append(L)
        S = R[:] + [0] * (len(L) - len(R))
        for i, cc in enumerate(L):
            S[i] -= cc
        S = trim(S)
        if val1(S) != 0:
            return f"not divisible at k={k}", None
        R = div_1mw(S) if S != [0] else [0]
        if R == [0]:
            Lams.extend([[0]] * (m - 1 - k))
            break
    if trim(R) != [0]:
        return "nonzero final residual", None
    return "OK", (a, Lams)


def verify(m, a, Lams, eps):
    """Check sum_k Lam_k (1-w)^k = eps w^{m-1} and the box condition."""
    acc = [0] * (2 * m + 2)
    for k, L in enumerate(Lams):
        pw = [0] * (k + 1)          # (1-w)^k
        for j in range(k + 1):
            pw[j] = comb(k, j) * (-1) ** j
        for i, c in enumerate(L):
            if c:
                for j, d in enumerate(pw):
                    acc[i + j] += c * d
    tgt = [0] * (2 * m + 2)
    tgt[m - 1] = eps
    if acc != tgt:
        return "identity FAILS"
    for k, L in enumerate(Lams):
        pos = sum(c for c in L if c > 0)
        neg = -sum(c for c in L if c < 0)
        if pos > 1 or neg > 1:
            return f"atom illegal at k={k}"
        if len(L) - 1 > a[k]:
            return f"degree overflow at k={k}"
    return "OK"


if __name__ == "__main__":
    print("carry ladder in the w = u/(1+u) coordinate")
    for m in (4, 8, 16, 32, 64, 128, 256):
        best = None
        for Cn, Cd in ((8, 5), (13, 8), (5, 3), (7, 4), (15, 8), (2, 1)):
            C = Fraction(Cn, Cd)
            for eps in (1, -1):
                st, res = run(m, C, eps)
                if st == "OK":
                    a, Lams = res
                    v = verify(m, a, Lams, eps)
                    if v == "OK":
                        best = (C, eps)
                        break
            if best:
                break
        if best:
            print(f"  m={m:4d}  smallest slope reached: C={best[0]} "
                  f"= {float(best[0]):.4f}  (eps={best[1]:+d})  VERIFIED")
        else:
            st, _ = run(m, Fraction(2, 1), 1)
            print(f"  m={m:4d}  none of the tested slopes worked; at C=2: {st}")
