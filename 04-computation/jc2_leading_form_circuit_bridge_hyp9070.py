"""Corrected bridge test.  Repo convention (THM-3000 sec 1):
   h_k = e_k / C(d,k),   R_k = h_k^2/(h_{k-1} h_{k+1}),   c_k = log(R_k/R_{k-1}).
THM-3003: R_k = R_{d-k} for all 1<=k<=d-1  <=>  {r_i} = {mu/r_i} (reciprocal
up to scaling)  <=>  log-root measure symmetric about its mean.
(My earlier test wrongly used R_k = h_k/h_{k-1} and so failed; redo.)"""
from itertools import combinations
import sympy as sp


def e_sym(roots):
    d = len(roots)
    return [sp.Integer(1)] + [sp.expand(sum(sp.prod(c) for c in combinations(roots, k)))
                              for k in range(1, d + 1)]


def h_seq(roots):
    d = len(roots)
    e = e_sym(roots)
    return [sp.nsimplify(e[k] / sp.binomial(d, k)) for k in range(d + 1)]


def R_seq(roots):
    """R_k = h_k^2/(h_{k-1}h_{k+1}), k = 1..d-1"""
    h = h_seq(roots)
    d = len(roots)
    return [sp.simplify(h[k] ** 2 / (h[k - 1] * h[k + 1])) for k in range(1, d)]


def antipalindromic(roots):
    R = R_seq(roots)
    d = len(roots)
    return all(sp.simplify(R[k - 1] - R[d - k - 1]) == 0 for k in range(1, d))


def reciprocal_closed(roots):
    d = len(roots)
    mu = sp.simplify(sp.prod(roots) ** sp.Rational(2, d))
    tgt = sorted([sp.simplify(mu / r) for r in roots], key=lambda z: float(sp.N(z)))
    have = sorted(roots, key=lambda z: float(sp.N(z)))
    return all(sp.simplify(a - b) == 0 for a, b in zip(have, tgt))


def sign_changes(roots):
    """circuit word: signs of c_k = log(R_k) - log(R_{k-1}), k=2..d-1"""
    R = R_seq(roots)
    w = []
    for k in range(1, len(R)):
        diff = sp.simplify(R[k] - R[k - 1])
        w.append('+' if sp.N(diff) > 0 else ('-' if sp.N(diff) < 0 else '0'))
    return "".join(w), sum(1 for i in range(len(w) - 1) if w[i] != w[i + 1] and '0' not in (w[i], w[i + 1]))


tests = [
    ("reciprocal {2,1/2}", [sp.Integer(2), sp.Rational(1, 2)]),
    ("reciprocal {3,1/3,1}", [sp.Integer(3), sp.Rational(1, 3), sp.Integer(1)]),
    ("reciprocal {2,1/2,5,1/5}", [sp.Integer(2), sp.Rational(1, 2), sp.Integer(5), sp.Rational(1, 5)]),
    ("scaled recip mu=4 {4,1,8,1/2}", [sp.Integer(4), sp.Integer(1), sp.Integer(8), sp.Rational(1, 2)]),
    ("generic {1,2,4}", [sp.Integer(1), sp.Integer(2), sp.Integer(4)]),
    ("generic {1,2,5,7}", [sp.Integer(1), sp.Integer(2), sp.Integer(5), sp.Integer(7)]),
    ("generic {1,2,3,9}", [sp.Integer(1), sp.Integer(2), sp.Integer(3), sp.Integer(9)]),
]
print("corrected THM-3003 bridge:  reciprocal-closed <=> antipalindromic circuit")
ok = True
for nm, r in tests:
    rc, ap = reciprocal_closed(r), antipalindromic(r)
    agree = (rc == ap)
    ok &= agree
    print(f"  {nm:32s} recip={rc!s:5s} antipal={ap!s:5s} {'AGREE' if agree else '*** DISAGREE ***'}")
print("\nBRIDGE VERIFIED" if ok else "\nBRIDGE STILL FAILING")

print("\nTHM-3004 cluster law check: sign changes <= 2K-3 (K = # distinct roots)")
cases = [
    ("(n+1)^2(n+3)^2(n+8)  [THM-3004 witness]",
     [sp.Integer(1), sp.Integer(1), sp.Integer(3), sp.Integer(3), sp.Integer(8)]),
    ("(n+1)^2(n+3)^2(n+9)^2", [sp.Integer(1)] * 2 + [sp.Integer(3)] * 2 + [sp.Integer(9)] * 2),
    ("(n+1)^2(n+3)^2(n+9)^2(n+27)^2",
     [sp.Integer(1)] * 2 + [sp.Integer(3)] * 2 + [sp.Integer(9)] * 2 + [sp.Integer(27)] * 2),
]
for nm, r in cases:
    w, sc = sign_changes(r)
    K = len(set(r))
    print(f"  {nm:42s} d={len(r)} K={K} word={w} changes={sc} bound 2K-3={2*K-3}"
          f" {'OK' if sc <= 2*K-3 else '*** VIOLATION ***'}")
