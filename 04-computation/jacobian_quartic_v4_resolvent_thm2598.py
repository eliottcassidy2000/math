# -*- coding: utf-8 -*-
"""Exact companion for THM-2598: quartic V4 resolvent transfer boundary."""

import sympy as sp


failures = []


def require(label, condition):
    """Optimization-safe exact check."""
    ok = bool(condition)
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        failures.append(label)
        raise RuntimeError(f"failed exact check: {label}")


def zero(label, expression):
    require(label, sp.expand(expression) == 0)


T, W, U, Z = sp.symbols("T W U Z")
A, B, C, D, E = sp.symbols("A B C D E")
p, q, r, s, t = sp.symbols("p q r s t")

print("THM-2598 QUARTIC V4 RESOLVENT AUDIT")

# General integral resolvent and invariant depression.
f = A*T**4 + B*T**3 + C*T**2 + D*T + E
R = W**3 - C*W**2 + (B*D - 4*A*E)*W + 4*A*C*E - B**2*E - A*D**2
disc_f = sp.discriminant(f, T)
disc_R = sp.discriminant(R, W)
zero("general integral resolvent has identical discriminant", disc_R - disc_f)

I = C**2 - 3*B*D + 12*A*E
J = 2*C**3 - 9*B*C*D - 72*A*C*E + 27*B**2*E + 27*A*D**2
zero("resolvent depression is Z^3-I Z/3-J/27",
     R.subs(W, Z + C/sp.Integer(3)) - (Z**3 - I*Z/3 - J/27))
zero("universal cusp identity 27 Disc=4I^3-J^2", 27*disc_f - 4*I**3 + J**2)

# Root provenance under Vieta, including the root-difference proof.
x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4")
roots = (x1, x2, x3, x4)
e1 = sum(roots)
e2 = sum(roots[i]*roots[j] for i in range(4) for j in range(i + 1, 4))
e3 = sum(roots[i]*roots[j]*roots[k]
         for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4))
e4 = sp.prod(roots)
vieta = {B: -A*e1, C: A*e2, D: -A*e3, E: A*e4}
betas = (
    A*(x1*x2 + x3*x4),
    A*(x1*x3 + x2*x4),
    A*(x1*x4 + x2*x3),
)
root_poly = sp.prod(W - beta for beta in betas)
zero("resolvent roots are the three complementary pair products",
     root_poly - R.subs(vieta))
zero("pairing difference factors into two quartic differences",
     betas[0] - betas[1] - A*(x1 - x4)*(x2 - x3))

# Depressed squared-pairing resolvent and homogeneous leading drop.
fA = A*T**4 + p*T**2 + q*T + r
SA = U**3 + 2*p*U**2 + (p**2 - 4*A*r)*U - A*q**2
zero("homogeneous squared-pairing resolvent has quartic discriminant",
     sp.discriminant(SA, U) - sp.discriminant(fA, T))
zero("leading drop is U(U+p)^2", SA.subs(A, 0) - U*(U + p)**2)
disc_expansion = A*(4*p**3*(4*p*r - q**2)
                    + A*(-128*p**2*r**2 + 144*p*q**2*r - 27*q**4)
                    + 256*A**2*r**3)
zero("exact leading-factor discriminant expansion",
     sp.discriminant(fA, T) - disc_expansion)

# V4 reconstruction on a fully explicit nondegenerate quartic.
Sf = U**3 - 14*U**2 + 49*U - 36
zero("f+ and f- have the same squared-pairing resolvent",
     Sf - (U - 1)*(U - 4)*(U - 9))
f_plus = T**4 - 7*T**2 + 6*T
f_minus = T**4 - 7*T**2 - 6*T
zero("f+ root factorization", f_plus - T*(T - 1)*(T - 2)*(T + 3))
zero("f- is the marked reversal f+(-T)", f_minus - f_plus.subs(T, -T))

sign_triples = ((1, 2, -3), (1, -2, 3), (-1, 2, 3), (-1, -2, -3))
reconstructed = []
for a, b, c in sign_triples:
    require("sign triple has fixed squares and product -q",
            (a*a, b*b, c*c, a*b*c) == (1, 4, 9, -6))
    reconstructed.append(((a+b+c)//2, (a-b-c)//2,
                          (-a+b-c)//2, (-a-b+c)//2))
expected = ((0, 1, 2, -3), (1, 0, -3, 2),
            (2, -3, 0, 1), (-3, 2, 1, 0))
require("four reconstruction sections are exactly the V4 relabellings",
        tuple(reconstructed) == expected)
zero("first biquadratic quotient algebra is split",
     (U**3 - 20*U**2 + 96*U) - U*(U - 8)*(U - 12))
zero("second biquadratic quotient algebra is split",
     (U**3 - 28*U**2 + 160*U) - U*(U - 8)*(U - 20))

# Tame inertia ledger: action on four roots and on three pairings.
pairings = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def apply_pairing(perm, pairing):
    return frozenset(frozenset((perm[a], perm[b])) for a, b in pairing)


def cycle_count(perm):
    seen = set()
    count = 0
    for start in range(len(perm)):
        if start not in seen:
            count += 1
            cur = start
            while cur not in seen:
                seen.add(cur)
                cur = perm[cur]
    return count


representatives = {
    "identity": (0, 1, 2, 3),
    "transposition": (1, 0, 2, 3),
    "double_transposition": (1, 0, 3, 2),
    "three_cycle": (1, 2, 0, 3),
    "four_cycle": (1, 2, 3, 0),
}
expected_ledger = {
    "identity": (0, 0, 0),
    "transposition": (1, 1, 0),
    "double_transposition": (2, 0, 1),
    "three_cycle": (2, 2, 0),
    "four_cycle": (3, 1, 1),
}
for name, perm in representatives.items():
    induced = tuple(pairings.index(apply_pairing(perm, pairing)) for pairing in pairings)
    d4 = 4 - cycle_count(perm)
    d3 = 3 - cycle_count(induced)
    require(f"inertia ledger {name}: (d4,d3,index tax)",
            (d4, d3, (d4 - d3)//2) == expected_ledger[name])

# Sharp double-transposition family and generic S4 specialization.
fst = (T**2 - 1)**2 - t*(T + s)
Rst = W**3 + 2*W**2 + (4*s*t - 4)*W + 8*s*t - t**2 - 8
Dst = -t**2*(256*s**3*t - 256*s**2 - 288*s*t + 27*t**2 + 256)
zero("double-transposition family quartic discriminant", sp.discriminant(fst, T) - Dst)
zero("double-transposition family resolvent discriminant", sp.discriminant(Rst, W) - Dst)
zero("t=0 resolvent is one simple plus one double root",
     Rst.subs(t, 0) - (W + 2)**2*(W - 2))

special = sp.Poly(fst.subs({s: 2, t: 1}), T)
degrees_mod2 = sorted(poly.degree() for poly, exponent in sp.factor_list(special, modulus=2)[1]
                      for _ in range(exponent))
degrees_mod3 = sorted(poly.degree() for poly, exponent in sp.factor_list(special, modulus=3)[1]
                      for _ in range(exponent))
require("S4 witness has Frobenius type 4 mod 2", degrees_mod2 == [4])
require("S4 witness has Frobenius type 1+3 mod 3", degrees_mod3 == [1, 3])

# Minimal inertia and finite-cusp controls.
zero("four-cycle control resolvent", (W**3 + 4*t*W) - W*(W**2 + 4*t))
zero("four-cycle control discriminant", sp.discriminant(T**4 - t, T) + 256*t**3)
R3 = W**3 - 3*t*W - t - t**2
zero("three-cycle control discriminant",
     sp.discriminant((T**3 - t)*(T - 1), T) + 27*t**2*(t - 1)**2)
zero("three-cycle resolvent discriminant",
     sp.discriminant(R3, W) + 27*t**2*(t - 1)**2)
u, v = sp.symbols("u v")
zero("finite cusp control has common discriminant",
     sp.discriminant(T**4 + u*T + v, T) - sp.discriminant(W**3 - 4*v*W - u**2, W))
zero("finite cusp equation", sp.discriminant(T**4 + u*T + v, T) - (256*v**3 - 27*u**4))

print("FAILED CHECKS:", failures if failures else "NONE")
