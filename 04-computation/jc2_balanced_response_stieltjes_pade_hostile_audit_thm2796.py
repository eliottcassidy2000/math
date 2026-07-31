#!/usr/bin/env python3
"""Independent hostile controls for the balanced THM-2796 theorem."""

import ast
from pathlib import Path

from itertools import product
from math import factorial, gcd

import sympy as sp


def require(ok, message):
    if not bool(ok):
        raise RuntimeError(message)


x = sp.symbols("x")


def cancel(expr):
    return sp.cancel(sp.together(expr))


def root_sums(poly, top):
    poly = sp.Poly(poly, x, domain=sp.QQ).monic()
    degree = poly.degree()
    coeff = poly.all_coeffs()
    out = [sp.Integer(degree)]
    for k in range(1, top + 1):
        if k <= degree:
            value = sum(coeff[j] * out[k - j] for j in range(1, k))
            value += k * coeff[k]
        else:
            value = sum(
                coeff[j] * out[k - j] for j in range(1, degree + 1)
            )
        out.append(sp.expand(-value))
    return out


def check_packet(S, E, points, parts, v, g, label):
    S = sp.Poly(S, x, domain=sp.QQ).monic().as_expr()
    E = sp.Poly(E, x, domain=sp.QQ).monic().as_expr()
    T = sp.prod(x - b for b in points)
    D = sp.prod((x - b) ** p for b, p in zip(points, parts))
    V = sp.expand(v * S * D * T**2)
    G = cancel(g * E / (D * T))
    F = cancel(V * G**2)
    lam = sp.cancel(v * g)
    M = sp.expand(S * E * T)
    C = sp.cancel(2 / lam)
    r = sp.degree(M, x) - 1
    require(r >= 1, f"{label}: illegal r=0 in balanced packet")

    identity = sp.expand(
        2 * S * T * sp.diff(E, x)
        + E * T * sp.diff(S, x)
        - E * S * sum(
            p * cancel(T / (x - b)) for b, p in zip(points, parts)
        )
    )
    require(sp.simplify(identity - C) == 0, f"{label}: factor identity")
    require(cancel(2 * V * sp.diff(G, x) + sp.diff(V, x) * G - 2) == 0,
            f"{label}: response ODE")
    require(cancel(F - V * G**2) == 0, f"{label}: square potential")
    require(cancel(M * sp.diff(F, x) - C * F) == 0,
            f"{label}: exact Pade identity")

    sums_s = root_sums(S, r)
    sums_e = root_sums(E, r)
    moments = []
    for k in range(r + 1):
        moment = (
            sums_s[k] + 2 * sums_e[k]
            - sum(p * b**k for b, p in zip(points, parts))
        )
        moments.append(sp.simplify(moment))
    require(all(value == 0 for value in moments[:r]),
            f"{label}: moment indexing")
    require(sp.simplify(moments[r] - C) == 0,
            f"{label}: first active moment")

    defect = sp.Poly(sp.expand(V - v * M**2), x)
    require(defect.degree() == r + 2, f"{label}: defect degree")
    require(sp.simplify(defect.LC() - v * C / r) == 0,
            f"{label}: defect leading coefficient")
    return F, V, M, C


# Literal hostile to an unconstrained converse: identity (2) alone does not
# force balanced degree in the exceptional r=0 configuration.
S_bad = sp.Integer(1)
E_bad = sp.Integer(1)
T_bad = x
D_bad = x
C_bad = -1
v_bad = 1
g_bad = -2
V_bad = x**3
G_bad = -2 / x**2
F_bad = cancel(V_bad * G_bad**2)
identity_bad = sp.expand(
    2 * S_bad * T_bad * sp.diff(E_bad, x)
    + E_bad * T_bad * sp.diff(S_bad, x)
    - E_bad * S_bad
)
require(identity_bad == C_bad, "unbalanced converse witness identity")
require(cancel(2 * V_bad * sp.diff(G_bad, x)
               + sp.diff(V_bad, x) * G_bad - 2) == 0,
        "unbalanced converse witness ODE")
require(F_bad == 4 / x and sp.limit(F_bad, x, sp.oo) == 0,
        "unbalanced converse witness did not expose zero infinity value")


packet_count = 0
for N in range(1, 9):
    # e=0, h=1.
    F, V, _M, _C = check_packet(
        x**N - 1, 1, (0,), (N,), sp.Rational(1, N**2), 2 * N,
        f"cyclic-{N}",
    )
    require(cancel(F - 4 * (1 - x ** (-N))) == 0, "cyclic formula")
    packet_count += 1

    if N >= 2:
        # e=h=1.
        B = sp.expand(x**N - N * x + N - 1)
        Q, rem = sp.div(B, (x - 1) ** 2)
        require(rem == 0, "symmetric double root")
        F, V, _M, _C = check_packet(
            Q, x - 1, (0,), (N,),
            sp.Rational(4, N**2 * (N - 1) ** 2),
            sp.Rational(N * (N - 1), 2),
            f"symmetric-{N}",
        )
        require(cancel(F - B / x**N) == 0, "symmetric formula")
        packet_count += 1

        # e=1, h=2 chord chamber.
        for d in range(1, N // 2 + 1):
            b = N - d
            gamma = sp.Rational(d, N)
            denominator = x**d * (x - 1) ** b
            c = denominator.subs(x, gamma)
            Q, rem = sp.div(sp.expand(denominator - c), (x - gamma) ** 2)
            require(rem == 0, "chord double root")
            F, _V, _M, _C = check_packet(
                Q, x - gamma, (0, 1), (d, b),
                sp.cancel(4 / (c**2 * N**2)),
                sp.cancel(c * N / 2),
                f"chord-{N}-{d}",
            )
            require(cancel(F - (1 - c / denominator)) == 0, "chord formula")
            packet_count += 1


def compose(a, b):
    return tuple(a[b[i]] for i in range(len(a)))


def inverse(a):
    out = [0] * len(a)
    for i, value in enumerate(a):
        out[value] = i
    return tuple(out)


def cycle_type(a):
    unseen = set(range(len(a)))
    cycles = []
    while unseen:
        current = min(unseen)
        length = 0
        while current in unseen:
            unseen.remove(current)
            current = a[current]
            length += 1
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def closure(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        a = frontier.pop()
        for b in generators:
            c = compose(a, b)
            if c not in group:
                group.add(c)
                frontier.append(c)
    return group


def swap(N, a, b):
    out = list(range(N))
    out[a], out[b] = out[b], out[a]
    return tuple(out)


chord_controls = 0
for N in range(2, 11):
    cycle = tuple(list(range(1, N)) + [0])
    adjacent_group = closure((cycle, swap(N, 0, 1)))
    require(len(adjacent_group) == factorial(N), f"S_N failure N={N}")
    for d in range(1, N // 2 + 1):
        chord = swap(N, 0, d)
        pole_perm = inverse(compose(chord, cycle))
        require(cycle_type(pole_perm)
                == tuple(sorted((d, N - d), reverse=True)),
                f"chord pole partition N={N},d={d}")
        group = closure((cycle, chord))
        g = gcd(N, d)
        m = N // g
        require(len(group) == g * factorial(m) ** g,
                f"chord group order N={N},d={d}")
        require((len(group) == factorial(N)) == (g == 1),
                f"chord S_N boundary N={N},d={d}")
        chord_controls += 1


# Explicit first passport-duplication hostile at N=5.
sigma_one = (1, 2, 3, 0, 4)
sigma_zero_a = (0, 2, 1, 4, 3)  # (1 2)(3 4)
sigma_zero_b = (0, 4, 3, 2, 1)  # (1 4)(2 3)
for sigma_zero in (sigma_zero_a, sigma_zero_b):
    sigma_inf = inverse(compose(sigma_zero, sigma_one))
    require(cycle_type(sigma_zero) == (2, 2, 1), "N5 zero passport")
    require(cycle_type(sigma_inf) == (4, 1), "N5 pole passport")
    require(len(closure((sigma_zero, sigma_one))) == 20, "N5 group order")
centralizer_orbit_a = set()
power = tuple(range(5))
for _ in range(4):
    centralizer_orbit_a.add(compose(compose(power, sigma_zero_a), inverse(power)))
    power = compose(power, sigma_one)
require(sigma_zero_b not in centralizer_orbit_a,
        "N5 duplicate dessins became isomorphic")


# Sharper THM-2245 response identities.
A, U, G, kappa = sp.symbols("A U G kappa", nonzero=True)
T_spectral = A**2 / U**2
q = A / U
R = U * G
K = sp.cancel(R / q)
V_symbol = U**2
require(sp.simplify(A * K - V_symbol * G) == 0, "AK=VG")
require(sp.simplify(V_symbol * T_spectral * K**2 - (A * K) ** 2) == 0,
        "V T K^2=(AK)^2")

T0, K0, Tp, Kp, A0 = sp.symbols("T0 K0 Tp Kp A0", nonzero=True)
log_derivative = sp.cancel((Tp * K0**2 + 2 * T0 * K0 * Kp)
                           / (T0 * K0**2))
Kp_from_one_form = sp.cancel(
    (2 * kappa * T0 / A0 - K0 * Tp) / (2 * T0)
)
require(sp.simplify(log_derivative.subs(Kp, Kp_from_one_form)
                    - 2 * kappa / (A0 * K0)) == 0,
        "Keller one-form logarithmic derivative")

source = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source)
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require(assert_nodes == 0, "truth-bearing Python assert found")


print("THM-2796 BALANCED RESPONSE INDEPENDENT HOSTILE AUDIT")
print("converse_without_balance=REFUTED_BY_S=E=1,T=x,p=1,F=4/x")
print("required_converse_repair=retain_sum_p=s+2e_and_nonzero_infinity")
print(f"balanced_symbolic_packets={packet_count}")
print("moment_indexing=m_0_through_m_(r-1)_zero;m_r=C")
print("square_defect=degree_rplus2;LC=vC/r")
print(f"chord_monodromy_controls={chord_controls}")
print("N5_duplicate_witness=(12)(34)_vs_(14)(23);passport=(4,1);orders=20")
print("AK_identity=AK=VG=lambda*M")
print("response_log_derivative=2kappa/(AK)=C/M")
print("squared_gate=V*T*K^2=(AK)^2=lambda^2*M^2")
print("THM2245_scope=retrospective_degree14_control_only")
print("live_scope=derive_degree26_K_then_test_AK_before_spectral_normalization")
print(f"assert_nodes={assert_nodes}")
print("verdict=REPAIR_CONVERSE_QUANTIFIER_AND_SCOPE;CORE_ALGEBRA_PASS")
