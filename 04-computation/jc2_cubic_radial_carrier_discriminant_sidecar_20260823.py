#!/usr/bin/env python3
"""Independent exact sidecar for the THM-3771 cubic radial-carrier family.

The canonical theorem carries the proof.  This assertion-based companion
checks the universal chart/Jacobian identities, rational primitive, vertical
component-value polynomial, three smooth controls, and all five hostile
boundaries.
"""

import hashlib
import json

import sympy as sp


X, T, z, S = sp.symbols("X T z S")
r, h = sp.symbols("r h")
c = sp.symbols("c", nonzero=True)
a = sp.symbols("a", nonzero=True)
b, d = sp.symbols("b d")
u0, u1, u2, u3 = sp.symbols("u0 u1 u2 u3")


def require(condition, label):
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def jacobian(f, g):
    return sp.factor(sp.diff(f, X) * sp.diff(g, T)
                     - sp.diff(f, T) * sp.diff(g, X))


def profile(expr):
    return u0 + u1 * expr + u2 * expr**2 + u3 * expr**3


def carrier(expr):
    return a * expr**2 + b * expr + d


def concrete_data(u_expr, phi_expr, r_value):
    z_xt = X * T
    U = sp.expand(X * u_expr.subs(z, z_xt))
    W = sp.expand(U + 3 * z_xt + r_value)
    Q = sp.expand(U * phi_expr.subs(S, W))
    V = sp.factor((S - r_value)
                  * u_expr.subs(z, (S - r_value) / 3) * phi_expr)
    return U, W, Q, V


def gradient_at(Q, x_value, t_value):
    return (
        sp.factor(sp.diff(Q, X).subs({X: x_value, T: t_value})),
        sp.factor(sp.diff(Q, T).subs({X: x_value, T: t_value})),
    )


def squarefree(poly):
    polynomial = sp.Poly(sp.factor(poly), S, domain=sp.QQ)
    return sp.gcd(polynomial, polynomial.diff()).degree() == 0


def main():
    z_xt = X * T
    u_xt = profile(z_xt)
    U = sp.expand(X * u_xt)
    W = sp.expand(U + 3 * z_xt + r)
    phi_W = carrier(W)
    Q = sp.expand(U * phi_W)

    require(sp.factor(jacobian(U, W) - 3 * U) == 0,
            "log-canonical Jacobian J(U,W)=3U")
    require(sp.factor(jacobian(W, Q) + 3 * Q) == 0,
            "Hamiltonian identity J(W,Q)=-3Q")
    P0 = -c * W / (3 * Q)
    require(sp.factor(jacobian(P0, Q) - c) == 0,
            "rational mate J(-cW/(3Q),Q)=c")

    # Birational inverse in the abstract (U,W)-chart.
    uu, ww = sp.symbols("U W", nonzero=True)
    z_inverse = (ww - uu - r) / 3
    u_inverse = profile(z_inverse)
    x_inverse = uu / u_inverse
    t_inverse = z_inverse * u_inverse / uu
    require(sp.factor(x_inverse * t_inverse - z_inverse) == 0,
            "inverse recovers z")
    require(sp.factor(x_inverse * profile(x_inverse * t_inverse) - uu) == 0,
            "inverse recovers U")
    require(sp.factor(uu + 3 * x_inverse * t_inverse + r - ww) == 0,
            "inverse recovers W")

    # Universal gradient identities in the etale chart (X,z), X!=0.
    x_chart = sp.symbols("x", nonzero=True)
    u_z = profile(z)
    U_chart = x_chart * u_z
    W_chart = U_chart + 3 * z + r
    phi_chart = carrier(W_chart)
    phi_prime_chart = 2 * a * W_chart + b
    F = phi_chart + U_chart * phi_prime_chart
    Q_chart = sp.expand(U_chart * phi_chart)
    require(sp.factor(sp.diff(Q_chart, x_chart) - u_z * F) == 0,
            "chart Q_x=u(phi+U phi')")
    require(
        sp.factor(
            sp.diff(Q_chart, z)
            - (x_chart * sp.diff(u_z, z) * F + 3 * U_chart * phi_prime_chart)
        ) == 0,
        "chart Q_z=xu'F+3Uphi'",
    )

    V_general = sp.factor((S - r) * profile((S - r) / 3) * carrier(S))
    require(sp.degree(V_general, S) == 6,
            "generic cubic-profile vertical polynomial degree deg(u)+3")

    # Smooth squarefree controls of profile degrees zero, one, and two.
    controls = (
        ("constant", sp.Integer(1), S * (S - 2), sp.Integer(1),
         S * (S - 1) * (S - 2)),
        ("linear", z - 1, (S - 1) * (S - 2), sp.Integer(0),
         S * (S - 1) * (S - 2) * (S - 3) / 3),
        ("quadratic", (z - 1) * (z - 2), (S - 4) * (S - 5), sp.Integer(0),
         S * (S - 3) * (S - 4) * (S - 5) * (S - 6) / 9),
    )
    control_records = []
    for name, u_control, phi_control, r_control, expected_V in controls:
        _, Wc, Qc, Vc = concrete_data(u_control, phi_control, r_control)
        require(sp.factor(Vc - expected_V) == 0, (name, "vertical polynomial"))
        require(squarefree(Vc), (name, "squarefree control"))
        profile_degree = sp.degree(u_control, z)
        if profile_degree == 0:
            expected_discriminant = (
                u_control**4
                * sp.discriminant(phi_control, S)
                * phi_control.subs(S, r_control)**2
            )
        else:
            expected_discriminant = (
                sp.Rational(1, 3)**(profile_degree * (profile_degree + 3))
                * sp.discriminant(u_control, z)
                * sp.discriminant(phi_control, S)
                * u_control.subs(z, 0)**2
                * phi_control.subs(S, r_control)**2
                * sp.resultant(
                    u_control,
                    phi_control.subs(S, 3 * z + r_control),
                    z,
                )**2
            )
        require(sp.factor(sp.discriminant(Vc, S) - expected_discriminant) == 0,
                (name, "vertical discriminant factorization"))
        gradient_basis = sp.groebner(
            [sp.diff(Qc, X), sp.diff(Qc, T)], X, T, order="grevlex"
        )
        require(
            len(gradient_basis.polys) == 1
            and gradient_basis.polys[0].as_expr() == 1,
            (name, "gradient ideal is the unit ideal"),
        )
        Pc = -c * Wc / (3 * Qc)
        require(sp.factor(jacobian(Pc, Qc) - c) == 0,
                (name, "direct rational mate"))
        line = h - c * S / 3
        remainder = sp.rem(sp.Poly(line, S), sp.Poly(Vc, S)).as_expr()
        require(sp.factor(remainder - line) == 0,
                (name, "linear equalizer cannot be divisible by V"))
        control_records.append((name, str(sp.factor(Vc))))

    # Five hostile boundaries, each with an explicit critical point.
    hostiles = (
        ("u_at_zero", z, (S - 1) * (S - 2), 0, (0, 0)),
        ("phi_at_r", sp.Integer(1), S * (S - 1), 0, (0, 0)),
        ("u_repeated_root", (z - 1)**2, (S - 4) * (S - 5), 0, (1, 1)),
        ("u_phi_collision", z - 1, (S - 3) * (S - 5), 0, (1, 1)),
        ("phi_repeated_root", sp.Integer(1), (S - 1)**2, 0, (1, 0)),
    )
    hostile_records = []
    for name, u_hostile, phi_hostile, r_hostile, point in hostiles:
        _, _, Qh, Vh = concrete_data(u_hostile, phi_hostile, r_hostile)
        require(gradient_at(Qh, *point) == (0, 0),
                (name, "explicit critical point"))
        require(not squarefree(Vh), (name, "V is not squarefree"))
        hostile_records.append((name, point, str(sp.factor(Vh))))

    semantic = {
        "jacobian_convention": "f_X*g_T-f_T*g_X",
        "log_jacobian": "3*U",
        "inverse": {
            "z": "(W-U-r)/3",
            "X": "U/u(z)",
            "T": "z*u(z)/U",
        },
        "mate": "-c*W/(3*Q)",
        "torsor": "P0+k(Q)",
        "vertical": "(S-r)*u((S-r)/3)*phi(S)",
        "controls": control_records,
        "hostiles": hostile_records,
    }
    digest = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("JC_CUBIC_RADIAL_CARRIER_DISCRIMINANT_SIDECAR")
    print("status=FINITE-EXACT_INDEPENDENT_SIDECAR;THM-3771_PROVED_AUDITED;JC2_OPEN")
    print("jacobian_convention=J(P,Q)=P_X*Q_T-P_T*Q_X")
    print("chart=z=XT;U=X*u(z);W=U+3*z+r;J(U,W)=3*U")
    print("inverse=z=(W-U-r)/3;X=U/u(z);T=z*u(z)/U")
    print("target=Q=U*phi(W);deg(phi)=2")
    print("rational_mate=P0=-c*W/(3*Q);J(P0,Q)=c")
    print("rational_torsor=P0+k(Q)")
    print("vertical_polynomial=V(S)=(S-r)*u((S-r)/3)*phi(S)")
    print("smoothness=Q_smooth_iff_V_squarefree")
    print("discV=3^(-d*(d+3))*disc(u)*disc(phi)*u(0)^2*phi(r)^2*Res(u,phi(3z+r))^2_for_d>=1")
    print("principal_equalizer=V(S)_divides_(h-c*S/3);impossible_degV=deg(u)+3>=3")
    print("smooth_controls=" + ";".join(name for name, _ in control_records))
    print("hostile_boundaries=" + ";".join(name for name, _, _ in hostile_records))
    print("semantic_sha256=" + digest)
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
