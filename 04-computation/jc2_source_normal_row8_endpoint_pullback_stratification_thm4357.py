#!/usr/bin/env python3
"""THM-4357 primary exact certificate.

Pull back the late exact-weight-twelve endpoint gates through THM-4308's
finite row-eight source-normal response.  This proves no all-row lift,
seam entry, Keller-pair existence, or Jacobian conjecture consequence.
"""

from fractions import Fraction as F
import sys


sys.stdout.reconfigure(newline="\n")


def q(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


DELTA = F(896, 15)
THETA = F(512, 75)
K = F(-32, 5)
E = F(-1376, 135)
UP5 = F(-731648, 2025)
TERMINAL_TANGENT_DIM = 7


def zeta(phi):
    return -F(3, 2) * phi


def U(xi):
    return (F(475515904) - F(109350) * xi) / F(200475)


def W(phi, xi):
    return -(
        F(4343625) * phi * phi
        - F(17172000) * xi
        + F(143826305024)
    ) / F(4009500)


def Z(phi, eta, xi):
    return (
        F(12506118074368)
        - F(173745000) * phi * phi
        - F(195463125) * phi * eta
        - F(926883000) * xi
    ) / F(108256500)


XI_U = F(475515904, 109350)
XI_W = F(143826305024, 17172000)
XI_Z = F(12506118074368, 926883000)
SLOPE_U = -F(6, 11)
SLOPE_W = F(424, 99)
SLOPE_Z = -F(22886, 2673)


def eta_for_Z0(phi, xi):
    if phi == 0:
        raise ValueError("Z=0 graph formula requires Phi!=0")
    return (
        F(12506118074368)
        - F(173745000) * phi * phi
        - F(926883000) * xi
    ) / (F(195463125) * phi)


def point(phi=0, eta=0, xi=0, alpha=0, beta=0):
    return {
        "phi": F(phi),
        "eta": F(eta),
        "xi": F(xi),
        "alpha": F(alpha),
        "beta": F(beta),
    }


def values(p):
    phi, eta, xi = p["phi"], p["eta"], p["xi"]
    return {
        **p,
        "zeta": zeta(phi),
        "U": U(xi),
        "W": W(phi, xi),
        "Z": Z(phi, eta, xi),
    }


def e5(v):
    return v["alpha"] ** 2 - 4 * v["W"] * UP5


def e3(v):
    return v["eta"] ** 2 - 4 * E * v["W"]


def gate_4327_u(v):
    implication_5 = (UP5 == 0) or (e5(v) != 0)
    antecedent_3 = UP5 == 0 and v["alpha"] == 0 and DELTA == 0
    implication_3 = (not antecedent_3) or (e3(v) != 0)
    return v["U"] == 0 and v["W"] * v["Z"] != 0 and implication_5 and implication_3


def gate_4327_z(v):
    return v["Z"] == 0 and v["U"] * v["beta"] * v["zeta"] != 0


def gate_4334(v):
    return (
        v["Z"] == 0
        and v["beta"] == 0
        and v["U"] * v["W"] * v["zeta"] != 0
    )


def gate_4337(v):
    return (
        v["Z"] == 0
        and v["zeta"] == 0
        and v["U"] * v["beta"] != 0
    )


def gate_4339(v):
    return (
        v["Z"] == 0
        and v["beta"] == 0
        and v["zeta"] == 0
        and v["U"] * K * v["W"] * (v["U"] + v["W"]) != 0
    )


def gate_4340(v):
    return v["U"] == 0 and v["W"] * v["Z"] != 0


def matrix_rank(rows):
    rows = [[F(x) for x in row] for row in rows]
    if not rows:
        return 0
    nrows, ncols = len(rows), len(rows[0])
    rank = col = 0
    while rank < nrows and col < ncols:
        pivot = next((r for r in range(rank, nrows) if rows[r][col] != 0), None)
        if pivot is None:
            col += 1
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [x / scale for x in rows[rank]]
        for r in range(nrows):
            if r != rank and rows[r][col] != 0:
                scale = rows[r][col]
                rows[r] = [a - scale * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
        col += 1
    return rank


def grad_Z_numerator(p):
    return [
        -F(347490000) * p["phi"] - F(195463125) * p["eta"],
        -F(195463125) * p["phi"],
        -F(926883000),
        F(0),
        F(0),
    ]


def cubic_discriminant(k0, theta, xi, w):
    return (
        xi * xi * theta * theta
        - 4 * w * theta**3
        - 4 * xi**3 * k0
        - 27 * w * w * k0 * k0
        + 18 * w * xi * theta * k0
    )


GATE_TEXT = (
    ("THM-4327-U", "U=0, WZ!=0; u!=0 implies E5!=0; u=a=d=0 implies E3!=0"),
    ("THM-4327-Z", "Z=0, U*beta*zeta!=0"),
    ("THM-4334", "Z=beta_11=0, U*W*zeta_3!=0"),
    (
        "THM-4337",
        "Z=0, zeta_3=0, U*beta_11!=0; W and Lambda=U+W arbitrary; K=[y^2]H arbitrary",
    ),
    (
        "THM-4339",
        "Z=0, beta_11=0, zeta_3=0; U*K*W*Lambda!=0, Lambda=U+W",
    ),
    ("THM-4340", "U=0, WZ!=0"),
    ("THM-4342", "Z=beta_11=zeta_3=K=0, U*W*(U+W)!=0"),
    ("THM-4344", "Z=beta_11=zeta_3=W=0, U*K*xi_10!=0"),
    (
        "THM-4350",
        "Z=beta_11=zeta_3=W=xi_10=0, U*K!=0, (alpha,Theta)!=(0,0)",
    ),
    (
        "THM-4351",
        "Z=beta_11=zeta_3=W=xi_10=alpha=Theta=0, U*K!=0",
    ),
    ("THM-4352", "formal-local law; no endpoint coefficient gate"),
    ("THM-4353", "Z=beta_11=zeta_3=W=xi_10=K=0, U!=0"),
    ("THM-4354", "Z=beta_11=zeta_3=W=xi_10=U=0, K*alpha!=0"),
    (
        "THM-4355",
        "Z=beta_11=zeta_3=W=xi_10=U=alpha=0, K!=0",
    ),
    (
        "THM-4356",
        "Z=beta_11=zeta_3=W=xi_10=U=K=0, Delta=D=5696/105",
    ),
)


def main():
    checks = 0

    def ck(condition):
        nonlocal checks
        if not condition:
            raise AssertionError(f"failed primary exact check {checks + 1}")
        checks += 1

    # THM-4308 constants and inherited seam relation.
    ck(DELTA == F(896, 15))
    ck(THETA == F(512, 75))
    ck(K == F(2848, 45) - F(7, 6) * DELTA)
    ck(K != 0 and THETA != 0 and UP5 != 0)
    ck(DELTA != F(5696, 105))
    ck(TERMINAL_TANGENT_DIM == 7)

    # Exact zeta-zero trident, tested as polynomial identities at enough
    # values for these affine functions and at every distinguished root.
    tests = (F(-1), F(0), F(1), XI_U, XI_W, XI_Z)
    for test_xi in tests:
        ck(U(test_xi) == SLOPE_U * (test_xi - XI_U))
        ck(W(F(0), test_xi) == SLOPE_W * (test_xi - XI_W))
        ck(Z(F(0), F(317), test_xi) == SLOPE_Z * (test_xi - XI_Z))
    ck(len({XI_U, XI_W, XI_Z}) == 3)
    ck(all(x != 0 for x in (XI_U, XI_W, XI_Z)))

    root_U = values(point(xi=XI_U))
    root_W = values(point(xi=XI_W))
    root_Z = values(point(xi=XI_Z))
    ck(root_U["U"] == 0 and root_U["W"] * root_U["Z"] != 0)
    ck(root_W["W"] == 0 and root_W["U"] * root_W["Z"] != 0)
    ck(root_Z["Z"] == 0 and root_Z["U"] * root_Z["W"] != 0)
    ck(root_Z["U"] + root_Z["W"] != 0)

    # The requested and stronger unit-ideal residuals.
    origin = values(point())
    ck(origin["zeta"] == 0 and origin["xi"] == 0 and origin["W"] != 0)
    ck(origin["W"] == F(-35956576256, 1002375))
    ck(root_W["Z"] == F(634780696576, 14488875) != 0)

    # Every proper pair among zeta,W,xi is nonempty over the algebraically
    # closed field: the third pair is represented by a nonzero Phi^2.
    phi2_w_xi = -F(143826305024, 4343625)
    ck(phi2_w_xi != 0)
    ck(
        -(
            F(4343625) * phi2_w_xi
            - F(17172000) * 0
            + F(143826305024)
        )
        / F(4009500)
        == 0
    )

    # Exact rational witnesses, evaluated against the full verbatim gates.
    witness_U = values(point(xi=XI_U))
    eta_generic = eta_for_Z0(F(1), F(0))
    witness_4327_Z = values(point(phi=1, eta=eta_generic, beta=1))
    witness_4334 = values(point(phi=1, eta=eta_generic, beta=0))
    witness_4337 = values(point(xi=XI_Z, beta=1))
    witness_4339 = values(point(xi=XI_Z, beta=0))
    ck(gate_4327_u(witness_U))
    ck(gate_4327_z(witness_4327_Z))
    ck(gate_4334(witness_4334))
    ck(gate_4337(witness_4337))
    ck(gate_4339(witness_4339))
    ck(gate_4340(witness_U))

    # The generic Z=0 graph identity and its active open inequalities.
    for test_phi, test_xi in ((F(1), F(0)), (F(-2), F(7)), (F(3, 5), XI_W)):
        graph_point = values(
            point(phi=test_phi, eta=eta_for_Z0(test_phi, test_xi), xi=test_xi)
        )
        ck(graph_point["Z"] == 0)
    ck(witness_4327_Z["U"] * witness_4327_Z["beta"] * witness_4327_Z["zeta"] != 0)
    ck(witness_4334["U"] * witness_4334["W"] * witness_4334["zeta"] != 0)

    # Exact codimensions in A^5(Phi,eta,xi,alpha,beta).  Nonzero gate
    # conditions are principal opens and hence preserve these dimensions.
    row_UN = [F(0), F(0), F(-109350), F(0), F(0)]
    row_phi = [F(1), F(0), F(0), F(0), F(0)]
    row_beta = [F(0), F(0), F(0), F(0), F(1)]
    rank_data = (
        ("4327-U/4340", [row_UN], witness_U, 1, 4),
        ("4327-Z", [grad_Z_numerator(witness_4327_Z)], witness_4327_Z, 1, 4),
        (
            "4334",
            [grad_Z_numerator(witness_4334), row_beta],
            witness_4334,
            2,
            3,
        ),
        ("4337", [row_phi, grad_Z_numerator(witness_4337)], witness_4337, 2, 3),
        (
            "4339",
            [row_phi, grad_Z_numerator(witness_4339), row_beta],
            witness_4339,
            3,
            2,
        ),
    )
    for _, rows, _, expected_rank, expected_dim in rank_data:
        ck(matrix_rank(rows) == expected_rank)
        ck(5 - expected_rank == expected_dim)

    # The smallest surviving late slice is uniformly off both Lambda=0 and
    # the cubic collision divisor because eta,alpha do not enter either.
    disc_4339 = cubic_discriminant(K, THETA, XI_Z, root_Z["W"])
    expected_disc = F(
        483125535259306642688385993911761371136,
        7776331997934642451171875,
    )
    ck(disc_4339 == expected_disc != 0)
    ck(root_Z["U"] + root_Z["W"] == F(17651227389952, 1042743375))

    # Root-uniqueness hostiles.
    ck(U(XI_U + 1) == SLOPE_U)
    ck(W(F(0), XI_W + 1) == SLOPE_W)
    ck(Z(F(0), F(999), XI_Z + 1) == SLOPE_Z)

    print("THM-4357 PRIMARY: SOURCE-NORMAL ROW-EIGHT ENDPOINT PULLBACK")
    print(f"fixed Delta={q(DELTA)} Theta={q(THETA)} K={q(K)} u={q(UP5)}")
    print(
        "zeta=0 trident: "
        f"U={q(SLOPE_U)}(xi-xi_U), "
        f"W={q(SLOPE_W)}(xi-xi_W), "
        f"Z={q(SLOPE_Z)}(xi-xi_Z)"
    )
    print(f"xi_U={q(XI_U)}")
    print(f"xi_W={q(XI_W)}")
    print(f"xi_Z={q(XI_Z)}")
    print(
        "root differences="
        f"{q(XI_U-XI_W)},{q(XI_U-XI_Z)},{q(XI_W-XI_Z)}"
    )
    print(f"unit residual <zeta,xi,W>: W(0,0)={q(origin['W'])}")
    print(f"unit residual <Z,zeta,W>: Z(Phi=0,xi=xi_W)={q(root_W['Z'])}")
    print(f"pairwise control W=xi=0: Phi^2={q(phi2_w_xi)}")
    print(
        "S_4339 constants: "
        f"U={q(root_Z['U'])} W={q(root_Z['W'])} "
        f"Lambda={q(root_Z['U']+root_Z['W'])}"
    )
    print(f"S_4339 cubic discriminant={q(disc_4339)}")
    print("VERBATIM GATES")
    for name, gate in GATE_TEXT:
        print(f"{name}: {gate}")
    print("PULLBACK CLASSIFICATION")
    print("dimension convention: source pullback in A^5; full row-eight jet adds A^7")
    print("NONEMPTY 4327-U source-dim=4 jet-dim=11 witness=(0,0,xi_U,0,0)")
    print("NONEMPTY 4327-Z source-dim=4 jet-dim=11 witness=(1,eta_g,0,0,1)")
    print("NONEMPTY 4334 source-dim=3 jet-dim=10 witness=(1,eta_g,0,0,0)")
    print("NONEMPTY 4337 source-dim=3 jet-dim=10 witness=(0,0,xi_Z,0,1)")
    print("NONEMPTY 4339 source-dim=2 jet-dim=9 witness=(0,0,xi_Z,0,0)")
    print("NONEMPTY 4340 source-dim=4 jet-dim=11 witness=(0,0,xi_U,0,0)")
    print(f"eta_g={q(eta_generic)}")
    print("EMPTY 4342: K=-32/5 and Delta!=5696/105")
    print("EMPTY 4344: <Z,zeta,W>=<1>")
    print("EMPTY 4350: <zeta,xi,W>=<1>")
    print("EMPTY 4351: Theta=512/75 and <zeta,xi,W>=<1>")
    print("N/A 4352: formal-local law, no endpoint coefficient gate")
    print("EMPTY 4353: K=-32/5 and <zeta,xi,W>=<1>")
    print("EMPTY 4354: <zeta,xi,W>=<1>")
    print("EMPTY 4355: <zeta,xi,W>=<1>")
    print("EMPTY 4356: K=-32/5, Delta!=5696/105, and <zeta,xi,W>=<1>")
    print("SCOPE finite-row-eight-only; no all-row lift, seam entry, Keller pair, JC(2), or DC(2)")
    print(f"PASS checks={checks}")


if __name__ == "__main__":
    main()
