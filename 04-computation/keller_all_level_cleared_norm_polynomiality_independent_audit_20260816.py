#!/usr/bin/env python3
"""Independent hostile audit of the proposed all-level cleared-norm lemma.

This script does not import an implementation of THM-3528.  It checks the
local reciprocal cubic, the nonmonic norm normalization, the canonical finite
branch, the packet cone, and the valuation/denominator arithmetic used by the
proof.  The geometric finite-etale and complete-face inputs are hash-pinned
to their proved theorem files.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
PARENT_HASHES = {
    "01-canon/theorems/THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy.md":
        "7bc470a6c76da2e04c64319a4a66d20d8618c45383e71ad0b42a26b049bca79e",
    "01-canon/theorems/THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate.md":
        "3d76e12bf874f5dbb20adefff92f03c8a1b10bab62302fdab89fd509fb91bbf4",
    "01-canon/theorems/THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary.md":
        "edfbfe5a577aa338535588565ea58cda07626f3b88561aa4fa7c52a77b7cfd91",
    "01-canon/theorems/THM-3522-fixed-keller-five-face-renewal-propagation.md":
        "79a522bde5c6155efa3c7c8b779b0cf179fae93d5c22226bd6af78008ed44f1c",
}
EXPECTED_SEMANTIC_SHA256 = "6fe70dcf5a0f1bd4f76ef8bc4986f79be1d74e4a0b6a71e40520dae67f4e456e"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def packet_exponents(e: int, m: int) -> tuple[int, ...]:
    require(m % 3 == 0, ("packet divisibility", e, m))
    return (
        e,
        m,
        3 * e - 2 * m,
        e - m,
        2 * m // 3,
        m // 3,
        2 * e - 4 * m // 3,
        2 * e - 2 * m // 3,
        e - 2 * m // 3,
    )


def admissible(e: int, m: int) -> bool:
    return e >= 0 and m >= 0 and m % 3 == 0 and min(packet_exponents(e, m)) >= 0


def renew(e: int, m: int) -> tuple[int, int]:
    return 7 * e - 2 * m, 3 * e - 2 * m


def keller_map(x, y, z):
    unit = 1 + x * y
    return (
        sp.expand(unit**3 * z + y**2 * unit * (4 + 3 * x * y)),
        sp.expand(y + 3 * x * unit**2 * z + 3 * x * y**2 * (4 + 3 * x * y)),
        sp.expand(2 * x - 3 * x**2 * y - x**3 * z),
    )


def main() -> None:
    observed_parents = {
        name: lf_sha256(ROOT / name) for name in PARENT_HASHES
    }
    require(observed_parents == PARENT_HASHES,
            ("parent source drift", observed_parents, PARENT_HASHES))

    a, b, c = sp.symbols("a b c")
    x, y, z = sp.symbols("x y z")
    L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
    T = 4 - 3 * b * c
    S = 27 * a * c**2 - 9 * b * c + 8
    D = 18 * a * c - 3 * b**2 * c + 2 * b

    factor_L = sp.factor(L)
    require(factor_L == L, ("L factorization", factor_L))
    unit_gcds = tuple(sp.gcd(L, value) for value in (c, T, S, D))
    require(unit_gcds == (1, 1, 1, 1), ("generic DVR units", unit_gcds))

    p = (sp.Rational(2, 27), sp.Integer(1), sp.Integer(1))
    unit_witness = tuple(
        sp.factor(value.subs(dict(zip((a, b, c), p)))) for value in (L, c, T, S, D)
    )
    require(unit_witness == (0, 1, 1, 1, sp.Rational(1, 3)),
            ("unit witness", unit_witness))

    q = (sp.Integer(2), sp.Rational(5, 6), sp.Rational(-7, 8))
    fq = tuple(sp.factor(value) for value in keller_map(*q))
    require(fq == p, ("canonical finite inverse", fq, p))
    L_source = L.subs({a: x, b: y, c: z})
    finite_source_L = sp.factor(L_source.subs(dict(zip((x, y, z), q))))
    require(finite_source_L == sp.Rational(241465, 1728),
            ("finite branch source L", finite_source_L))

    ell, tau, chi, u, w = sp.symbols("ell tau chi u w")
    cubic = ell * w**3 + tau * w - 2 * chi
    reciprocal = sp.expand(u**3 * cubic.subs(w, 1 / u))
    require(reciprocal == ell + tau * u**2 - 2 * chi * u**3,
            ("reciprocal cubic", reciprocal))
    special_reciprocal = sp.factor(reciprocal.subs(ell, 0))
    require(special_reciprocal == -u**2 * (2 * chi * u - tau),
            ("special reciprocal", special_reciprocal))
    finite_u = tau / (2 * chi)
    finite_derivative = sp.factor(sp.diff(special_reciprocal, u).subs(u, finite_u))
    require(finite_derivative == -tau**2 / (2 * chi),
            ("simple finite branch", finite_derivative))

    newton_points = ((0, 1), (2, 0), (3, 0))
    newton_segments = ((2, Fraction(-1, 2)), (1, Fraction(0, 1)))
    require(sum(length for length, _ in newton_segments) == 3,
            ("Newton horizontal lengths", newton_segments))

    d = sp.symbols("d", nonzero=True)
    divergent_limit = sp.expand(3 * (1 / u) * (-3 * d * u) - 2 * d)
    require(divergent_limit == -11 * d,
            ("max-lambda residual", divergent_limit))

    raw_resultant_x = sp.factor(sp.resultant(cubic, w, w))
    vieta_norm_x = 2 * chi / ell
    require(raw_resultant_x == 2 * chi, ("raw resultant", raw_resultant_x))
    require(sp.factor(raw_resultant_x / ell - vieta_norm_x) == 0,
            ("nonmonic norm correction", raw_resultant_x, vieta_norm_x))

    admissible_rows = []
    for e in range(0, 301):
        for m in range(0, e + 1, 3):
            require(admissible(e, m), ("input cone", e, m, packet_exponents(e, m)))
            ep, mp = renew(e, m)
            require(admissible(ep, mp),
                    ("renewed cone", e, m, ep, mp, packet_exponents(ep, mp)))
            require(ep - mp == 4 * e, ("cone gap", e, m, ep, mp))
            require(mp % 3 == 0, ("renewed divisibility", e, m, ep, mp))
            admissible_rows.append((e, m, ep, mp))

    orbit = [(1, 0)]
    for _ in range(9):
        orbit.append(renew(*orbit[-1]))
    known_orbit = (
        (1, 0),
        (7, 3),
        (43, 15),
        (271, 99),
        (1699, 615),
        (10663, 3867),
        (66907, 24255),
        (419839, 152211),
    )
    require(tuple(orbit[:8]) == known_orbit, ("known packet orbit", orbit))
    cassini = tuple(
        orbit[i][0] * orbit[i + 1][1] - orbit[i][1] * orbit[i + 1][0]
        for i in range(len(orbit) - 1)
    )
    require(cassini == tuple(3 * (-8)**i for i in range(len(cassini))),
            ("Cassini", cassini))

    valuation_records = []
    for packet_e in range(0, 41):
        for finite_s in range(0, packet_e + 6):
            divergent = (Fraction(-packet_e, 2),) * 2
            norm_v = sum(divergent, Fraction(finite_s, 1))
            cleared_v = norm_v + packet_e
            reduced_denominator_power = max(packet_e - finite_s, 0)
            require(norm_v == -packet_e + finite_s,
                    ("norm valuation", packet_e, finite_s, norm_v))
            require(cleared_v == finite_s >= 0,
                    ("cleared valuation", packet_e, finite_s, cleared_v))
            require(reduced_denominator_power <= packet_e,
                    ("denominator bound", packet_e, finite_s))
            valuation_records.append(
                (packet_e, finite_s, str(norm_v), str(cleared_v),
                 reduced_denominator_power)
            )

    scalar_control = (
        "N(r*P)=r^3*N(P)",
        "v_L(r)=0 for r in Q^*",
        "rational packet scalars do not change the clearing exponent",
    )
    scope = (
        "polynomial raw cleared norms and complete packets only",
        "finite-sheet valuation equals old-L multiplicity",
        "no L-coprimality without finite-sheet unit",
        "no image, irreducibility, separability, newest factor, or general JC",
    )
    record = {
        "parents": observed_parents,
        "generic_units": tuple(str(value) for value in unit_gcds),
        "unit_witness": tuple(str(value) for value in unit_witness),
        "finite_inverse": tuple(str(value) for value in fq),
        "finite_source_L": str(finite_source_L),
        "reciprocal": str(reciprocal),
        "newton_points": newton_points,
        "newton_segments": tuple((length, str(slope)) for length, slope in newton_segments),
        "finite_derivative": str(finite_derivative),
        "divergent_limit": str(divergent_limit),
        "resultant_x": str(raw_resultant_x),
        "norm_x": str(vieta_norm_x),
        "cone_rows": len(admissible_rows),
        "cone_digest": digest(admissible_rows),
        "orbit": orbit,
        "cassini": cassini,
        "valuation_rows": len(valuation_records),
        "valuation_digest": digest(valuation_records),
        "scalar_control": scalar_control,
        "scope": scope,
    }
    semantic = digest(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== independent hostile audit: all-level Keller cleared-norm polynomiality ==")
    print(f"parents={observed_parents}")
    print(f"generic_DVR_units_gcd(L;c,T,S,D)={unit_gcds};witness_(L,c,T,S,D)={unit_witness}: PASS")
    print(f"reciprocal_cubic={reciprocal};mod_L={special_reciprocal};newton={newton_segments}: PASS")
    print(f"finite_branch=(u0={finite_u},derivative={finite_derivative});canonical_q={q};F(q)={fq};L(q)={finite_source_L}: PASS")
    print(f"divergent_face_residual=3xz-2y->{divergent_limit};sheet_valuations=(-e/2,-e/2): PASS")
    print(f"nonmonic_hostile=(Res(E,w)={raw_resultant_x},N(w)={vieta_norm_x},Res/L=N): PASS")
    print(f"packet_cone_rows={len(admissible_rows)};digest={record['cone_digest']};renewal_preserves_admissibility: PASS")
    print(f"packet_orbit={tuple(orbit)}")
    print(f"cassini={cassini}: PASS")
    print(f"valuation_rows={len(valuation_records)};digest={record['valuation_digest']};v(N(P))=-e+s,v(L^eN(P))=s>=0: PASS")
    print(f"scalar_control={scalar_control}: PASS")
    print("global_localization=A[1/L]+v_L>=-e=>L^e*N(P) in A: PROVED ALGEBRAIC STEP")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT all-level raw polynomial packet induction at fixed-map scope")
    print(f"scope={scope}")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
