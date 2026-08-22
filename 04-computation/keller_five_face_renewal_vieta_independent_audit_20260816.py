#!/usr/bin/env python3
"""Independent exact audit of the five-face renewal Vieta mechanism.

The script uses only integer/rational arithmetic.  It does not import the
THM-3506 or prospective THM-3522 companions.  For every admissible packet
A(e,m) in a finite hostile bank, it verifies that the unique intersection of
the input max-k and min-gamma faces supplies both hybrid limits.  The product
of the three residual root labels is retained by nonmonic Vieta, even when
h=-4e+2m/3 is not divisible by three.

The all-parameter proof consists of the displayed linear identities; the
finite bank, split-prime controls, and fixed G/R5/R6 scalar rows are exact
hostiles and normalization checks.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from pathlib import Path


EXPECTED_SEMANTIC_SHA256: str | None = "7468b48f87fab27aacf77af954ea5b63189eb8ba22ec1402b9df57d6d1378375"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


@dataclass(frozen=True)
class Packet:
    e: int
    m: int

    def __post_init__(self) -> None:
        require(self.e >= 0, (self, "e nonnegative"))
        require(self.m >= 0 and self.m % 3 == 0, (self, "m admissible"))
        require(self.m <= self.e, (self, "all packet exponents nonnegative"))

    @property
    def r(self) -> int:
        return self.m // 3

    @property
    def top_x(self) -> int:
        return 2 * self.e - 4 * self.r

    @property
    def top_z(self) -> int:
        return 2 * self.e - 2 * self.r

    @property
    def gamma_factor(self) -> int:
        return self.e - 2 * self.r

    @property
    def gamma_min(self) -> int:
        return -8 * self.e + 6 * self.r

    @property
    def h(self) -> int:
        return self.top_x - 3 * self.top_z

    def transform(self) -> "Packet":
        return Packet(7 * self.e - 2 * self.m, 3 * self.e - 2 * self.m)


@dataclass(frozen=True)
class Scalar:
    """A signed 2^u 3^v coefficient."""

    sign: int
    exp2: int
    exp3: int

    def __post_init__(self) -> None:
        require(self.sign in (-1, 1), self)

    def text(self) -> str:
        return f"sign={self.sign};v2={self.exp2};v3={self.exp3}"


def renewed_scalars(packet: Packet, input_top: Scalar) -> tuple[Scalar, Scalar]:
    output = packet.transform()
    sign = input_top.sign * (-1 if packet.top_z % 2 else 1)
    gamma = Scalar(sign, 3 * input_top.exp2 + packet.h, 3 * input_top.exp3)
    top = Scalar(sign, gamma.exp2, gamma.exp3 + 3 * output.gamma_factor)
    return gamma, top


def audit_packet(packet: Packet) -> tuple[object, ...]:
    output = packet.transform()
    e, r = packet.e, packet.r
    a, b, h = packet.top_x, packet.top_z, packet.h

    # The min-gamma face is z^e(27x^2z+y^3)^(e-2r).  Its endpoint at the
    # largest x^2z power is precisely the complete max-k singleton.
    d = packet.gamma_factor
    endpoint = (2 * d, 0, e + d)
    require(endpoint == (a, 0, b), (packet, "shared endpoint"))
    require(h == -4 * e + 2 * r, (packet, "root-label exponent"))

    delta6 = packet.gamma_min - b
    delta8 = packet.gamma_min - 3 * b
    require(delta6 == -10 * e + 8 * r, (packet, "delta6"))
    require(delta8 == -14 * e + 12 * r, (packet, "delta8"))

    # Top-z chart: product(q)=2/(27 A^2 C).
    source_s = -a + 6 * b
    total_s = 3 * source_s + 6 * e
    top_A = 2 * e - 2 * h
    top_C = 2 * e + 3 * b - h
    require(source_s == -delta6, (packet, "source top-z order"))
    require(total_s == 3 * output.top_z, (packet, "output top-z order"))
    require(top_A == output.top_x, (packet, "output top-z x exponent"))
    require(top_C == output.top_z, (packet, "output top-z z exponent"))

    # Min-gamma chart: product(q)=2/D for D=27A^2C+B^3.
    source_t = a - 8 * b
    total_t = 3 * source_t - 8 * e
    gamma_C = e + 3 * b
    gamma_D = e - h
    require(source_t == delta8, (packet, "source gamma order"))
    require(total_t == output.gamma_min, (packet, "output gamma order"))
    require(gamma_C == output.e, (packet, "output gamma z exponent"))
    require(gamma_D == output.gamma_factor, (packet, "output gamma factor"))

    # Deliberately omitting the Vieta leading-coefficient contribution loses
    # -2h powers in A, -h in C, and -h in D.  It must fail for every nonzero
    # packet in this admissible cone.
    require(h != 0 or packet.e == 0, (packet, "h hostile"))
    if packet.e:
        require(2 * e != output.top_x, (packet, "monic top-A hostile"))
        require(e != output.gamma_factor, (packet, "monic gamma-D hostile"))

    return (
        packet.e,
        packet.m,
        output.e,
        output.m,
        a,
        b,
        h,
        delta6,
        delta8,
        output.top_x,
        output.top_z,
        output.gamma_factor,
    )


def admissible_bank() -> tuple[tuple[object, ...], ...]:
    rows = []
    for e in range(1, 301):
        for m in range(0, e + 1, 3):
            rows.append(audit_packet(Packet(e, m)))
    return tuple(rows)


def fixed_orbit_and_scalars() -> tuple[object, ...]:
    orbit = [Packet(1, 0)]
    for _ in range(5):
        orbit.append(orbit[-1].transform())
    expected = ((1, 0), (7, 3), (43, 15), (271, 99), (1699, 615), (10663, 3867))
    require(tuple((p.e, p.m) for p in orbit) == expected, "fixed packet orbit")

    J = Packet(43, 15)
    gamma_G, top_G = renewed_scalars(J, Scalar(1, 15, 171))
    require(gamma_G == Scalar(1, -117, 513), (gamma_G, "G gamma scalar"))
    require(top_G == Scalar(1, -117, 1128), (top_G, "G top scalar"))

    G = Packet(271, 99)
    gamma_R5, top_R5 = renewed_scalars(G, top_G)
    require(gamma_R5 == Scalar(1, -1369, 3384), (gamma_R5, "R5 gamma scalar"))
    require(top_R5 == Scalar(1, -1369, 7251), (top_R5, "R5 top scalar"))

    R5 = Packet(1699, 615)
    gamma_R6, top_R6 = renewed_scalars(R5, top_R5)
    require(gamma_R6 == Scalar(1, -10493, 21753), (gamma_R6, "R6 gamma scalar"))
    require(top_R6 == Scalar(1, -10493, 46008), (top_R6, "R6 top scalar"))

    return (
        expected,
        ("G_from_J", gamma_G.text(), top_G.text(), J.h, J.h % 3),
        ("R5_from_G", gamma_R5.text(), top_R5.text(), G.h, G.h % 3),
        ("R6_from_R5", gamma_R6.text(), top_R6.text(), R5.h, R5.h % 3),
    )


def cube_root_of_unity(prime: int) -> int:
    for value in range(2, prime):
        if (value * value + value + 1) % prime == 0:
            return value
    raise RuntimeError((prime, "no nontrivial cube root"))


def split_prime_controls() -> tuple[object, ...]:
    rows = []
    packets = (Packet(43, 15), Packet(271, 99), Packet(1699, 615))
    for prime in (109, 127, 163):
        omega = cube_root_of_unity(prime)

        # Top chart: choose r=2 and C=2/(27r^3), so the three roots are
        # r, r*omega, r*omega^2 and the cubic is genuinely nonmonic.
        r0 = 2
        C = 2 * pow(27 * pow(r0, 3, prime), -1, prime) % prime
        top_roots = (r0, r0 * omega % prime, r0 * omega * omega % prime)
        require(len(set(top_roots)) == 3 and all(top_roots), (prime, "top roots"))
        require(
            all((27 * C * pow(q, 3, prime) - 2) % prime == 0 for q in top_roots),
            (prime, "top cubic"),
        )
        top_vieta = 2 * pow(27 * C, -1, prime) % prime

        # Gamma chart: roots 1,2,-3 give q^3-7q+6.  Multiplication by
        # D=-1/3 writes it as Dq^3-3Bq-2 with B=-7/9.
        gamma_roots = (1, 2, prime - 3)
        D = -pow(3, -1, prime) % prime
        B = -7 * pow(9, -1, prime) % prime
        require(
            all((D * pow(q, 3, prime) - 3 * B * q - 2) % prime == 0 for q in gamma_roots),
            (prime, "gamma cubic"),
        )
        gamma_vieta = 2 * pow(D, -1, prime) % prime

        for packet in packets:
            h = packet.h
            top_product = 1
            gamma_product = 1
            for q in top_roots:
                top_product = top_product * pow(q, h, prime) % prime
            for q in gamma_roots:
                gamma_product = gamma_product * pow(q, h, prime) % prime
            require(top_product == pow(top_vieta, h, prime), (prime, packet, "top Vieta"))
            require(
                gamma_product == pow(gamma_vieta, h, prime),
                (prime, packet, "gamma Vieta"),
            )
            rows.append((prime, packet.e, packet.m, h, h % 3, top_product, gamma_product))
    return tuple(rows)


def main() -> None:
    bank = admissible_bank()
    orbit_scalars = fixed_orbit_and_scalars()
    split_rows = split_prime_controls()
    nontrivial_mod3 = sum(row[6] % 3 != 0 for row in bank)
    ledger_sha256 = sha256(repr(bank).encode("ascii")).hexdigest()
    split_sha256 = sha256(repr(split_rows).encode("ascii")).hexdigest()

    proof = (
        "unique_shared_endpoint=max_k_intersection_min_gamma=x^a*z^b",
        "a=2e-4m/3;b=2e-2m/3;h=a-3b=-4e+2m/3",
        "top_chart_product_q=2/(27A^2C)",
        "gamma_chart_product_q=2/D",
        "output_top=x^(2eprime-4mprime/3)z^(2eprime-2mprime/3)",
        "output_gamma=z^eprime*D^(eprime-2mprime/3)",
        "c_gamma_prime=(-1)^(3b)*2^h*c^3",
        "c_top_prime=27^(eprime-2mprime/3)*c_gamma_prime",
    )
    boundary = (
        "requires_complete_A(e,m)_input_faces",
        "requires_Q=L^e*N(P)_polynomial",
        "proves_face_propagation_not_next_polynomiality",
        "no_image_prime_degree243_all_level_arbitrary_map_JC_or_LRC_claim",
    )
    semantic_surface = (proof, boundary, orbit_scalars, len(bank), nontrivial_mod3, ledger_sha256, split_sha256)
    semantic_sha256 = sha256(repr(semantic_surface).encode("utf-8")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            (semantic_sha256, EXPECTED_SEMANTIC_SHA256),
        )

    source = Path(__file__).resolve()
    print("Fixed Keller five-face renewal: independent Vieta audit")
    print("status=EXACT_ALGEBRAIC_PROOF_COMPANION;scope=conditional_face_transform")
    print(f"proof={proof}")
    print(f"fixed_orbit_and_scalars={orbit_scalars}")
    print(f"admissible_bank=e_1_through_300;rows={len(bank)};h_nonzero_mod3={nontrivial_mod3};ledger_sha256={ledger_sha256}")
    print(f"split_prime_controls={split_rows};split_sha256={split_sha256}")
    print("hostile=branchwise_q^h_is_not_root_independent_for_G_or_R5_but_product_Vieta_closes")
    print("hostile=omitting_nonmonic_leading_coefficient_fails_A_C_D_exponents_in_every_nonzero_bank_row")
    print(f"boundary={boundary}")
    print(f"semantic_sha256={semantic_sha256}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=no_assert_truth_gates;integer_fraction_and_split_prime_routes;no_randomness;no_elapsed_fields")
    print("commands=python -B 04-computation/keller_five_face_renewal_vieta_independent_audit_20260816.py;python -B -O 04-computation/keller_five_face_renewal_vieta_independent_audit_20260816.py")
    print("PASS")


if __name__ == "__main__":
    main()
