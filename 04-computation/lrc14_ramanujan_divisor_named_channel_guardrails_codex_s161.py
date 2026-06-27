#!/usr/bin/env python3
"""HYP-2978: Ramanujan/divisor quotient guardrails for LRC14.

This script turns the divisor-function / Ramanujan-sum prompt into an
executable admissibility audit.

The guiding question is not "which scalar predicts LRC14?"  It is:

    what is a quotient allowed to forget before it stops preserving the
    proof predicate needed by the next implication?

Arithmetic seed:
  * sigma_k = id_k * 1 and phi = mu * id by Dirichlet convolution.
  * psi = id * |mu| keeps squarefree support with plus signs.
  * Jordan J_k = mu * id_k keeps k-dimensional unit capacity.
  * Ramanujan c_q(n) keeps the primitive q-th root shell:
        c_q(n) = sum_{d | gcd(q,n)} d mu(q/d).

LRC14 seed:
  * scalar qdiv, residue profiles, and raw divisor signatures are useful
    only if they do not collapse rows with different LRC routes.
  * the admissible quotient is therefore a labelled packet:
        q/Farey predicate + unit capacity + gcd profile + Ramanujan shell
        + Haar open/boundary state + endpoint/state-lift labels.

Tournament Analysis:
  Vertices are quotient channels, not runners.  Pair A -> B iff A preserves
  more of the declared LRC proof payload; ties follow the explicit Hamiltonian
  path in CHANNEL_ORDER.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, isqrt
from operator import mul
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
THRESHOLD = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s146 = load_module(
    "s161_haar_baire_boundary",
    REPO / "04-computation" / "lrc14_haar_baire_taut_boundary_s146.py",
)


def section(title: str) -> None:
    print("\n" + "=" * 96)
    print(title)
    print("=" * 96)


def prod(vals) -> int:
    return reduce(mul, vals, 1)


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors(n: int) -> list[int]:
    out: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def fmt_factor(n: int) -> str:
    fs = factor(n)
    if not fs:
        return "1"
    return "*".join(str(p) if e == 1 else f"{p}^{e}" for p, e in fs.items())


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def row_lcm(speeds: tuple[int, ...]) -> int:
    out = 1
    for v in speeds:
        out = lcm(out, v)
    return out


def radical(n: int) -> int:
    return prod(factor(n).keys()) if n > 1 else 1


def mobius(n: int) -> int:
    fs = factor(n)
    if any(e > 1 for e in fs.values()):
        return 0
    return -1 if len(fs) % 2 else 1


def phi(n: int) -> int:
    out = n
    for p in factor(n):
        out = out // p * (p - 1)
    return out


def sigma(n: int, k: int = 1) -> int:
    return sum(d**k for d in divisors(n))


def tau(n: int) -> int:
    return len(divisors(n))


def psi(n: int) -> int:
    return sum(d * abs(mobius(n // d)) for d in divisors(n))


def jordan(n: int, k: int) -> int:
    return sum(mobius(d) * (n // d) ** k for d in divisors(n))


def ramanujan_sum(q: int, n: int) -> int:
    return sum(d * mobius(q // d) for d in divisors(gcd(q, n)))


def circular_dist_num_mod(num: int, den: int) -> Fraction:
    r = num % den
    return Fraction(min(r, den - r), den)


def qdiv(speeds: tuple[int, ...], cap: int = 240) -> int:
    for q in range(2, cap + 1):
        if all(v % q for v in speeds):
            return q
    return cap + 1


def replace_one(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - {drop}) | {add}))


def replace_many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(drops)) | set(adds)))


@dataclass(frozen=True)
class Row:
    name: str
    route: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class RowAudit:
    row: Row
    q_threshold: int
    safe_mass: Fraction
    component_count: int
    lcm_value: int
    tau_lcm: int
    phi_lcm: int
    psi_lcm: int
    ram14: tuple[tuple[int, int], ...]
    ram27: tuple[tuple[int, int], ...]
    ram41: tuple[tuple[int, int], ...]
    safe_units_14: int
    safe_units_27: int
    safe_units_41: int

    @property
    def open_state(self) -> str:
        if self.safe_mass == 0:
            return "zero-open"
        return "positive-open"


def named_rows() -> list[Row]:
    return [
        Row("AP", "BOUNDARY-AP-GW", AP, "consecutive equality atom"),
        Row("GW 12->24", "BOUNDARY-AP-GW", replace_one(12, 24), "Goddyn-Wong equality atom"),
        Row("residue liar 12->26", "Q-WITNESS", replace_one(12, 26), "same mod-14 residues as AP, but q=12 witness"),
        Row("near 12->36", "K33-STATE-LIFT", replace_one(12, 36), "K33 / nonunit near miss"),
        Row("petal 10->20", "BOUNDARY-PETAL", replace_one(10, 20), "unit-petal positive row"),
        Row("petal 13->26", "BOUNDARY-PETAL", replace_one(13, 26), "unit-petal positive row"),
        Row(
            "P10+GW",
            "BOUNDARY-PETAL",
            replace_many((10, 12), (20, 24)),
            "unit-petal plus GW splice",
        ),
        Row(
            "P10+K33",
            "K33-STATE-LIFT",
            replace_many((10, 12), (20, 36)),
            "unit-petal plus K33 splice",
        ),
        Row("covering 12->84", "COVERING-MOMENT", replace_one(12, 84), "lcm-tail covering row"),
        Row("covering 12->168", "COVERING-MOMENT", replace_one(12, 168), "same lcm-tail family"),
        Row("covering 6->98", "COVERING-MOMENT", replace_one(6, 98), "HYP-2965 small covering mass"),
    ]


def ram_profile(speeds: tuple[int, ...], q: int) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(ramanujan_sum(q, v) for v in speeds).items()))


def safe_unit_count(speeds: tuple[int, ...], q: int) -> int:
    count = 0
    for a in range(1, q + 1):
        if gcd(a, q) != 1:
            continue
        if all(circular_dist_num_mod(v * a, q) >= THRESHOLD for v in speeds):
            count += 1
    return count


def audit_row(row: Row) -> RowAudit:
    comps = s146.safe_open_components(row.speeds)
    mass = s146.interval_measure(comps)
    L = row_lcm(row.speeds)
    return RowAudit(
        row=row,
        q_threshold=qdiv(row.speeds),
        safe_mass=mass,
        component_count=len(comps),
        lcm_value=L,
        tau_lcm=tau(L),
        phi_lcm=phi(L),
        psi_lcm=psi(L),
        ram14=ram_profile(row.speeds, 14),
        ram27=ram_profile(row.speeds, 27),
        ram41=ram_profile(row.speeds, 41),
        safe_units_14=safe_unit_count(row.speeds, 14),
        safe_units_27=safe_unit_count(row.speeds, 27),
        safe_units_41=safe_unit_count(row.speeds, 41),
    )


def verify_arithmetic_identities(limit: int) -> None:
    section("ARITHMETIC IDENTITIES AND WHAT THEY ARE ALLOWED TO FORGET")
    errors: dict[str, list[int]] = defaultdict(list)
    for n in range(1, limit + 1):
        if sigma(n, 1) != sum(d * 1 for d in divisors(n)):
            errors["sigma"].append(n)
        if sigma(n, 0) != tau(n):
            errors["tau=sigma0"].append(n)
        if sum(phi(d) for d in divisors(n)) != n:
            errors["sum phi"].append(n)
        if sum(mobius(d) * (n // d) for d in divisors(n)) != phi(n):
            errors["phi=mu*id"].append(n)
        if sum((n // d) * abs(mobius(d)) for d in divisors(n)) != psi(n):
            errors["psi=id*|mu|"].append(n)
        if sum(jordan(d, 2) for d in divisors(n)) != n * n:
            errors["sum J2"].append(n)
        for q in range(1, 25):
            direct = sum(
                complex(0)  # placeholder to keep direct formula out of float arithmetic
                for _ in ()
            )
            if direct != 0:
                raise AssertionError("unreachable")
            # Integer formula for Ramanujan sums avoids numerical roots of unity.
            c = ramanujan_sum(q, n)
            alt = mobius(q // gcd(q, n)) * phi(q) // phi(q // gcd(q, n))
            if c != alt:
                errors["ramanujan"].append(1000 * q + n)
    print(f"checked n<= {limit}, q<=24 for Ramanujan formula")
    for key in ("tau=sigma0", "sum phi", "phi=mu*id", "psi=id*|mu|", "sum J2", "ramanujan"):
        print(f"  {key:<14} mismatches={errors.get(key, [])[:8]}")
    print()
    print("Guardrail reading:")
    print("  tau/sigma forget residue order and phase; use them only for divisor-size capacity.")
    print("  phi/Jordan retain primitive unit capacity; use them for exact-period resource counts.")
    print("  psi retains squarefree plus-support; use it as a support capacity, not as a route label.")
    print("  Ramanujan c_q keeps the primitive cyclotomic shell but still forgets endpoint owners.")


def print_row_audit(audits: list[RowAudit]) -> None:
    section("NAMED LRC14 ROWS THROUGH DIVISOR AND RAMANUJAN QUOTIENTS")
    header = (
        f"{'row':<20} {'route':<18} {'qdiv':>4} {'safe_mu':>12} "
        f"{'lcm':>10} {'tau':>5} {'phi/lcm':>10} {'psi/lcm':>10} "
        f"{'U14':>4} {'U27':>4} {'U41':>4}"
    )
    print(header)
    for a in audits:
        print(
            f"{a.row.name:<20} {a.row.route:<18} {a.q_threshold:>4} {str(a.safe_mass):>12} "
            f"{fmt_factor(a.lcm_value):>10} {a.tau_lcm:>5} "
            f"{str(Fraction(a.phi_lcm, a.lcm_value)):>10} "
            f"{str(Fraction(a.psi_lcm, a.lcm_value)):>10} "
            f"{a.safe_units_14:>4} {a.safe_units_27:>4} {a.safe_units_41:>4}"
        )
    print()
    for q, attr in ((14, "ram14"), (27, "ram27"), (41, "ram41")):
        print(f"Ramanujan primitive-shell profiles c_{q}(v) over speeds:")
        for a in audits:
            print(f"  {a.row.name:<20} {getattr(a, attr)}")
        print()


def signature_collisions(audits: list[RowAudit]) -> None:
    section("QUOTIENT COLLISION AUDIT")
    channels = {
        "qdiv_only": lambda a: a.q_threshold,
        "open_state_only": lambda a: a.open_state,
        "mod14_residue_multiset": lambda a: tuple(sorted(v % 14 for v in a.row.speeds)),
        "ramanujan_14_profile": lambda a: a.ram14,
        "unit_counts_14_27_41": lambda a: (a.safe_units_14, a.safe_units_27, a.safe_units_41),
        "divisor_lcm_scalars": lambda a: (radical(a.lcm_value), a.tau_lcm, a.phi_lcm, a.psi_lcm),
        "guarded_packet_signature": lambda a: (
            a.q_threshold,
            a.open_state,
            a.ram14,
            a.ram27,
            a.ram41,
            a.safe_units_14,
            a.safe_units_27,
            a.safe_units_41,
            a.row.route,
        ),
    }
    for name, fn in channels.items():
        groups: dict[object, list[RowAudit]] = defaultdict(list)
        for a in audits:
            groups[fn(a)].append(a)
        bad = [
            group
            for group in groups.values()
            if len(group) > 1 and len({g.row.route for g in group}) > 1
        ]
        print(f"{name}: route-mixing collisions={len(bad)}")
        for group in bad[:5]:
            labels = ", ".join(f"{g.row.name}:{g.row.route}" for g in group)
            print(f"  collision -> {labels}")
    print()
    print("Readout: scalar channels are features, not proof quotients, whenever they mix routes.")
    print("The guarded packet signature is intentionally over-labelled; its job is to say what")
    print("has not yet been proved safe to forget.")


@dataclass(frozen=True)
class Channel:
    name: str
    keeps: tuple[str, ...]
    forgets: tuple[str, ...]
    vector: tuple[int, int, int, int, int, int, int]


CHANNELS = (
    Channel(
        "raw_divisor_counts",
        ("factor exponents", "tau/sigma capacity"),
        ("unit residues", "primitive phase", "endpoint owners", "route labels"),
        (1, 0, 0, 0, 0, 0, 0),
    ),
    Channel(
        "squarefree_psi_support",
        ("radical support", "plus squarefree capacity"),
        ("prime powers", "unit residues", "endpoint owners", "signs"),
        (1, 1, 0, 0, 0, 0, 0),
    ),
    Channel(
        "totient_jordan_unit_capacity",
        ("primitive exact-period counts", "CRT unit capacity"),
        ("unit arrangement", "endpoint owners", "K33/petal route labels"),
        (1, 1, 1, 0, 0, 0, 0),
    ),
    Channel(
        "gcd_strata",
        ("which divisors of q each speed hits", "imprimitive/unit split"),
        ("cyclic phase inside each stratum", "endpoint owners"),
        (1, 1, 1, 1, 0, 0, 0),
    ),
    Channel(
        "ramanujan_primitive_shell",
        ("primitive cyclotomic phase by q", "orthogonal shell projection"),
        ("individual primitive residues", "endpoint owners", "state-lift debt"),
        (1, 1, 1, 1, 1, 0, 0),
    ),
    Channel(
        "exact_period_packet",
        ("q/Farey predicate", "unit residues", "Ramanujan shell"),
        ("continuous endpoint order unless attached", "proof-family names"),
        (1, 1, 1, 1, 1, 1, 0),
    ),
    Channel(
        "labelled_lrc_packet_sheaf",
        ("q/Farey", "Haar open/boundary", "endpoint owners", "C27/K33/state labels", "dual certificates"),
        ("raw runner names when irrelevant",),
        (1, 1, 1, 1, 1, 1, 1),
    ),
)
CHANNEL_ORDER = tuple(c.name for c in CHANNELS)


def tournament_analysis() -> None:
    section("TOURNAMENT ANALYSIS ON QUOTIENT CHANNELS")
    names = [c.name for c in CHANNELS]
    order_index = {name: i for i, name in enumerate(CHANNEL_ORDER)}
    score = Counter({name: 0 for name in names})
    c3 = 0
    edges: dict[tuple[str, str], str] = {}
    for a, b in combinations(CHANNELS, 2):
        av = sum(a.vector)
        bv = sum(b.vector)
        if av > bv or (av == bv and order_index[a.name] < order_index[b.name]):
            winner, loser = a.name, b.name
        else:
            winner, loser = b.name, a.name
        score[winner] += 1
        edges[(winner, loser)] = "retains_more_payload"

    for x, y, z in combinations(names, 3):
        xy = (x, y) in edges
        yz = (y, z) in edges
        zx = (z, x) in edges
        if xy and yz and zx:
            c3 += 1
        yx = (y, x) in edges
        zy = (z, y) in edges
        xz = (x, z) in edges
        if yx and zy and xz:
            c3 += 1

    print("vertices:")
    for c in CHANNELS:
        print(f"  {c.name}: keeps={'; '.join(c.keeps)}")
        print(f"    forgets={'; '.join(c.forgets)}")
    print()
    print(f"score_hist={dict(Counter(score.values()))}")
    print(f"directed_3_cycles={c3}")
    print("Hamiltonian path:")
    print("  " + " > ".join(sorted(names, key=lambda n: -score[n])))
    print()
    print("Assumption challenge:")
    print("  considered vertices: runners, residues, divisors, unit twists, gcd strata,")
    print("    primitive Fourier modes, endpoint owners, C27/K33 packets, and proof obligations.")
    print("  chosen vertices here: quotient channels.")
    print("  preserved predicate: whether a channel retains enough payload to support the next")
    print("    LRC implication without silently mixing AP/GW, q-witness, K33, petal, and covering routes.")
    print("  destroyed information: raw runner identity when it is not part of the next predicate.")


def proof_target() -> None:
    section("HYP-2978 PROOF TARGET")
    print("Ramanujan-divisor admissibility theorem target:")
    print()
    print("  A quotient Q used in an LRC14 proof step is admissible only if every")
    print("  coordinate it forgets is one of:")
    print("    (1) invariant under the step's allowed moves;")
    print("    (2) reconstructible from retained q/Farey, unit, gcd, and Ramanujan data;")
    print("    (3) annihilated by a proved orthogonality or dual certificate; or")
    print("    (4) placed in an explicit residual bucket with endpoint/state-lift labels.")
    print()
    print("Concrete LRC14 route:")
    print("  Use divisor/totient/Jordan/psi data as capacity ledgers.")
    print("  Use Ramanujan c_q profiles as primitive-shell Fourier guards.")
    print("  Then hand off to HYP-2974/2977 harmonic duals or HYP-2965/2969 endpoint")
    print("  packets before claiming that a qdiv>14 zero-open packet cannot exist.")
    print()
    print("Most important new falsifier:")
    print("  two rows sharing a scalar divisor/Ramanujan signature but taking different")
    print("  LRC proof routes.  The audit already exhibits this for qdiv/open-state/")
    print("  residue-shell quotients, so those channels must remain labelled features,")
    print("  not final proof quotients.")


def main() -> None:
    verify_arithmetic_identities(80)
    audits = [audit_row(row) for row in named_rows()]
    print_row_audit(audits)
    signature_collisions(audits)
    tournament_analysis()
    proof_target()


if __name__ == "__main__":
    main()
