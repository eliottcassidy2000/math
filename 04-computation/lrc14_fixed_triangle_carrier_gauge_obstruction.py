#!/usr/bin/env python3
"""Exact THM-2763 quotient-gauge audit for the proved THM-2749 carrier."""

from hashlib import sha256
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
HELPER = HERE / "lrc14_extended_carrier_endpoint_lib.py"
HELPER_SHA256 = "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57"
if sha256(HELPER.read_bytes().replace(b"\r\n", b"\n")).hexdigest() != HELPER_SHA256:
    raise RuntimeError("extended-carrier helper changed")
import lrc14_extended_carrier_endpoint_lib as base


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def q_value(present, carrier):
    embedding = (base.endpoint.MODS[0],)
    restricted = base.intersect_sorted(present, carrier)
    return base.endpoint.endpoint_sum(restricted, -base.Y0, embedding)[0]


def pq_value(present, carrier, terminal, q_starts, tabs):
    prime = base.endpoint.MODS[0][0]
    embedding = (base.endpoint.MODS[0],)
    restricted = base.intersect_sorted(present, carrier)
    left, _overlap = base.endpoint.x_sweep(
        restricted, terminal, q_starts, base.X0, embedding, tabs,
    )
    right = base.endpoint.endpoint_sum(
        restricted, -base.Y0, embedding,
    )
    return (
        left[0],
        right[0],
        left[0] * right[0] % prime,
        base.interval_mass(restricted),
    )


def main():
    source, target, terminal = base.build_carriers()["base"]
    carriers = {"source": source, "target": target}
    prime_root = (base.endpoint.MODS[0],)
    q_starts = [left for left, _right in terminal]
    tabs = base.endpoint.make_tabs(terminal, base.X0, prime_root)
    addresses = ((0, 0), (1, 0), (0, 1), (3, 7))
    gauge_shift = -(base.T // base.P)
    rows = []
    bare_controls = []
    for address in addresses:
        ell = base.REPS[address]
        shifted_ell = tuple(
            (ell[index] + base.WMOD[index]) % base.P
            for index in range(9)
        )
        present = tuple(base.endpoint.build_set(base.PAT_E3, ell))
        shifted_present = tuple(
            base.endpoint.build_set(base.PAT_E3, shifted_ell)
        )
        require(
            shifted_present
            == base.marked.clutch.shift_union(present, gauge_shift),
            "ell+W physical shift law",
        )

        bare = base.endpoint.endpoint_sum(
            present, -base.Y0, (base.endpoint.MODS[0],)
        )[0]
        bare_shifted = base.endpoint.endpoint_sum(
            shifted_present, -base.Y0, (base.endpoint.MODS[0],)
        )[0]
        require(bare == bare_shifted, "bare quotient-gauge descent")
        bare_controls.append((address, bare))

        for name, carrier in carriers.items():
            transported_carrier = base.marked.clutch.shift_union(
                carrier, gauge_shift
            )
            original = q_value(present, carrier)
            fixed = q_value(shifted_present, carrier)
            transported = q_value(shifted_present, transported_carrier)
            require(original == transported,
                    f"transported carrier gauge covariance: {name}")
            rows.append((address, name, original, fixed,
                         original == fixed))

    failures = tuple(row for row in rows if not row[-1])
    expected_bare_controls = (
        ((0, 0), 272457584061297438),
        ((1, 0), 12491018700458795),
        ((0, 1), 244255279417483343),
        ((3, 7), 146843097480870949),
    )
    expected_rows = (
        ((0, 0), "source", 287659712270709994, 0, False),
        ((0, 0), "target", 248235870634784933, 0, False),
        ((1, 0), "source", 287659712270709994, 0, False),
        ((1, 0), "target", 248235870634784933, 0, False),
        ((0, 1), "source", 69469253769051213, 69469253769051213, True),
        ((0, 1), "target", 248235870634784933, 248235870634784933, True),
        ((3, 7), "source", 153998791602873972, 153998791602873972, True),
        ((3, 7), "target", 289285119589719712, 289285119589719712, True),
    )
    require(tuple(bare_controls) == expected_bare_controls,
            "bare control values changed")
    require(tuple(rows) == expected_rows,
            "fixed-carrier gauge rows changed")
    require(len(failures) == 4, "fixed-carrier failure count changed")

    ell0 = base.REPS[(0, 0)]
    ell0_plus_w = tuple(
        (ell0[index] + base.WMOD[index]) % base.P
        for index in range(9)
    )
    present0 = tuple(base.endpoint.build_set(base.PAT_E3, ell0))
    present0_plus_w = tuple(
        base.endpoint.build_set(base.PAT_E3, ell0_plus_w)
    )
    full_rows = []
    for name, carrier in carriers.items():
        transported_carrier = base.marked.clutch.shift_union(
            carrier, gauge_shift
        )
        original = pq_value(
            present0, carrier, terminal, q_starts, tabs,
        )
        fixed = pq_value(
            present0_plus_w, carrier, terminal, q_starts, tabs,
        )
        transported = pq_value(
            present0_plus_w, transported_carrier, terminal, q_starts, tabs,
        )
        require(original == transported,
                f"full transported-current gauge covariance: {name}")
        full_rows.append((name, original, fixed))
    expected_full_rows = (
        (
            "source",
            (189041250036777056, 287659712270709994,
             65280867241115379, 6320326320),
            (0, 0, 0, 0),
        ),
        (
            "target",
            (218344733173586894, 248235870634784933,
             209108233808250489, 6320326320),
            (0, 0, 0, 0),
        ),
    )
    require(tuple(full_rows) == expected_full_rows,
            "full P/Q/H gauge hostile changed")
    print("THM-2763 fixed-triangle quotient-gauge audit")
    print(f"helper_sha256={HELPER_SHA256}")
    print(f"gauge=ell->ell+W; support_shift={gauge_shift}/T=-1/13")
    print(f"bare_Q_descent_controls={tuple(bare_controls)}")
    print(f"fixed_carrier_rows={tuple(rows)}")
    print(f"fixed_carrier_descent_failures={len(failures)}/{len(rows)}")
    print(f"full_PQH_mass_hostile_at_(0,0)={tuple(full_rows)}")
    print("transported_carrier_descent=PASS")
    print("CONCLUSION: the fixed collapsed common carrier breaks the <W> representative descent. Its 169 canonical-representative values are not yet a lawful THM-2334 quotient bank. Transporting the carrier with the gauge repairs descent, but exact-address use still requires the carrier harmonic/factor sidecar.")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
