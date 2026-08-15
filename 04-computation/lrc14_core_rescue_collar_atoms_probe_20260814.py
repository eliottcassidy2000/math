#!/usr/bin/env python3
"""Exact collar anatomy of all twelve q<=7 core rescues.

The THM-3387 rescue ledger is reconstructed from four transverse cover atoms.
Their full-cover phase sets are translated pairs of rational collars.  The
q=3 and q=6 atoms share the same radial collar; the q=6 atom fills the missing
parity centres.  Runtime gates survive python -O.

This is an unnumbered FINITE-EXACT sidecar, not a theorem promotion.
"""

from __future__ import annotations

import ast
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "THM-3387-script",
        ROOT / "04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py",
        "9b0b46874a569d674b937b37cf74a8985fca2b77e3e480a75fb4924ea602f25a",
    ),
    (
        "THM-3387-output",
        ROOT / "05-knowledge/results/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.out",
        "b4d9ce439bab4501bfd5e2cf13eb0b0e3685b7364f30e43b7d5ca9138d25cb5c",
    ),
    (
        "q3-script",
        ROOT / "04-computation/lrc14_q3_phase_triangle_clutter_thm3388.py",
        "5323346310a9a6b188caa0131b177b2ae8e23c7113808cda8955f89828e62154",
    ),
    (
        "q3-output",
        ROOT / "05-knowledge/results/lrc14_q3_phase_triangle_clutter_thm3388.out",
        "5a32319fb8a91b476d292da292ae3cc9933f5f94aad7eb0e834f49e52252c535",
    ),
    (
        "q5-script",
        ROOT / "04-computation/lrc14_q5_singleton_cover_clutter_probe_20260814.py",
        "96b5addc012546b051987d623eff01581f6616219228524f9465be63f09776a0",
    ),
    (
        "q5-output",
        ROOT / "05-knowledge/results/lrc14_q5_singleton_cover_clutter_probe_20260814.out",
        "83a76d6d04be5a7721d070f923612c9f034ed3dac72b5f84636d8632796facfe",
    ),
    (
        "q6-script",
        ROOT / "04-computation/lrc14_q6_typed_cover_clutter_probe_20260814.py",
        "e16ae18816b0620dde6ae05db01846ce49a37c292c7cccc76532ec2fc9c728d6",
    ),
    (
        "q6-output",
        ROOT / "05-knowledge/results/lrc14_q6_typed_cover_clutter_probe_20260814.out",
        "3628260d408258b1f9ad385c23c4d4bdea8cf0c7bad65f02dbbff9dc0e5de889",
    ),
)

MODULE_PATHS = (
    ("thm3387", ROOT / "04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py"),
    ("q3", ROOT / "04-computation/lrc14_q3_phase_triangle_clutter_thm3388.py"),
    ("q5", ROOT / "04-computation/lrc14_q5_singleton_cover_clutter_probe_20260814.py"),
    ("q6", ROOT / "04-computation/lrc14_q6_typed_cover_clutter_probe_20260814.py"),
)

Q3_EDGE = (8, 11, 13)
Q5_EDGE_A = (6, 11, 12, 13, 14)
Q5_EDGE_B = (7, 9, 11, 12, 13)
Q6_EDGE = (3, 8, 11, 13)
Q6_INERT_EXTENSIONS = (1, 5, 7, 9, 10, 14)
Q6_ACTIVE_EXTENSIONS = (2, 4)
EXPECTED_BOUNDARY_ANATOMY = (
    ("q3", ((11, 0),), ((8, 1),), ((8, (1,)), (11, (0,)), (13, (2,)))),
    (
        "q5a",
        ((11, 4),),
        ((6, 0),),
        ((6, (0,)), (11, (4,)), (12, (2,)), (13, (3,)), (14, (1,))),
    ),
    (
        "q5b",
        ((7, 3),),
        ((12, 2),),
        ((7, (3,)), (9, (0,)), (11, (4,)), (12, (2,)), (13, (1,))),
    ),
    (
        "q6",
        ((11, 1),),
        ((8, 0), (8, 3)),
        ((3, (0, 2, 4)), (8, (0, 3)), (11, (1,)), (13, (5,))),
    ),
)
EXPECTED_SEMANTIC_DIGEST = "3a1d4bba94086d82c12c94225fc0403be165cc6e7bb5245fde07d4bc60d7c022"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module spec", name))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


def strict_danger(numerator, denominator):
    residue = numerator % denominator
    return 14 * min(residue, denominator - residue) < denominator


def circle_norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def boundary_anatomy(name, q, edge, centre, inner, outer):
    critical = []
    for radius in (inner, outer):
        point = centre + radius
        critical.append(
            tuple(
                (speed, sheet)
                for speed in edge
                for sheet in range(q)
                if circle_norm(speed * (point + Q(sheet, q))) == Q(1, 14)
            )
        )
    midpoint = centre + (inner + outer) / 2
    blocks = tuple(
        (
            speed,
            tuple(
                sheet
                for sheet in range(q)
                if circle_norm(speed * (midpoint + Q(sheet, q))) < Q(1, 14)
            ),
        )
        for speed in edge
    )
    return name, critical[0], critical[1], blocks


def transverse_danger(q, speed, sample, scale, sheet):
    return strict_danger(
        speed * (q * sample + 2 * scale * sheet),
        2 * q * scale,
    )


def core_danger(q, clock, sample, scale):
    return strict_danger(q * clock * sample, 2 * scale)


def event_points(q, transverse, core=()):
    scale = 14 * q * lcm(*transverse, *core)
    events = {0}
    for speed in transverse:
        for sheet in range(q):
            for tooth in range(speed):
                for sign in (-1, 1):
                    events.add(
                        (
                            scale * tooth // speed
                            - scale * sheet // q
                            + sign * scale // (14 * speed)
                        )
                        % scale
                    )
    for clock in core:
        for tooth in range(q * clock):
            for sign in (-1, 1):
                events.add(
                    (
                        scale * tooth // (q * clock)
                        + sign * scale // (14 * q * clock)
                    )
                    % scale
                )
    return scale, tuple(sorted(events))


def full_cover(q, transverse, sample, scale):
    return all(
        any(
            transverse_danger(q, speed, sample, scale, sheet)
            for speed in transverse
        )
        for sheet in range(q)
    )


def cover_components(q, transverse):
    scale, events = event_points(q, transverse)
    count = len(events)
    live = []
    for index, left in enumerate(events):
        right = events[(index + 1) % count]
        if index + 1 == count:
            right += scale
        live.append(full_cover(q, transverse, (left + right) % (2 * scale), scale))
    require(not all(live), ("whole circle covered", q, transverse))

    components = []
    for index in range(count):
        if not live[index] or live[(index - 1) % count]:
            continue
        stop = index
        while live[stop % count]:
            stop += 1
        left = events[index]
        right = events[stop % count] + (scale if stop >= count else 0)
        components.append((Q(left, scale), Q(right, scale)))
    return tuple(components)


def leak_samples(q, transverse, core):
    scale, events = event_points(q, transverse, core)
    count = len(events)
    samples = tuple(2 * event for event in events) + tuple(
        (
            events[index]
            + events[(index + 1) % count]
            + (scale if index + 1 == count else 0)
        )
        % (2 * scale)
        for index in range(count)
    )
    return tuple(
        Q(sample, 2 * scale)
        for sample in samples
        if full_cover(q, transverse, sample, scale)
        and not any(core_danger(q, clock, sample, scale) for clock in core)
    )


def collars(centres, inner, outer):
    intervals = []
    for centre in centres:
        for left, right in (
            (centre - outer, centre - inner),
            (centre + inner, centre + outer),
        ):
            while left < 0:
                left += 1
                right += 1
            while left >= 1:
                left -= 1
                right -= 1
            if right <= 1:
                intervals.append((left, right))
            else:
                intervals.append((left, Q(1)))
                intervals.append((Q(0), right - 1))
    return tuple(sorted(intervals))


def minimal_cover_atom(q, edge):
    return bool(cover_components(q, edge)) and all(
        not cover_components(q, edge[:index] + edge[index + 1 :])
        for index in range(len(edge))
    )


def main():
    for name, path, expected in DEPENDENCIES:
        require(lf_hash(path) == expected, ("dependency changed", name, path))
    modules = {name: load_module(name, path) for name, path in MODULE_PATHS}

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    atom_specs = (
        (
            "q3",
            3,
            Q3_EDGE,
            tuple(Q(2 * index + 1, 6) for index in range(3)),
            Q(2, 231),
            Q(1, 112),
            (2,),
        ),
        (
            "q5a",
            5,
            Q5_EDGE_A,
            tuple(Q(index, 5) for index in range(5)),
            Q(9, 770),
            Q(1, 84),
            (1,),
        ),
        (
            "q5b",
            5,
            Q5_EDGE_B,
            tuple(Q(2 * index + 1, 10) for index in range(5)),
            Q(1, 245),
            Q(1, 168),
            (2,),
        ),
        (
            "q6",
            6,
            Q6_EDGE,
            tuple(Q(index, 6) for index in range(6)),
            Q(2, 231),
            Q(1, 112),
            (1,),
        ),
    )

    atoms = []
    anatomies = []
    for name, q, edge, centres, inner, outer, core in atom_specs:
        actual = cover_components(q, edge)
        expected = collars(centres, inner, outer)
        require(actual == expected, ("collar formula", name, actual, expected))
        require(minimal_cover_atom(q, edge), ("not minimal", name, edge))
        require(not leak_samples(q, edge, core), ("core leak", name, edge, core))
        require(
            outer < Q(1, 14 * q * core[0]),
            ("collar not inside core radius", name),
        )
        anatomies.append(boundary_anatomy(name, q, edge, centres[0], inner, outer))
        atoms.append((name, q, edge, centres, inner, outer, actual, core))
    require(tuple(anatomies) == EXPECTED_BOUNDARY_ANATOMY, anatomies)

    q3_components = atoms[0][6]
    q6_components = atoms[3][6]
    require(
        q6_components
        == tuple(
            sorted(
                q3_components
                + tuple(
                    (
                        (left - Q(1, 6)) % 1,
                        (right - Q(1, 6)) % 1,
                    )
                    for left, right in q3_components
                )
            )
        ),
        "q3/q6 parity fill",
    )

    inert = []
    active_controls = []
    for extension in Q6_INERT_EXTENSIONS + Q6_ACTIVE_EXTENSIONS:
        transverse = tuple(sorted(Q6_EDGE + (extension,)))
        components = cover_components(6, transverse)
        leaks = leak_samples(6, transverse, (1,))
        if components == q6_components and not leaks:
            inert.append(extension)
        else:
            active_controls.append((extension, Q(1, 84) in leaks, leaks[:4]))
    require(tuple(inert) == Q6_INERT_EXTENSIONS, inert)
    require(
        tuple(item[0] for item in active_controls) == Q6_ACTIVE_EXTENSIONS
        and all(item[1] for item in active_controls),
        active_controls,
    )

    generated_records = []
    for core in combinations((1, 2, 3, 4), 3):
        if 2 in core:
            generated_records.append((3, core, Q3_EDGE))
    generated_records.extend(
        (
            (5, (1,), Q5_EDGE_A),
            (5, (2,), Q5_EDGE_B),
            (6, (1, 2), Q6_EDGE),
        )
    )
    generated_records.extend(
        (6, (1,), tuple(sorted(Q6_EDGE + (extension,))))
        for extension in Q6_INERT_EXTENSIONS
    )
    generated_records = tuple(sorted(generated_records))

    expected_records = tuple(
        sorted(
            (q, tuple(core), tuple(transverse))
            for _, _, _, q, core, transverse in modules["thm3387"].EXPECTED_CORE_RESCUES
        )
    )
    require(generated_records == expected_records, (generated_records, expected_records))
    require(
        set(q for q, _, _ in generated_records) == {3, 5, 6}
        and len(generated_records) == 12,
        "rescue support",
    )
    require(
        tuple(modules["q3"].EXPECTED_CORE_RESCUES)
        == tuple((core, transverse) for q, core, transverse in generated_records if q == 3),
        "q3 rescue module",
    )
    require(
        tuple(modules["q5"].EXPECTED_CORE_RESCUES)
        == tuple((core, transverse) for q, core, transverse in generated_records if q == 5),
        "q5 rescue module",
    )
    require(
        tuple(modules["q6"].EXPECTED_CORE_RESCUES)
        == tuple((core, transverse) for q, core, transverse in generated_records if q == 6),
        "q6 rescue module",
    )

    semantic = ExactDigest()
    semantic.add(("atoms", tuple(atoms)))
    semantic.add(("boundary_anatomy", tuple(anatomies)))
    semantic.add(("q6_extensions", tuple(inert), tuple(active_controls)))
    semantic.add(("records", generated_records))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("LRC14 CORE-RESCUE COLLAR ATOMS EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in DEPENDENCIES)}")
    print("status=FINITE-EXACT unnumbered sidecar;four collar atoms generate all twelve q<=7 core rescues")
    print("q3_atom=edge:(8,11,13);centres:(2j+1)/6,j=0..2;radial_collar:(2/231,1/112);core_clock=2")
    print("q5_atom_a=edge:(6,11,12,13,14);centres:j/5,j=0..4;radial_collar:(9/770,1/84);core_clock=1")
    print("q5_atom_b=edge:(7,9,11,12,13);centres:(2j+1)/10,j=0..4;radial_collar:(1/245,1/168);core_clock=2")
    print("q6_atom=edge:(3,8,11,13);centres:j/6,j=0..5;radial_collar:(2/231,1/112);core_clock=1")
    print("branch_transplant=q3_odd_sixth_collars_plus_even_sixth_fill_equals_q6_all_sixth_collars")
    print(f"boundary_anatomy={tuple(anatomies)}")
    print(f"q6_inert_extensions={tuple(inert)};active_extensions={tuple(active_controls)}")
    print(f"rescue_records={generated_records}")
    print("rescue_count=3(q3)+2(q5)+7(q6)=12;q2_and_q4_have_none")
    print("scope=exact_phase_cell_anatomy;no_new_row,no_refined_decrement,no_LRC14")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
