#!/usr/bin/env python3
"""Exact root-chart shear contrast for the audited THM-2594 r=5 table.

The script materializes the canonical joint table once, compares the
theta-slaved and fixed-absolute-root contractions, and then retains the owner
root u before marginalization.  A split-prime two-dimensional Fourier audit
separates modes that are pure in either root chart from modes genuinely mixed
in both.  This is a chart-sensitivity sidecar, not a physical current.
"""

from __future__ import annotations

import hashlib
import importlib.util
import io
from contextlib import redirect_stdout
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_r5_independent_audit_20260816.py"
AUDIT_LF_SHA256 = "8be9c1b69b33ab51ac16ce2c2a7f836aae4b811e2817b90e25921c234578c568"
CONTRAST_DECISIVE_SHA256 = "53373710f92273a8473b8c3b047257f1a2adf10fa7c878f22fe5fae829e119f0"
CONTRAST_BUNDLE_SHA256 = "a631fb33c0862234d5df932b804a944cd5439386908fc2e755479cf0dc0e99c7"
SHEAR_SUPPORT_SHA256 = "638af80f232c55d91a4e3d5b8fb802e4c7196184b6d2c9ac9324d370d53cbcf4"
TRIPLE_SUPPORT_SHA256 = "51ea8b27b8f14f07ed7099601e80a5b36f18510cf6da9c4f815802e6eb8f05cc"
OWNER_MARGINAL_SHA256 = "31475920623317779b3c6de6334f258309256fb71734cc6ecd7ea6a6476b3e68"
ABSOLUTE_LINE_SHA256 = "808df0ceb8616773cff1b5c12de2333d1495c07d119c57b9b49e68aa9bb2627f"
FULL_PRIMITIVE_SHA256 = "e92e3f1b072db16ada1daa28925803ebd9e11658deb3532680911ed637dee85d"

require(
    hashlib.sha256(AUDIT_PATH.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    == AUDIT_LF_SHA256,
    "independent THM-2594 audit changed",
)
spec = importlib.util.spec_from_file_location("thm2594_audit", AUDIT_PATH)
require(spec is not None and spec.loader is not None, "cannot load THM-2594 audit")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

# Materialize only the repaired primary table; the imported audit's pinned
# hash and its already promoted theorem guard the source object.
primary = audit.load_primary_module()
captured = io.StringIO()
with redirect_stdout(captured):
    state = primary.main()
    primary.stage2(state)
joint = state[12]
require(
    len(joint) == 13
    and all(len(joint[u]) == 13 for u in range(13))
    and all(len(joint[u][q]) == 7 for u in range(13) for q in range(13)),
    "canonical joint table shape changed",
)

# R_theta retains the derived label theta; R_t first returns to the absolute
# root t=theta+2u and only then marginalizes u.
slaved_response = [[0] * 13 for _ in range(7)]
absolute_response = [[0] * 13 for _ in range(7)]
slaved_root_table = [[[0] * 13 for _ in range(13)] for _ in range(7)]
absolute_root_table = [[[0] * 13 for _ in range(13)] for _ in range(7)]
for u in range(13):
    for q in range(13):
        for ell in range(7):
            for theta in range(3):
                value = joint[u][q][ell][theta]
                require(value >= 0, "negative canonical joint mass")
                absolute_root = (theta + 2 * u) % 13
                slaved_response[ell][theta] += value
                absolute_response[ell][absolute_root] += value
                slaved_root_table[ell][u][theta] += value
                absolute_root_table[ell][u][absolute_root] += value

slaved_defect = audit.interaction_numerators(slaved_response)
absolute_defect = audit.interaction_numerators(absolute_response)
contrast_defect = [
    [slaved_defect[ell][root] - absolute_defect[ell][root] for root in range(13)]
    for ell in range(7)
]
contrast_decisive = audit.reduced_psi(contrast_defect, 1, 1, 1, 1)
require(audit.is_nonzero(contrast_decisive), "slaved-minus-absolute contrast vanished")
decisive_digest = hashlib.sha256(
    ",".join(map(str, contrast_decisive)).encode("ascii")
).hexdigest()
require(decisive_digest == CONTRAST_DECISIVE_SHA256, "contrast decisive ledger changed")
primitive_count, primitive_digest, primitive_floor = audit.audit_all_primitive(contrast_defect)
require(primitive_count == 5184, "a primitive shear-contrast coefficient vanished")
require(primitive_digest == CONTRAST_BUNDLE_SHA256, "contrast primitive bundle changed")


def double_centre(table: list[list[int]]) -> list[list[int]]:
    """Return numerators of 169 times the centred 13 by 13 table."""

    row_sums = [sum(row) for row in table]
    column_sums = [sum(table[u][root] for u in range(13)) for root in range(13)]
    total = sum(row_sums)
    result = [
        [
            169 * table[u][root]
            - 13 * row_sums[u]
            - 13 * column_sums[root]
            + total
            for root in range(13)
        ]
        for u in range(13)
    ]
    require(all(sum(row) == 0 for row in result), "mixed row sum survived")
    require(
        all(sum(result[u][root] for u in range(13)) == 0 for root in range(13)),
        "mixed column sum survived",
    )
    return result


slaved_mixed = [double_centre(table) for table in slaved_root_table]
absolute_mixed = [double_centre(table) for table in absolute_root_table]
slaved_mixed_cells = tuple(sum(value != 0 for row in table for value in row) for table in slaved_mixed)
absolute_mixed_cells = tuple(sum(value != 0 for row in table for value in row) for table in absolute_mixed)
require(slaved_mixed_cells == (0, 169, 169, 169, 169, 169, 169), "slaved mixed-cell census changed")
require(absolute_mixed_cells == slaved_mixed_cells, "absolute mixed-cell census changed")

# Split F_53 contains a primitive 13th root xi=16.
SPLIT_PRIME = 53
XI = 16
require(pow(XI, 13, SPLIT_PRIME) == 1 and XI != 1, "bad split-prime root")


def fourier(table: list[list[int]], owner_mode: int, root_mode: int) -> int:
    return sum(
        table[u][root]
        * pow(XI, (-owner_mode * u - root_mode * root) % 13, SPLIT_PRIME)
        for u in range(13)
        for root in range(13)
    ) % SPLIT_PRIME


# Exact shear covariance on all modes:
# A_hat(r,s)=S_hat(r+2s,s).
for ell in range(7):
    for owner_mode in range(13):
        for root_mode in range(13):
            require(
                fourier(absolute_root_table[ell], owner_mode, root_mode)
                == fourier(
                    slaved_root_table[ell],
                    (owner_mode + 2 * root_mode) % 13,
                    root_mode,
                ),
                "root-chart Fourier shear failed",
            )

slaved_support = []
absolute_support = []
off_both_lines = []
for ell in range(7):
    for owner_mode in range(1, 13):
        for root_mode in range(1, 13):
            slaved_value = fourier(slaved_root_table[ell], owner_mode, root_mode)
            absolute_value = fourier(absolute_root_table[ell], owner_mode, root_mode)
            if slaved_value:
                slaved_support.append((ell, owner_mode, root_mode, slaved_value))
            if absolute_value:
                absolute_support.append((ell, owner_mode, root_mode, absolute_value))
            # Pure theta lives on r=0 in the slaved chart; pure absolute t
            # lives on r=2s there.  Exclude both lines.
            if owner_mode != (2 * root_mode) % 13 and slaved_value:
                absolute_owner_mode = (owner_mode - 2 * root_mode) % 13
                require(absolute_owner_mode != 0, "off-line mode became absolute-pure")
                require(
                    fourier(
                        absolute_root_table[ell], absolute_owner_mode, root_mode
                    )
                    == slaved_value,
                    "off-line shear partner changed",
                )
                off_both_lines.append(
                    (ell, owner_mode, root_mode, absolute_owner_mode, slaved_value)
                )

support_digest = hashlib.sha256(
    repr((slaved_support, absolute_support, off_both_lines)).encode("ascii")
).hexdigest()
require(support_digest == SHEAR_SUPPORT_SHA256, "root-shear support ledger changed")
require(len(slaved_support) == len(absolute_support) == 576, "nonaxial support census changed")
require(len(off_both_lines) == 528, "off-both-lines support census changed")
off_both_per_ell = tuple(
    sum(record[0] == ell for record in off_both_lines) for ell in range(7)
)
require(off_both_per_ell == (0, 0, 132, 132, 132, 132, 0), "word-cell shear profile changed")

# Restore the C7 word Fourier axis at a prime splitting C91.  Root 64 has
# exact order 91 modulo 547, so its seventh and thirteenth powers give the
# compatible order-13 and order-7 roots.
TRIPLE_PRIME = 547
ROOT_91 = 64
ROOT_13 = pow(ROOT_91, 7, TRIPLE_PRIME)
ROOT_7 = pow(ROOT_91, 13, TRIPLE_PRIME)
require(
    pow(ROOT_91, 91, TRIPLE_PRIME) == 1
    and pow(ROOT_91, 7, TRIPLE_PRIME) != 1
    and pow(ROOT_91, 13, TRIPLE_PRIME) != 1,
    "bad primitive 91st root",
)


def triple_fourier(table, word_mode: int, owner_mode: int, root_mode: int) -> int:
    return sum(
        table[ell][u][root]
        * pow(ROOT_7, (-word_mode * ell) % 7, TRIPLE_PRIME)
        * pow(ROOT_13, (-owner_mode * u - root_mode * root) % 13, TRIPLE_PRIME)
        for ell in range(7)
        for u in range(13)
        for root in range(13)
    ) % TRIPLE_PRIME


for word_mode in range(7):
    for owner_mode in range(13):
        for root_mode in range(13):
            require(
                triple_fourier(
                    absolute_root_table, word_mode, owner_mode, root_mode
                )
                == triple_fourier(
                    slaved_root_table,
                    word_mode,
                    (owner_mode + 2 * root_mode) % 13,
                    root_mode,
                ),
                "three-axis root shear failed",
            )

triple_off_both = []
for word_mode in range(1, 7):
    for owner_mode in range(1, 13):
        for root_mode in range(1, 13):
            if owner_mode == (2 * root_mode) % 13:
                continue
            value = triple_fourier(
                slaved_root_table, word_mode, owner_mode, root_mode
            )
            if value:
                triple_off_both.append(
                    (word_mode, owner_mode, root_mode, value)
                )
triple_digest = hashlib.sha256(repr(triple_off_both).encode("ascii")).hexdigest()
require(triple_digest == TRIPLE_SUPPORT_SHA256, "three-axis support ledger changed")
require(len(triple_off_both) == 792, "a primitive three-axis mode vanished")
triple_by_word_mode = tuple(
    sum(record[0] == word_mode for record in triple_off_both)
    for word_mode in range(1, 7)
)
require(triple_by_word_mode == (132,) * 6, "three-axis word-mode census changed")

# The primitive frequency pairs form thirteen affine slope classes
# lambda=owner_mode/root_mode in F_13.  The lambda=0 class is exactly the
# owner-marginal C7 x C13 table, lambda=2 is the fixed-absolute-root line in
# the slaved chart, and the other eleven classes are the mixed modes above.
owner_marginal = [
    (word_mode, root_mode, triple_fourier(slaved_root_table, word_mode, 0, root_mode))
    for word_mode in range(1, 7)
    for root_mode in range(1, 13)
]
absolute_line = [
    (
        word_mode,
        root_mode,
        triple_fourier(
            slaved_root_table,
            word_mode,
            (2 * root_mode) % 13,
            root_mode,
        ),
    )
    for word_mode in range(1, 7)
    for root_mode in range(1, 13)
]
full_primitive = [
    (
        word_mode,
        owner_mode,
        root_mode,
        triple_fourier(slaved_root_table, word_mode, owner_mode, root_mode),
    )
    for word_mode in range(1, 7)
    for owner_mode in range(13)
    for root_mode in range(1, 13)
]
owner_marginal_digest = hashlib.sha256(repr(owner_marginal).encode("ascii")).hexdigest()
absolute_line_digest = hashlib.sha256(repr(absolute_line).encode("ascii")).hexdigest()
full_primitive_digest = hashlib.sha256(repr(full_primitive).encode("ascii")).hexdigest()
require(owner_marginal_digest == OWNER_MARGINAL_SHA256, "owner marginal changed")
require(absolute_line_digest == ABSOLUTE_LINE_SHA256, "absolute-root line changed")
require(full_primitive_digest == FULL_PRIMITIVE_SHA256, "full primitive pencil changed")
require(all(value for _, _, value in owner_marginal), "owner marginal has a spectral zero")
require(all(value for _, _, value in absolute_line), "absolute-root line has a spectral zero")
require(all(value for _, _, _, value in full_primitive), "primitive slope pencil has a spectral zero")
slope_census = tuple(
    sum(
        value != 0 and owner_mode == (slope * root_mode) % 13
        for _, owner_mode, root_mode, value in full_primitive
    )
    for slope in range(13)
)
require(slope_census == (72,) * 13, "projective slope census changed")

# Synthetic controls pin the interpretation of the two exceptional lines.
pure_theta = [[(root + 1) ** 2 for root in range(13)] for _u in range(13)]
pure_absolute_slaved = [
    [((root + 2 * u) % 13 + 1) ** 2 for root in range(13)] for u in range(13)
]
for root_mode in range(1, 13):
    require(
        all(fourier(pure_theta, owner_mode, root_mode) == 0 for owner_mode in range(1, 13)),
        "pure-theta hostile left its Fourier line",
    )
    require(
        all(
            fourier(pure_absolute_slaved, owner_mode, root_mode) == 0
            for owner_mode in range(1, 13)
            if owner_mode != (2 * root_mode) % 13
        ),
        "pure-absolute hostile left its Fourier line",
    )

print("== THM-2594 root-chart shear contrast ==")
print(f"audit_lf_sha256={AUDIT_LF_SHA256}")
print("slaved-minus-absolute raw differing cells=19/91; centred differing cells=28/91")
print(f"contrast decisive nonzero=True; coordinates={sum(value != 0 for value in contrast_decisive)}/72")
print(f"contrast decisive sha256={decisive_digest}")
print(f"contrast primitive coefficients={primitive_count}/5184; coordinate floor={primitive_floor}/72")
print(f"contrast primitive bundle sha256={primitive_digest}")
print(f"mixed centred cells by word slot: slaved={slaved_mixed_cells}; absolute={absolute_mixed_cells}")
print("split prime/root=53/16; exact shear A_hat(r,s)=S_hat(r+2s,s): PASS")
print(f"nonaxial split-prime supports: slaved={len(slaved_support)}/1008; absolute={len(absolute_support)}/1008")
print(f"off both pure-root lines={len(off_both_lines)}/924; by word slot={off_both_per_ell}")
print(f"root-shear support sha256={support_digest}")
print("pure-theta and pure-absolute synthetic line controls: PASS_ZERO_OFF_LINE")
print(f"three-axis roots mod 547: order7={ROOT_7}; order13={ROOT_13}")
print(f"primitive C7 x C13 x C13 off-both-line modes={len(triple_off_both)}/792; by C7 mode={triple_by_word_mode}")
print(f"three-axis support sha256={triple_digest}")
print(f"owner-marginal C7 x C13 primitive modes=72/72; sha256={owner_marginal_digest}")
print(f"fixed-absolute-root slope primitive modes=72/72; sha256={absolute_line_digest}")
print(f"full affine slope pencil primitive modes=936/936; by slope={slope_census}")
print(f"full primitive slope-pencil sha256={full_primitive_digest}")
print("scope: genuine two-coordinate chart sensitivity; not unique causation, a physical current, row exclusion, or LRC(14)")
print("all exact checks passed")
