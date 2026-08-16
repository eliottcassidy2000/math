#!/usr/bin/env python3
"""Pull the actual U_full endpoint pair function through source-time labels.

The source dependency realizes THM-2471 (34): both 39-atom labels are read at
endpoint/source time, carried through the K=2 inverse-branch sidecar, and met
on one first-collision base before marginalization.  The target dependency
independently reconstructs THM-3514's actual endpoint pair function E.

Finite linearity then realizes

    sum_(omega,nu,c) M(omega,nu,c) E(omega,nu) zeta^(-kc)

as one simple-kernel integral on that ancestry base.  This closes the
arrival-versus-source atom-label mismatch for the finite address function.
It does not identify the AX/BY endpoint integrations with the ancestry
integrand, construct a THM-2334 relation current, or exclude an LRC row.
"""

from __future__ import annotations

from contextlib import redirect_stdout
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib.util
import io
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
P = 13
SOURCE_PATH = ROOT / "04-computation/lrc_r5_source_aligned_guard_atom_branch_sidecar_probe_20260816.py"
SOURCE_SHA256 = "22c5c748392817ccc36889a007c65bd5f44b26c10638df6f6aac48e917547f41"
SOURCE_SEMANTIC = "31e9e90c63053944b590195555be07ccbf84fd4c7abc2101de6a2a3562202de6"
TARGET_PATH = ROOT / "04-computation/lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_independent_audit_20260816.py"
TARGET_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
TARGET_BUCKET_SHA256 = "553cfe7289b0556a19a8bcd1a0382dc1372545358feac62e5229adca315f8a26"
EXPECTED_SEMANTIC_SHA256 = "430aa2201ee8f624de66f034bdf5e31b11ec6ecbe1f3f0b9cef1342b2b1bc2ff"
_TARGET_WORKER_MODULE = None


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    body = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return sha256(body).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    require(lf_sha256(path) == expected_hash, (name, "source hash drift"))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "module loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def target_worker(alpha: int):
    global _TARGET_WORKER_MODULE
    if _TARGET_WORKER_MODULE is None:
        _TARGET_WORKER_MODULE = load_module(
            TARGET_PATH, "source_aligned_endpoint_target_worker", TARGET_SHA256
        )
    return _TARGET_WORKER_MODULE.worker(alpha)


def build_pair_weight(target) -> tuple[list[list[int]], tuple[object, ...]]:
    """Independently reconstruct E_(omega,nu) from THM-3514 primitives."""
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(target_worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")
    tables = tuple(row for chunk in chunks for row in chunk[1])
    require(len(tables) == P**2, "endpoint table count")

    (
        _word,
        _t_den,
        _nn,
        prime,
        _root,
        zeta,
        _q_intervals,
        _q_starts,
        _embeddings,
        _tabs,
        _atom_intervals,
    ) = target.context()
    require(pow(zeta, P, prime) == 1 and zeta != 1, "target zeta order")

    h_kernel = [[0] * len(target.ATOMS) for _left in target.ATOMS]
    q5_kernel = [[0] * len(target.ATOMS) for _left in target.ATOMS]
    for left_index, (left_sheet, left_chamber) in enumerate(target.ATOMS):
        for right_index, (right_sheet, right_chamber) in enumerate(target.ATOMS):
            products = tuple(
                target.safe(left_chamber, left_sheet + tau)
                * target.safe(right_chamber, right_sheet + tau)
                for tau in range(P)
            )
            q5_kernel[left_index][right_index] = sum(products) % prime
            h_kernel[left_index][right_index] = sum(
                value * pow(zeta, -tau % P, prime)
                for tau, value in enumerate(products)
            ) % prime

    pair_weight = [[0] * len(target.ATOMS) for _left in target.ATOMS]
    table_index = 0
    for alpha in range(P):
        alpha_phase = pow(zeta, -alpha % P, prime)
        for beta in range(P):
            stored_beta, ax_values, by_values = tables[table_index]
            table_index += 1
            require(stored_beta == beta, (alpha, beta, stored_beta))
            phase = pow(zeta, beta, prime) * alpha_phase % prime
            for left_index, left_value in enumerate(ax_values):
                if not left_value:
                    continue
                for right_index, right_value in enumerate(by_values):
                    if not right_value:
                        continue
                    kernel = (
                        h_kernel[left_index][right_index]
                        - q5_kernel[left_index][right_index]
                    ) % prime
                    pair_weight[left_index][right_index] = (
                        pair_weight[left_index][right_index]
                        + phase * left_value % prime * right_value % prime * kernel
                    ) % prime
    normalizer = pow(P**3, -1, prime)
    pair_weight = [[value * normalizer % prime for value in row]
                   for row in pair_weight]
    bridge = sum(sum(row) for row in pair_weight) % prime
    require(bridge == target.EXPECTED_BRIDGE, ("endpoint bridge", bridge))
    record = (
        prime,
        zeta,
        bridge,
        sum(value != 0 for row in pair_weight for value in row),
        digest_json(pair_weight),
    )
    return pair_weight, record


def decompose_pair_function(
    pair_weight: list[list[int]], prime: int
) -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    """Return left-only, right-only, and doubly-centred ANOVA pieces."""
    size = len(pair_weight)
    inv_size = pow(size, -1, prime)
    inv_square = inv_size * inv_size % prime
    row_sums = [sum(row) % prime for row in pair_weight]
    col_sums = [sum(pair_weight[left][right] for left in range(size)) % prime
                for right in range(size)]
    total = sum(row_sums) % prime
    grand = total * inv_square % prime
    left_only = [[(row_sums[left] * inv_size - grand) % prime] * size
                 for left in range(size)]
    right_only = [[(col_sums[right] * inv_size - grand) % prime
                   for right in range(size)] for _left in range(size)]
    interaction = [[
        (pair_weight[left][right]
         - row_sums[left] * inv_size
         - col_sums[right] * inv_size
         + grand) % prime
        for right in range(size)
    ] for left in range(size)]
    require(all(sum(row) % prime == 0 for row in interaction),
            "interaction row centering")
    require(all(sum(interaction[left][right] for left in range(size)) % prime == 0
                for right in range(size)), "interaction column centering")
    return left_only, right_only, interaction


def pull_profile(
    pair_gauge_offset: list[list[list[int]]],
    pair_function: list[list[int]],
    denominator: int,
    zeta: int,
    prime: int,
) -> tuple[int, ...]:
    inverse_denominator = pow(denominator, -1, prime)
    answer = []
    for owner_frequency in range(P):
        value = 0
        for left in range(len(pair_function)):
            for right in range(len(pair_function)):
                source_weight = sum(
                    pair_gauge_offset[left][right][offset]
                    * pow(zeta, -owner_frequency * offset % P, prime)
                    for offset in range(P)
                ) % prime
                value = (
                    value + source_weight * pair_function[left][right]
                ) % prime
        answer.append(value * inverse_denominator % prime)
    return tuple(answer)


def main() -> None:
    source = load_module(SOURCE_PATH, "source_aligned_guard_ancestry", SOURCE_SHA256)
    with redirect_stdout(io.StringIO()):
        source_data = source.main()
    require(source_data["semantic_sha256"] == SOURCE_SEMANTIC, "source semantic drift")

    target = load_module(TARGET_PATH, "source_aligned_endpoint_target", TARGET_SHA256)
    pair_weight, target_record = build_pair_weight(target)
    prime, zeta, bridge, endpoint_pair_support, pair_digest = target_record
    require(pair_digest == "c2d5911b287510335edc6aefa6d3b865c982568f678bb89ee9b82ee211962df1",
            ("endpoint pair digest", pair_digest))

    pair_gauge_offset = source_data["pair_gauge_offset"]
    denominator = source_data["denominator"]
    require(denominator % prime != 0, "source denominator at target prime")
    source_pair_support = sum(
        any(pair_gauge_offset[left][right])
        for left in range(len(pair_gauge_offset))
        for right in range(len(pair_gauge_offset))
    )
    weighted_pair_support = sum(
        any(pair_gauge_offset[left][right]) and pair_weight[left][right] != 0
        for left in range(len(pair_gauge_offset))
        for right in range(len(pair_gauge_offset))
    )
    require((source_pair_support, weighted_pair_support) == (362, 362),
            ("pair support intersection", source_pair_support, weighted_pair_support))

    full_profile = pull_profile(
        pair_gauge_offset, pair_weight, denominator, zeta, prime
    )
    left_only, right_only, interaction = decompose_pair_function(pair_weight, prime)
    left_profile = pull_profile(
        pair_gauge_offset, left_only, denominator, zeta, prime
    )
    right_profile = pull_profile(
        pair_gauge_offset, right_only, denominator, zeta, prime
    )
    interaction_profile = pull_profile(
        pair_gauge_offset, interaction, denominator, zeta, prime
    )
    constant_function = [[1] * len(pair_weight) for _left in pair_weight]
    flat_profile = pull_profile(
        pair_gauge_offset, constant_function, denominator, zeta, prime
    )
    require(all(value != 0 for value in full_profile),
            ("source-aligned endpoint pullback", full_profile))
    require(all(value != 0 for value in left_profile),
            ("left-only endpoint ANOVA profile", left_profile))
    require(all(value != 0 for value in right_profile),
            ("right-only endpoint ANOVA profile", right_profile))
    require(all(value != 0 for value in interaction_profile),
            ("doubly-centred endpoint interaction profile", interaction_profile))
    require(all(value != 0 for value in flat_profile),
            ("flat simple-kernel control", flat_profile))

    record = (
        SOURCE_SHA256,
        SOURCE_SEMANTIC,
        TARGET_SHA256,
        TARGET_BUCKET_SHA256,
        target_record,
        source_pair_support,
        weighted_pair_support,
        full_profile,
        tuple(int(value != 0) for value in full_profile),
        left_profile,
        right_profile,
        interaction_profile,
        flat_profile,
        digest_json(interaction),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 source-aligned actual endpoint simple-kernel probe ==")
    print(f"dependencies=((source,{SOURCE_SHA256},{SOURCE_SEMANTIC}),(target,{TARGET_SHA256},{TARGET_BUCKET_SHA256}))")
    print(f"target_split_field=(prime={prime},zeta13={zeta})")
    print(f"target_endpoint_pair=(support={endpoint_pair_support}/1521,bridge={bridge},sha256={pair_digest})")
    print(f"source_common_gauge_pair_support={source_pair_support}/1521; actual_endpoint_weighted_support={weighted_pair_support}/1521")
    print("typing=both guard labels read at endpoint/source time via f_omega^src=1_Q P^K(e P_omega), then meet on one THM-2471 base")
    print(f"actual_endpoint_simple_kernel_by_owner_frequency={full_profile}")
    print(f"actual_endpoint_simple_kernel_support={sum(value != 0 for value in full_profile)}/13")
    print(f"endpoint_ANOVA_left_only_profile={left_profile}")
    print(f"endpoint_ANOVA_right_only_profile={right_profile}")
    print(f"endpoint_ANOVA_doubly_centred_interaction_profile={interaction_profile}")
    print("endpoint_ANOVA_profile_support=(left_only=13/13,right_only=13/13,doubly_centred_interaction=13/13)")
    print(f"flat_simple_kernel_profile={flat_profile}")
    print(f"semantic_sha256={semantic}")
    print("scope=one-base source-time finite simple-kernel transplant; AX/BY remain preintegrated scalars, so no physical relation current,grouped coefficient,row exclusion,LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
