#!/usr/bin/env python3
"""Exact carry-response and scalar-amplitude audit for THM-3657's quotient."""

from __future__ import annotations

import ast
from collections import defaultdict
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT_SCRIPT = ROOT / "04-computation/lrc_two_current_quotient_address_reversal_gate_thm3657.py"
PARENT_OUTPUT = ROOT / "05-knowledge/results/lrc_two_current_quotient_address_reversal_gate_thm3657.out"
EXPECTED_PARENT_HASHES = (
    "f0323550a039bd3c59bc3367a9a48503ff10db2b642e99e21d9499e984492ccd",
    "fdebbbba161149edec8c86b9f382ede76153e85aa01da635f0900075fe751628",
)
EXPECTED_SEMANTIC_SHA256 = "31a8f50ea2d5a9dbbc52403a50dd49594addc45baa3de9c0043752db64a1fac8"

P = 13
N = 169
MOD = 755373809845391722745761
OMEGA = 298763986285447441216949
INTERVALS = ((13, 24), (26, 47), (49, 70), (72, 77),
             (91, 96), (98, 119), (121, 142), (144, 155))
EXCEPTIONAL_A = frozenset((12, 25, 48, 71, 97, 120, 143, 156))


def require(condition: bool, payload: object) -> None:
    if condition is not True:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def load_parent():
    hashes = (lf_sha256(PARENT_SCRIPT), lf_sha256(PARENT_OUTPUT))
    require(hashes == EXPECTED_PARENT_HASHES, ("THM-3657 parent hashes", hashes))
    spec = importlib.util.spec_from_file_location("thm3657_parent", PARENT_SCRIPT)
    require(spec is not None and spec.loader is not None, "THM-3657 loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module, hashes


def add(left, right):
    return tuple((a + b) % MOD for a, b in zip(left, right))


def subtract(*vectors):
    first, *rest = vectors
    return tuple((first[index] - sum(vector[index] for vector in rest)) % MOD
                 for index in range(len(first)))


def split_add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def cyclic_add(left, right):
    total = (left[0] + P * left[1] + right[0] + P * right[1]) % N
    return total % P, total // P


def polynomial_degree(values):
    current = tuple(value % MOD for value in values)
    for order in range(len(current)):
        if len(set(current)) <= 1:
            return order
        current = tuple((current[index + 1] - current[index]) % MOD
                        for index in range(len(current) - 1))
    return len(values) - 1


def cyclic_support(values):
    support = []
    for frequency in range(N):
        phase = pow(OMEGA, -frequency % N, MOD)
        running = 1
        total = 0
        for value in values:
            total = (total + value * running) % MOD
            running = running * phase % MOD
        if total:
            support.append(frequency)
    return tuple(support)


def main() -> None:
    T, parent_hashes = load_parent()
    M = T.M
    two_tensor, reconstruction = M.reconstruct_two_current()
    require(reconstruction == M.EXPECTED_PARENT_RECONSTRUCTION_SHA256[:2],
            "two-current reconstruction drift")
    raw = M.flatten_two_current(two_tensor)
    correction = tuple(T.state_side(row) for row in raw)
    basis = M.canonical_row_basis(correction)
    require(len(basis) == 2, "correction rank")
    basis_digest = M.rowspace_digest(basis)
    require(basis_digest == T.G.EXPECTED_DIGESTS["state_side"],
            "correction basis digest")
    coordinates = tuple(T.coordinates_in_rref(basis, row) for row in correction)
    labels = tuple((r0, r1) for r0 in range(P) for r1 in range(P))
    by_label = dict(zip(labels, coordinates))

    split_defects = []
    cyclic_defects = []
    carry_defects = []
    carry_bits = []
    pair_records = []
    for left in labels:
        for right in labels:
            split_sum = split_add(left, right)
            cyclic_sum = cyclic_add(left, right)
            split_defect = subtract(by_label[split_sum], by_label[left], by_label[right])
            cyclic_defect = subtract(by_label[cyclic_sum], by_label[left], by_label[right])
            carry_defect = subtract(cyclic_defect, split_defect)
            carry = int(left[0] + right[0] >= P)
            require((cyclic_sum == split_sum) == (carry == 0),
                    ("carry/sum equivalence", left, right))
            if carry == 0:
                require(carry_defect == (0, 0), ("zero-carry defect", left, right))
            else:
                expected_sum = (split_sum[0], (split_sum[1] + 1) % P)
                require(cyclic_sum == expected_sum, ("carry increment", left, right))
                require(carry_defect
                        == subtract(by_label[expected_sum], by_label[split_sum]),
                        ("carry response factor", left, right))
            split_defects.append(split_defect)
            cyclic_defects.append(cyclic_defect)
            carry_defects.append(carry_defect)
            carry_bits.append(carry)
            pair_records.append((left, right, split_sum, cyclic_sum, carry_defect))

    def defect_record(rows):
        distinct = tuple(sorted(set(rows)))
        return (
            M.rank_mod(tuple(rows)),
            sum(row != (0, 0) for row in rows),
            len(distinct),
            digest_json(tuple(rows)),
            digest_json(distinct),
        )

    split_record = defect_record(split_defects)
    cyclic_record = defect_record(cyclic_defects)
    carry_record = defect_record(carry_defects)
    carry_one = tuple(row for row, bit in zip(carry_defects, carry_bits) if bit)
    carry_one_record = defect_record(carry_one)
    require(split_record[:3] == (2, 26474, 549), ("split record", split_record))
    require(cyclic_record[:3] == (2, 26079, 542), ("cyclic record", cyclic_record))
    require(carry_one_record[:3] == (2, 11856, 55),
            ("carry-one record", carry_one_record))

    # The carry response factors through the split sum t=(t0,t1), with t0<12.
    response_table = tuple(
        ((t0, t1), subtract(by_label[(t0, (t1 + 1) % P)], by_label[(t0, t1)]))
        for t0 in range(P - 1) for t1 in range(P)
    )
    response_rows = tuple(row for _label, row in response_table)
    response_record = defect_record(response_rows)
    require(response_record[0] == 2 and response_record[2] == 55,
            ("response table record", response_record))
    require(set(response_rows) == set(carry_one), "carry response image")
    weighted_nonzero = sum(
        13 * (12 - t0) * int(row != (0, 0))
        for (t0, _t1), row in response_table
    )
    require(weighted_nonzero == carry_one_record[1],
            ("carry multiplicity invoice", weighted_nonzero))

    # The generic +1 line has a finite 16-value scalar automaton, but no
    # low-degree or cyclic-Fourier sparsity in the assembled coordinate.
    assembled = tuple(by_label[(a % P, a // P)] for a in range(N))
    zero_a = frozenset(range(0, 12)) | frozenset(range(78, 91)) | frozenset(range(157, 169))
    generic_a = tuple(a for a in range(N) if a not in zero_a and a not in EXCEPTIONAL_A)
    require(generic_a == tuple(a for left, right in INTERVALS
                               for a in range(left, right + 1)),
            "generic interval assembly")
    amplitudes = tuple(assembled[a][0] for a in generic_a)
    require(len(set(amplitudes)) == 16, "generic amplitude count")
    require(all(assembled[a][0] == assembled[168 - a][0] for a in generic_a),
            "generic amplitude reversal")
    interval_degrees = tuple(
        polynomial_degree(tuple(assembled[a][0] for a in range(left, right + 1)))
        for left, right in INTERVALS
    )
    require(interval_degrees == (11, 21, 21, 5, 5, 21, 21, 11),
            ("generic interpolation degrees", interval_degrees))
    coordinate_sequences = tuple(tuple(row[index] for row in assembled)
                                 for index in range(2))
    cyclic_supports = tuple(cyclic_support(sequence) for sequence in coordinate_sequences)
    require(tuple(len(support) for support in cyclic_supports) == (169, 169),
            "cyclic Fourier support")
    amplitude_groups = defaultdict(list)
    for a in generic_a:
        amplitude_groups[assembled[a][0]].append(abs(a - 84))
    amplitude_group_record = tuple(
        (value, tuple(sorted(set(distances))))
        for value, distances in sorted(amplitude_groups.items(), key=lambda item: min(item[1]))
    )

    semantic = digest_json((
        MOD, P, N, OMEGA, parent_hashes, reconstruction, basis_digest,
        split_record, cyclic_record, carry_record, carry_one_record,
        response_record, digest_json(response_table), weighted_nonzero,
        interval_degrees, tuple(map(len, cyclic_supports)),
        digest_json(cyclic_supports), amplitude_group_record,
        digest_json(pair_records),
    ))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic, EXPECTED_SEMANTIC_SHA256))

    source = Path(__file__).resolve()
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source.read_text(encoding="utf-8")))),
            "Python assert node present")
    print("== THM-3659 LRC mod-169 carry-response cochain obstruction ==")
    print(f"field=p:{MOD};addresses:169;ordered_pairs:28561")
    print(f"parent_sha256_lf={parent_hashes};reconstruction_sha256={reconstruction}")
    print(f"correction_plane=rank:2;rref_sha256:{basis_digest}")
    print(f"split_defect=(rank,nonzero,distinct,ledger_sha256,set_sha256)={split_record}")
    print(f"cyclic_defect=(rank,nonzero,distinct,ledger_sha256,set_sha256)={cyclic_record}")
    print(f"carry_defect_all=(rank,nonzero,distinct,ledger_sha256,set_sha256)={carry_record}")
    print(f"carry_defect_kappa1=(rank,nonzero,distinct,ledger_sha256,set_sha256)={carry_one_record}")
    print(f"carry_response_table=(rank,nonzero,distinct,ledger_sha256,set_sha256)={response_record};table_sha256={digest_json(response_table)}")
    print(f"carry_factor=Gamma(x,y)=e(t0,t1+1)-e(t0,t1);weighted_nonzero_invoice={weighted_nonzero}")
    print(f"generic_scalar_automaton=addresses:124;distinct_values:16;reversal_even:True;group_sha256:{digest_json(amplitude_group_record)}")
    print(f"generic_interval_degrees={interval_degrees};cyclic_fourier_support_sizes={tuple(map(len, cyclic_supports))};support_sha256={digest_json(cyclic_supports)}")
    print("one_bit_carry_sidecar=REFUTED;carry1 response has 55 values and rank 2")
    print(f"semantic_sha256={semantic}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("scope=static finite-field response cochain;not current/chronology/entry/characteristic-zero/LRC14")
    print("PASS")


if __name__ == "__main__":
    main()
