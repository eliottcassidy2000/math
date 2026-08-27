#!/usr/bin/env python3
"""Postprocess THM-4278's 72 repairs against THM-4271's full prefix."""

from collections import Counter
from math import gcd
from pathlib import Path
import re
import sys


POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)


def require(condition, message="postprocess failure"):
    if not condition:
        raise RuntimeError(message)


def parse_hex_line(text, label, separator=None):
    match = re.search(rf"(?m)^{re.escape(label)} (.+)$", text)
    require(match is not None, f"missing {label}")
    payload = match.group(1)
    words = payload.split(separator) if separator is not None else payload.split()
    return tuple(int(word, 16) for word in words if word)


def labels(mask):
    return tuple(POOL[bit] for bit in range(30) if mask & (1 << bit))


def main():
    if len(sys.argv) == 1:
        literal_path = Path(
            "05-knowledge/results/"
            "lrc14_endpoint_520_688_minimum_one_atom_literal_wall_audit_thm4278.out"
        )
        prefix_path = Path(
            "05-knowledge/results/lrc14_fourth_round_learned_carrier_thm4271/"
            "full_discovery_520_688_O3.semantic.out"
        )
    elif len(sys.argv) == 3:
        literal_path = Path(sys.argv[1])
        prefix_path = Path(sys.argv[2])
    else:
        raise RuntimeError("usage: postprocess [LITERAL_OUT THM4271_DISCOVERY_OUT]")

    literal = literal_path.read_text(encoding="utf-8")
    discovery = prefix_path.read_text(encoding="utf-8")
    repairs = parse_hex_line(literal, "COMMON_ACTIVE_MASKS_HEX")
    prefix = parse_hex_line(discovery, "PREFIX_MASKS_HEX", separator=",")

    require(len(repairs) == len(set(repairs)) == 72, "repair family changed")
    require(len(prefix) == len(set(prefix)) == 5398, "THM-4271 prefix changed")
    require(all(mask.bit_count() == 8 for mask in repairs), "non-rank-eight repair")
    overlap = sorted(set(repairs) & set(prefix))
    require(not overlap, "targeted repairs entered THM-4271 prefix")
    require(repairs[0] == 0x00048EC1, "least repair changed")

    gcd_profile = Counter()
    odd_profile = Counter()
    residue_profile = Counter()
    for mask in repairs:
        values = labels(mask)
        common_gcd = 0
        for value in values:
            common_gcd = gcd(common_gcd, value)
        gcd_profile[common_gcd] += 1
        odd_profile[sum(value & 1 for value in values)] += 1
        residue_profile[sum(values) % 14] += 1

    require(set(gcd_profile) == {1, 2}, "unexpected gcd profile")
    require(set(odd_profile) == {0, 1, 2, 3}, "unexpected odd-count profile")
    require(set(residue_profile) == set(range(14)), "sum residues not exhaustive")

    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    print("lrc14_endpoint_520_688_one_atom_overlap_postprocess_thm4278")
    print(f"common_active_repairs={len(repairs)}")
    print(f"thm4271_full_prefix={len(prefix)}")
    print(f"prefix_overlap={len(overlap)}")
    print("thm4271_novel_overlap=0")
    print(f"least_repair=0x{repairs[0]:08x}")
    print(f"gcd_profile={dict(sorted(gcd_profile.items()))}")
    print(f"odd_label_count_profile={dict(sorted(odd_profile.items()))}")
    print(f"sum_mod14_profile={dict(sorted(residue_profile.items()))}")
    print("no_single_gcd_parity_or_sum_mod14_classifier=PASS")
    print("current_proof_graph_contribution_after_thm4271=0")
    print("status=PASS")


if __name__ == "__main__":
    main()
