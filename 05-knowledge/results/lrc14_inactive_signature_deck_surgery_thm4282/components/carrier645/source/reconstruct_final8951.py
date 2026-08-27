#!/usr/bin/env python3
"""Independent Python reconstruction of the ordered THM-4281 carrier."""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import re


FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
EXPECTED_PAIRS = [
    (616,755),(616,756),(616,757),(616,758),(616,759),(616,760),
    (616,761),(616,762),(616,763),(616,764),(616,765),(616,766),
    (616,767),(616,768),(698,755),(698,757),(704,755),(704,757),
    (704,758),(704,759),(704,761),(704,762),(704,763),(704,764),
    (704,765),(721,755),(721,757),(721,758),(721,759),(721,761),
    (721,762),(721,763),(721,764),(721,765),(721,766),(721,767),
    (721,768),(726,755),(726,757),(726,758),(726,761),(732,755),
    (732,757),(732,761),(732,762),(732,763),(744,762),(744,763),
    (744,765),(744,766),(744,768),(750,762),(750,763),(750,765),
    (750,766),(750,768),(765,766),(765,768),(766,768),
]
PREFIX_CONTROLS = {
    "prefix704": ((416,704), 2608, 0x18FF663A123E684E),
    "prefix416": ((416,700), 5894, 0x701C0233C8F8ABEB),
    "prefix520": ((520,700), 4557, 0xD8466558BCD8EF9D),
    "prefix384": ((384,694), 5307, 0x7B9A46C60E8514F8),
    "prefix688": ((520,688), 5398, 0x6AB471DA88C8E1D1),
}
THM4276 = [0x00289285,0x0260812C,0x18689040,
           0x20C0C124,0x302C1006,0x30580888]
THM4281_CONTINUATION = [0x003884C8,0x00C4C124,0x02A05206,
                        0x10E05240,0x0004409F,0x00C0C125]


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fnv(values: list[int]) -> int:
    state = FNV_OFFSET
    for value in values:
        need(0 <= value < 1 << 64, "FNV value escaped u64")
        for byte in value.to_bytes(8, "little"):
            state ^= byte
            state = state * FNV_PRIME & ((1 << 64) - 1)
    return state


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_probe(path: pathlib.Path) -> tuple[tuple[int, int], list[int], int]:
    pair: tuple[int, int] | None = None
    stated_size: int | None = None
    stated_fnv: int | None = None
    masks: list[int] | None = None
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("PAIR "):
            match = re.match(r"PAIR ([0-9]+),([0-9]+) ", line)
            need(match is not None, f"malformed PAIR row in {path}")
            pair = (int(match.group(1)), int(match.group(2)))
        elif line.startswith("PREFIX_CERT "):
            fields = line.split()
            need(len(fields) >= 5 and fields[:2] == ["PREFIX_CERT", "SIZE"]
                 and fields[3] == "FNV", f"malformed PREFIX_CERT in {path}")
            stated_size = int(fields[2])
            stated_fnv = int(fields[4], 16)
        elif line.startswith("PREFIX_MASKS_HEX "):
            need(masks is None, f"duplicate PREFIX_MASKS_HEX in {path}")
            masks = [int(token, 16) for token in
                     line.removeprefix("PREFIX_MASKS_HEX ").split(",")]
    need(pair is not None and stated_size is not None and stated_fnv is not None
         and masks is not None, f"incomplete probe transcript {path}")
    need(len(masks) == stated_size and len(set(masks)) == len(masks),
         f"prefix size/distinctness changed in {path}")
    need(all(mask < 1 << 30 and mask.bit_count() == 8 for mask in masks),
         f"prefix has invalid rank-eight mask in {path}")
    need(fnv(masks) == stated_fnv, f"prefix FNV changed in {path}")
    return pair, masks, stated_fnv


def read_joint(path: pathlib.Path) -> list[int]:
    masks = [int(line, 16) for line in path.read_text(encoding="ascii").splitlines()
             if line]
    need(len(masks) == 421 and len(set(masks)) == 421 and
         all(mask < 1 << 30 and mask.bit_count() == 8 for mask in masks) and
         fnv(masks) == 0x20D63DD42FE8150E,
         "joint deck identity changed")
    return masks


def append_unseen(target: list[int], seen: set[int], source: list[int]) -> None:
    for mask in source:
        if mask not in seen:
            seen.add(mask)
            target.append(mask)


def stage(name: str, carrier: list[int], count: int, expected_fnv: int) -> str:
    actual = fnv(carrier)
    need(len(carrier) == count and actual == expected_fnv,
         f"{name} count/FNV changed")
    return f"STAGE {name} COUNT {count} FNV {actual:016x}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--band", type=pathlib.Path, required=True)
    for name in PREFIX_CONTROLS:
        parser.add_argument(f"--{name}", type=pathlib.Path, required=True)
    parser.add_argument("--joint", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--compare", type=pathlib.Path)
    args = parser.parse_args()

    expected_names = {f"{q}_{r}.out" for q, r in EXPECTED_PAIRS}
    actual_names = {path.name for path in args.band.glob("*.out")}
    need(actual_names == expected_names, "59-file replay-band identity changed")
    carrier: list[int] = []
    seen: set[int] = set()
    for expected in EXPECTED_PAIRS:
        pair, masks, _ = parse_probe(args.band / f"{expected[0]}_{expected[1]}.out")
        need(pair == expected, "replay-band pair/path mismatch")
        append_unseen(carrier, seen, masks)
    rows = [stage("BASE4254", carrier, 4675, 0xCE4E76EC11DF057C)]

    parsed: dict[str, list[int]] = {}
    for name, (expected_pair, expected_count, expected_fnv) in PREFIX_CONTROLS.items():
        pair, masks, stated_fnv = parse_probe(getattr(args, name))
        need(pair == expected_pair and len(masks) == expected_count and
             stated_fnv == expected_fnv, f"{name} identity changed")
        parsed[name] = masks

    append_unseen(carrier, seen, parsed["prefix704"])
    rows.append(stage("ROUND1", carrier, 4733, 0xA7B046289655C733))
    round1 = set(seen)
    novel416 = [mask for mask in parsed["prefix416"] if mask not in round1]
    novel520 = [mask for mask in parsed["prefix520"] if mask not in round1]
    append_unseen(carrier, seen, novel416)
    append_unseen(carrier, seen, novel520)
    rows.append(stage("ROUND2", carrier, 7986, 0xBAEF1D2F49444638))
    append_unseen(carrier, seen, parsed["prefix384"])
    rows.append(stage("ROUND3", carrier, 8319, 0xE08B227730F6793C))
    append_unseen(carrier, seen, parsed["prefix688"])
    rows.append(stage("THM4271", carrier, 8518, 0x1603E3FE970F8428))
    append_unseen(carrier, seen, THM4276)
    rows.append(stage("THM4276", carrier, 8524, 0x5DDB84A44F5D2AD7))

    joint = read_joint(args.joint)
    need(not (set(joint) & seen), "joint deck overlaps THM-4276 carrier")
    append_unseen(carrier, seen, joint)
    rows.append(stage("JOINT421", carrier, 8945, 0x3212EFA05DD18C00))
    need(not (set(THM4281_CONTINUATION) & seen) and
         len(set(THM4281_CONTINUATION)) == 6 and
         all(mask.bit_count() == 8 for mask in THM4281_CONTINUATION),
         "THM-4281 continuation identity/novelty changed")
    append_unseen(carrier, seen, THM4281_CONTINUATION)
    rows.append(stage("FINAL4281", carrier, 8951, 0x188F82AB9DD1695A))
    need(len(seen) == len(carrier), "final carrier repeats a mask")

    args.output.write_text("".join(f"{mask:08x}\n" for mask in carrier),
                           encoding="ascii", newline="\n")
    output_sha = sha256(args.output)
    need(output_sha == "a5dac3c7e5a2715e4c9ef8bb1b54bc98792e904f7b0b5ef55e4dd4313ebc87f6",
         "reconstructed carrier SHA changed")
    if args.compare is not None:
        need(args.output.read_bytes() == args.compare.read_bytes(),
             "reconstructed carrier differs from dormant supplied copy")

    print("THM4282_FINAL8951_INDEPENDENT_PYTHON_RECONSTRUCTION_V1")
    print(*rows, sep="\n")
    print(f"FINAL_SHA256 {output_sha} COMPARE_IDENTICAL {int(args.compare is not None)}")
    print("VERDICT PASS EXACT_ORDERED_CARRIER_RECONSTRUCTION")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        import sys
        print(f"RECONSTRUCTION_ERROR {error}", file=sys.stderr)
        raise SystemExit(1)
