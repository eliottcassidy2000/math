"""Export and pin THM-4231's exact 181,194-edge patched remainder."""

import contextlib
import io
import pathlib
import runpy
import sys

SOURCE = pathlib.Path(
    "04-computation/lrc14_literal_boundary_45_residual_postprocess_thm4231.py"
)


def fnv_add(value: int, word: int) -> int:
    for shift in range(0, 64, 8):
        value ^= (word >> shift) & 0xFF
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


def edge_fnv(edges: list[tuple[int, int]]) -> int:
    value = 0xCBF29CE484222325
    for q, r in edges:
        value = fnv_add(value, q)
        value = fnv_add(value, r)
    return value


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: export_thm4231_remainder.py OUTPUT")
    sink = io.StringIO()
    with contextlib.redirect_stdout(sink):
        namespace = runpy.run_path(str(SOURCE))
    remainder = namespace["remainder"]
    assert len(remainder) == 181_194
    assert edge_fnv(remainder) == 0x3874FECAC4ECBD8A
    assert max(r for _, r in remainder) == 769
    output = pathlib.Path(sys.argv[1])
    output.write_text("".join(f"{q},{r}\n" for q, r in remainder), encoding="ascii")
    print(
        "THM4231_REMAINDER_EXPORT PASS "
        f"count={len(remainder)} fnv={edge_fnv(remainder):016x} "
        f"max_endpoint={max(r for _, r in remainder)}"
    )


if __name__ == "__main__":
    main()
