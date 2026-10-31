"""Explicit uniform-weight examples and independent extension DP controls."""

import hashlib
import json
from pathlib import Path

from cut_response import compile_response, extension
from verify import direct_h


def main():
    rows = []
    direct_checks = 0
    for n, expected_h, expected_w in [(7, 189, 39), (11, 95095, 14421)]:
        residues = {x * x % n for x in range(1, n)}
        out = tuple(sum(1 << v for v in range(n) if (v - u) % n in residues)
                    for u in range(n))
        response = compile_response(out)
        if response.h != expected_h or direct_h(out) != expected_h:
            raise RuntimeError(("base count", n))
        if any(row.bit_count() != (n - 1) // 2 for row in out):
            raise RuntimeError(("regularity", n))
        if any(response.weights[u][v] != expected_w
               for u in range(n) for v in range(u + 1, n)):
            raise RuntimeError(("constant weight", n))
        values = response.all_values()
        for cut, value in enumerate(values):
            size = cut.bit_count()
            if value != expected_h + expected_w * size * (n - size):
                raise RuntimeError(("size response", n, cut))
            if n == 7 or cut == (1 << size) - 1:
                if value != direct_h(extension(out, cut)):
                    raise RuntimeError(("independent extension count", n, cut))
                direct_checks += 1
        if 2 * expected_h - 1 in values or 2 * expected_h + 1 in values:
            raise RuntimeError(("binary target exclusion", n))
        rows.append({"n": n, "out": out, "h": expected_h, "pair_weight": expected_w,
                     "gain_gcd": response.gain_content(),
                     "extension_values": sorted(set(values)),
                     "minimum_positive": response.minimum_positive()})
    print(json.dumps({"status": "PASS", "examples": rows,
                      "independent_extension_dp_checks": direct_checks,
                      "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                     sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
