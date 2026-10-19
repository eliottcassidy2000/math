#!/usr/bin/env python3
"""Independent direct-leg and direct-root audit; standard library only.

The two censuses below use neither Euclid parameters nor prime-factor counts.
All checks remain active under python -O. Importing the companion creates no
bytecode cache, and running this file writes nothing.
"""

import importlib.util
import json
from math import gcd, isqrt
from pathlib import Path
import sys

sys.dont_write_bytecode = True
source = Path(__file__).with_name("endpoint_pythagorean_fibers_20261019.py")
spec = importlib.util.spec_from_file_location("endpoint_bridge", source)
bridge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)

CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def run(limit=500):
    direct_legs = {}
    direct_ap = {}
    for c in range(1, limit + 1):
        triples = set()
        for a in range(1, c):
            b = isqrt(c * c - a * a)
            if a <= b and a * a + b * b == c * c and gcd(a, b) == 1:
                triples.add((a, b, c) if a % 2 else (b, a, c))
        direct_legs[c] = triples

        # Directly enumerate the pair-sum roots p < q < r, with q=2C.
        # The label gcd two selects exactly the primitive/base scale.
        labels = set()
        q = 2 * c
        for p in range(1, q):
            r = isqrt(2 * q * q - p * p)
            if p * p + r * r != 2 * q * q:
                continue
            numerators = (p * p + q * q - r * r,
                          p * p + r * r - q * q,
                          q * q + r * r - p * p)
            if any(v % 2 for v in numerators):
                continue
            row = tuple(v // 2 for v in numerators)
            if gcd(*row) == 2:
                labels.add(row)
        direct_ap[c] = labels

    euclid = bridge.euclid_fibers(limit)
    admissible_norms = 0
    for c in range(1, limit + 1):
        check(direct_legs[c] == euclid[c], f"direct-leg/Euclid census at {c}")
        try:
            data = bridge.fiber_data(c)
        except ValueError:
            check(not direct_legs[c], f"excluded norm {c}")
        else:
            admissible_norms += 1
            gaussian = {bridge.ppt_from_alpha(bridge.alpha_from_mask(data, mask))
                        for mask in range(1 << len(data))}
            check(gaussian == direct_legs[c], f"direct-leg/Gaussian census at {c}")
            for triple in direct_legs[c]:
                address = bridge.canonical_mask_from_ppt(data, triple)
                check(address & 1, "canonical complement representative")
                check(bridge.ppt_from_alpha(bridge.alpha_from_mask(data, address)) == triple,
                      "inverse round trip from independently found triple")

        expected_ap = set()
        for a, b, _ in direct_legs[c]:
            row = (2 * (c * c - 4 * a * b), 2 * c * c,
                   2 * (c * c + 4 * a * b))
            x, y, z = row
            check(x < y < z and x + z == 2 * y, "strict signed AP")
            check((x + y, x + z, y + z) ==
                  (4 * (a - b) ** 2, 4 * c * c, 4 * (a + b) ** 2),
                  "all pair-sum squares")
            check(gcd(x, y) == gcd(x, z) == gcd(y, z) == 2,
                  "base scale and pairwise label gcd")
            expected_ap.add(row)
        check(expected_ap == direct_ap[c], f"independent signed-AP root census at {c}")
        check(len(expected_ap) == len(direct_legs[c]), "injective primitive AP map")

    return {"status": "PASS", "checks": CHECKS, "hypotenuse_limit": limit,
            "admissible_norms": admissible_norms,
            "direct_leg_primitive_triples": sum(map(len, direct_legs.values())),
            "direct_root_base_signed_ap_triangles": sum(map(len, direct_ap.values()))}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
