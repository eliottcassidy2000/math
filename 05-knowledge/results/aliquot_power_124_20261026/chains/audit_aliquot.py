#!/usr/bin/env python3
"""Independent addition-sieve audit of the neighboring aliquot receipt."""
import hashlib
import json
from math import isqrt
from pathlib import Path

base = Path(__file__).resolve().parent.parent / "aliquot"
producer = json.loads((base / "results.normal.json").read_text())
bound = 3481
sums = [0] * (bound + 1)
prime = [True] * (bound + 1)
prime[0] = prime[1] = False
for d in range(1, bound // 2 + 1):
    for n in range(2 * d, bound + 1, d):
        sums[n] += d
        if d > 1:
            prime[n] = False

gates = 0


def require(test, message):
    global gates
    gates += 1
    if not test:
        raise RuntimeError(message)


fibres = {str(m): [] for m in range(2, 61)}
for n in range(1, bound + 1):
    m = sums[n]
    if str(m) in fibres:
        fibres[str(m)].append(n)
    require((m == 0) == (n == 1), ("zero", n))
    require((m == 1) == prime[n], ("prime", n))
    if m > 1:
        require(n <= (m - 1) ** 2, ("inverse-bound", n))
        root = isqrt(n)
        square_prime = root * root == n and prime[root]
        require((n == (m - 1) ** 2) == square_prime, ("equality", n))
        require(square_prime or 4 * n <= (m - 1) ** 2, ("strong-bound", n))
for m in fibres:
    require(fibres[m] == producer["fibres"][m], ("fibre", m))
require([int(m) for m, f in fibres.items() if not f] == [2, 5, 52], "untouchables")
require(fibres["6"] == [6, 25] and sums[6] == 6, "unique6-control")
require(fibres["60"][-1] == 3481, "late-source")
require(sums[729] == 364 != 28, "composite27-hostile")
require(sums[9] == 4 != 7, "repeated-prime-hostile")

print(json.dumps({"status": "PASS", "gates": gates,
                  "method": "add each proper divisor to its strict multiples; primality from composite marking",
                  "source_bound": bound, "target_bound": 60,
                  "producer_receipt_sha256": hashlib.sha256((base / "results.normal.json").read_bytes()).hexdigest(),
                  "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                 indent=2, sort_keys=True))
