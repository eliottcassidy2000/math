#!/usr/bin/env python3
"""Primary exact certificate for proposed THM-4374.

The script rebuilds THM-4365's centered-residue metric exit and exhausts the
THM-4367 active scale fibres under P -> P+2.  It checks the sharp horizon
filtration, the global 17-step phase cover, and the horizon-16 hostile.  It
imports no repository computation and uses an explicit check function that
remains active under ``python -O``.
"""

from collections import Counter
from fractions import Fraction
from math import gcd
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


A = 3371
S = 1303
M = 14 * A
N = M // 2
TAIL = 11019
UNITS14 = (1, 3, 5, 9, 11, 13)

CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError("check failed: " + message)


def centered_rho(P: int) -> int:
    remainder = (S * P) % M
    if remainder > M // 2:
        remainder -= M
    return remainder


def metric_exit(P: int) -> Fraction:
    rho_value = centered_rho(P)
    if abs(rho_value) < A:
        c = A - rho_value
        return Fraction(S * P + c, M * P)
    return Fraction(S, M)


def phase_is_active(x: int) -> bool:
    """For x=c/2 modulo N, strict activity is 1,...,A-1."""
    return 1 <= x % N <= A - 1


def first_hit(x: int) -> int:
    for time in range(17):
        if phase_is_active(x - S * time):
            return time
    raise RuntimeError(f"phase {x} has no active hit through time 16")


def representative_class(a: int, g0: int) -> tuple[int, int]:
    """Choose one large admissible (b,kappa) for a scale residue."""
    k0 = pow(g0, -1, 14)
    b = (pow(S, -1, M) * (A * k0 - a)) % M
    if b == 0:
        b = M
    while b < TAIL or gcd(a, b) != 1:
        b += M
    numerator = a + S * b
    require(numerator % A == 0, "normal-form divisibility")
    kappa = numerator // A
    require(b % 2 == 1 and gcd(a, b) == 1, "primitive parity")
    require((g0 * kappa) % 14 == 1, "scale residue")
    return b, kappa


require(M == 47194 and N == 23597 and N == 7 * A, "constants")
require(gcd(S, M) == 1, "rotation is primitive")

# Exact first-hit clock on the full half-phase circle.
hit_counts = Counter()
for x in range(N):
    hit = first_hit(x)
    hit_counts[hit] += 1
    require(0 <= hit <= 16, "phase hitting bound")
    require(phase_is_active(x - S * hit), "reported hit is active")
    for time in range(hit):
        require(not phase_is_active(x - S * time), "reported hit is first")

expected_hits = {0: A - 1, 16: 682}
expected_hits.update({time: S for time in range(1, 16)})
require(dict(sorted(hit_counts.items())) == expected_hits,
        "exact hit distribution")

# Every positive-time first entry is followed by a second active output.
for x in range(N):
    hit = first_hit(x)
    if hit:
        entry = (x - S * hit) % N
        require(2068 <= entry <= 3370, "first-entry half-coordinate")
        require(phase_is_active(entry - S), "first entry followed by activity")

# Exhaust every possible pair of scale addresses in a current active metric
# fibre.  Each structural (a,g0) class gets a literal admissible (b,kappa).
pair_count = 0
first_split_counts = Counter()
max_classes = Counter()

for a in range(2, 6741, 2):
    scale_cap = 6740 // a
    for g0 in UNITS14:
        if g0 > scale_cap:
            continue
        scales = tuple(range(g0, scale_cap + 1, 14))
        if not scales:
            continue
        b, kappa = representative_class(a, g0)
        words: dict[int, tuple[Fraction, ...]] = {}
        for scale in scales:
            P = b * scale
            require(P >= TAIL and P % 2 == 1, "tail point")
            require(centered_rho(P) == A - a * scale,
                    "literal centered cell")
            require(metric_exit(P) == Fraction(kappa, 14 * b),
                    "common metric exit")
            words[scale] = tuple(metric_exit(P + 2 * time)
                                 for time in range(18))

        for horizon in range(18):
            buckets = Counter(word[:horizon + 1] for word in words.values())
            max_classes[horizon] = max(max_classes[horizon],
                                       max(buckets.values()))

        for index, scale1 in enumerate(scales):
            for scale2 in scales[index + 1:]:
                pair_count += 1
                word1, word2 = words[scale1], words[scale2]
                require(word1[0] == word2[0], "same current fibre")
                split = next(time for time in range(1, 18)
                             if word1[time] != word2[time])
                c1, c2 = a * scale1, a * scale2
                if c2 >= 2608:
                    predicted = 1
                elif c1 >= 1244:
                    predicted = 17
                else:
                    predicted = 16
                require(split == predicted, "exact first-split clock")
                first_split_counts[split] += 1

require(pair_count == 281073, "pair universe size")
require(set(first_split_counts) == {1, 16, 17}, "only three split times")
require(max_classes[0] == 241, "sharp current metric multiplicity")
for horizon in range(1, 16):
    require(max_classes[horizon] == 94,
            f"sharp horizon-{horizon} multiplicity")
require(max_classes[16] == 49, "sharp horizon-16 multiplicity")
require(max_classes[17] == 1, "horizon-17 active-fibre injectivity")

# One class simultaneously attains all three nontrivial maxima.
a = 2
b = 47595
kappa = 18397
scales = tuple(1 + 14 * index for index in range(241))
require(a + S * b == A * kappa, "sharp class equation")
require(all((scale * kappa) % 14 == 1 for scale in scales),
        "sharp class residues")
require(sum(a * scale <= 2606 for scale in scales) == 94,
        "sharp low bucket")
require(sum(1244 <= a * scale <= 2606 for scale in scales) == 49,
        "sharp time-16 bucket")

# A literal pair proves that horizon 16 is insufficient.
scale1, scale2 = 631, 645
P1, P2 = 401 * scale1, 401 * scale2
word1 = tuple(metric_exit(P1 + 2 * time) for time in range(18))
word2 = tuple(metric_exit(P2 + 2 * time) for time in range(18))
require(P1 >= TAIL and P2 >= TAIL, "sharp pair in tail")
require(word1[:17] == word2[:17], "horizon 16 is not injective")
require(word1[17] != word2[17], "sharp pair splits at 17")
require(word1[0] == Fraction(155, 5614), "sharp current exit")
require(all(value == Fraction(S, M) for value in word1[1:17]),
        "sharp baseline block 1")
require(all(value == Fraction(S, M) for value in word2[1:17]),
        "sharp baseline block 2")
require(word1[17] == Fraction(97819, 3542910), "sharp first return 1")
require(word2[17] == Fraction(99989, 3621506), "sharp first return 2")

# Literal multi-period smoke census for the symbolic global proof.
seen: dict[tuple[Fraction, ...], int] = {}
sample_start = TAIL
sample_stop = TAIL + 2 * M
for P in range(sample_start, sample_stop, 2):
    word = tuple(metric_exit(P + 2 * time) for time in range(18))
    require(word not in seen, "literal horizon-17 collision")
    seen[word] = P
require(len(seen) == M, "two full odd-residue periods sampled")

print("THM-4374 seventeen-step metric-exit observability primary: PASS")
print(f"checks={CHECKS}")
print(f"phases={N}")
print("first_hit=" + ",".join(f"{time}:{hit_counts[time]}"
                              for time in range(17)))
print(f"active_fibre_scale_pairs={pair_count}")
print("first_split=" + ",".join(
    f"{time}:{first_split_counts[time]}" for time in sorted(first_split_counts)
))
print("max_horizon_fibre=" + ",".join(
    f"{horizon}:{max_classes[horizon]}" for horizon in range(18)
))
print(f"sharp_pair={P1},{P2}")
print(f"sharp_word_0={word1[0]}")
print(f"sharp_word_1_to_16={Fraction(S, M)}")
print(f"sharp_word_17={word1[17]},{word2[17]}")
print(f"global_sample=[{sample_start},{sample_stop})_odd,size={len(seen)}")
