"""Select a unique-start receiver and suffix together for the eight-vertex switch."""
from argparse import ArgumentParser
from math import gcd
from pathlib import Path
import sys

from cm_crossings import CHARTS
from e15_pair import p, tagged_data
from replay import inverse_one_sided

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "chanpak_e31_physical_20261122"))
from check import two_slots


def first_crossing(midpoint):
    core = next(chart for chart in CHARTS if chart.discriminant == 67)
    h, _, c, _, e, _, k = tagged_data(core.lower, core.tags)
    target = midpoint - 1
    if midpoint <= 0 or midpoint % 2:
        return None
    for v in range(1, target//h + 1):
        bits = bin(v)[2:]
        untagged_end = (1 << len(bits)) - int(bits[::-1], 2)
        even = int(v % 2 == 0)
        H, Q, U = 1+p(v), two_slots(v), v+1
        base = H*k + h*(U-even) + c*(Q-H-untagged_end) + e*(H-1)
        slope = c*Q + h*U
        weight = h*H
        for parity in (0, 1):
            intercept = base + parity*h*even
            if target == intercept and parity == 0:
                return v, 0
            if (target-intercept-parity*slope) % gcd(2*weight, 2*slope):
                continue
            matches = inverse_one_sided(target, intercept, weight, slope, parity)
            if matches:
                return v, matches[0]
    return None


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("midpoint", type=int, nargs="?", default=5280**3+744)
    target = parser.parse_args().midpoint
    solution = first_crossing(target)
    if target == 5280**3+744:
        assert solution == (1634,650)
    print("midpoint", target, "receiver and suffix", solution)
