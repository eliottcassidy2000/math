"""Exact finite certificate for Conway's classical 1/20 tail.

Run this file with Python 3 (also under -O). Optional --plot uses matplotlib.
The proof transferring the finite window to all integers is in the paired
result note. Q is tested only on the explicitly generated finite prefix.
"""
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

CHECKS = 0


def check(value, message):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(message)


def conway(limit):
    if type(limit) is not int or limit < 2:
        raise ValueError("exact integer limit >= 2 required")
    c = [0, 1, 1]
    for n in range(3, limit + 1):
        j = c[n - 1]
        k = n - j
        if not (1 <= j < n and 1 <= k < n):
            raise ValueError("illegal recursive address")
        c.append(c[j] + c[k])
    return c


def window_end(n0, q):
    if type(n0) is not int or n0 < 2 or not isinstance(q, Fraction):
        raise ValueError("integer N>=2 and an exact rational slope required")
    if not Fraction(1, 2) < q < 1:
        raise ValueError("slope outside (1/2,1)")
    threshold = (n0 - 1 - q) / (1 - q)
    ceiling = -(-threshold.numerator // threshold.denominator)
    return max(2 * n0, ceiling)


def authenticate_tail(values, n0, q):
    """Authenticate every source value, then the finite induction boundary."""
    stop = window_end(n0, q)
    if len(values) < stop or values[:3] != [0, 1, 1]:
        raise ValueError("missing prefix or wrong initial values")
    for n in range(1, stop):
        if type(values[n]) is not int:
            raise ValueError("nonintegral sequence value")
        if n >= 3:
            j = values[n - 1]
            k = n - j
            if not (1 <= j < n and 1 <= k < n):
                raise ValueError("illegal address")
            if values[n] != values[j] + values[k]:
                raise ValueError("unauthenticated recurrence value")
        if n >= n0 and q.denominator * values[n] >= q.numerator * n:
            raise ValueError("boundary violates strict barrier")
    return {"start": n0, "exclusive_end": stop,
            "slope": str(q), "checked_boundary_values": stop - n0,
            "largest_integer_barrier": max(q.denominator * values[n] -
                q.numerator * n for n in range(n0, stop))}


def expect_rejection(call, name):
    try:
        call()
    except (ValueError, TypeError):
        check(True, name)
    else:
        check(False, "accepted hostile: " + name)


def main():
    c = conway(4096)
    qseq = [0, 1, 1]
    leaves = [(0, 0), (1, 0), (0, 1)]
    first_q_down = None
    for n in range(3, 4097):
        j, k = c[n - 1], n - c[n - 1]
        check(1 <= j < n and 1 <= k < n, "Conway legality")
        check(n <= 2 * c[n] <= 2 * (n - 1), "coupled bounds")
        check(c[n] - c[n - 1] in (0, 1), "slow growth")
        one = leaves[j][0] + leaves[k][0]
        two = leaves[j][1] + leaves[k][1]
        leaves.append((one, two))
        check(one + two == c[n] and one + 2 * two == n,
              "terminal count and conserved address mass")
        check(one == 2 * c[n] - n, "error equals singleton leaves")
        a, b = n - qseq[n - 1], n - qseq[n - 2]
        check(1 <= a < n and 1 <= b < n, "finite Q legality")
        qseq.append(qseq[a] + qseq[b])
        check(2 * qseq[n] - n == (2 * qseq[a] - a) +
              (2 * qseq[b] - b) + a + b - n, "Q source defect")
        if first_q_down is None and qseq[n] < qseq[n - 1]:
            first_q_down = n

    # An independent demand-driven recursion authenticates the finite prefix.
    @lru_cache(None)
    def recursive_c(n):
        if n < 3:
            return 1
        j = recursive_c(n - 1)
        return recursive_c(j) + recursive_c(n - j)

    for n in range(1, 4097):
        check(recursive_c(n) == c[n], "independent memoized recursion")
    for power in range(1, 13):
        check(c[1 << power] == 1 << (power - 1), "dyadic endpoints")
    for power in range(1, 12):
        edge = 1 << power
        for n in range(edge, 2 * edge + 1):
            j, k = c[n - 1], n - c[n - 1]
            check(edge // 2 <= min(j, k) <= max(j, k) <= edge, "closed previous arc")
        old_peak = max(Fraction(c[n], n) for n in range(edge // 2, edge + 1))
        new_peak = max(Fraction(c[n], n) for n in range(edge, 2 * edge + 1))
        check(new_peak <= old_peak, "nonincreasing closed-arc maximum")
    # Incoming rounding-to-one claim has actual exceptions at the boundary.
    c20 = conway((1 << 21) - 1)
    edge = 1 << 20
    outside_half_open = [n for n in range(edge, 2 * edge)
                         if not (edge // 2 <= c20[n - 1] < edge and
                                 edge // 2 <= n - c20[n - 1] < edge)]
    check(outside_half_open == list(range(2097133, 2097152)), "nineteen real boundary exceptions")
    n = 3 * (1 << 18)
    lagged_selector_bill = c20[n] - 2 * c20[n // 2]
    check(lagged_selector_bill == -1430, "lagged generation selector is a separate approximation")

    slope = Fraction(11, 20)
    receipt = authenticate_tail(c, 1490, slope)
    check(receipt["exclusive_end"] == 3308, "exact tail window")
    check(receipt["largest_integer_barrier"] == -3, "strict window margin")
    violations = [n for n in range(1024, 2049) if 20 * c[n] >= 11 * n]
    check(violations == [1487, 1489], "complete arc violations")
    check(20 * c[1489] - 11 * 1489 == 1, "last violation exact excess")
    expect_rejection(lambda: authenticate_tail(c, 1489, slope), "last bad index")
    expect_rejection(lambda: authenticate_tail(c, 1025, slope), "zero endpoint is insufficient")
    expect_rejection(lambda: authenticate_tail(c[:3307], 1490, slope), "truncated boundary")
    wrong = c.copy()
    wrong[1487] -= 1
    expect_rejection(lambda: authenticate_tail(wrong, 1490, slope), "hidden corrupt ancestor")
    expect_rejection(lambda: window_end(1490, .55), "floating slope")
    expect_rejection(lambda: window_end(True, slope), "boolean index")
    expect_rejection(lambda: window_end(1490, Fraction(1, 2)), "zero-width target")

    # Additional exact barriers probe the general consumer, not just 11/20.
    other_receipts = []
    for threshold in (Fraction(2, 3), Fraction(3, 5), Fraction(4, 7)):
        last = max(n for n in range(1, 4097) if c[n] >= threshold * n)
        record = authenticate_tail(c, last + 1, threshold)
        other_receipts.append(record)
        check(record["exclusive_end"] <= len(c), "independent slope closure")

    arcs = []
    for start in (512, 1024, 2048):
        peak = max(range(start, 2 * start + 1), key=lambda n: Fraction(c[n], n))
        bad = [n for n in range(start, 2 * start + 1) if 20 * c[n] >= 11 * n]
        arcs.append({"start": start, "end": 2 * start, "peak_index": peak,
                     "peak_value": c[peak], "peak_error": str(Fraction(c[peak], peak) - Fraction(1, 2)),
                     "violation_count": len(bad), "last_violation": max(bad) if bad else None})
    check(c[20] == qseq[20] == 12 and c[1024] == 512 and qseq[1024] == 557,
          "20 and 1024 distinguish and coincide")
    check(first_q_down == 16, "Q is not slow growing")
    report = {"status": "FINITE-EXACT boundary; PROVED global consequence in paired note",
              "universe": "Conway/Q identities 1..4096; closed arcs through4096; Conway k20 boundary probe through2097151; Q finite only",
              "certificate": receipt, "last_one_twentieth_violation": 1489,
              "last_violation_value": c[1489], "arc_1024_violations": violations,
              "arcs": arcs, "other_slopes": other_receipts,
              "sequence_controls": {"C20": c[20], "Q20": qseq[20], "C1024": c[1024],
                                    "Q1024": qseq[1024], "first_Q_decrease": first_q_down},
              "incoming_controls": {"k20_half_open_exceptions": outside_half_open,
                                    "k19_at_x_three_halves_n_times_generation_difference": lagged_selector_bill},
              "prefix_sha256": hashlib.sha256(json.dumps(c[:3308], separators=(",", ":")).encode()).hexdigest(),
              "checks": CHECKS}
    print(json.dumps(report, indent=2, sort_keys=True))
    if "--plot" in sys.argv:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.55, 1]})
        for ax, lo, hi in ((axes[0], 512, 4096), (axes[1], 1478, 1500)):
            xs = list(range(lo, hi + 1))
            ys = [c[n] / n - .5 for n in xs]
            ax.plot(xs, ys, color="#155e75", lw=1.2, marker="o" if lo == 1478 else None, ms=3)
            ax.axhline(.05, color="#b45309", linestyle="--", label="1/20")
            ax.set_xlabel("index n")
            ax.grid(alpha=.18)
        axes[0].set_ylabel("C(n)/n - 1/2")
        axes[0].set_xticks([512, 1024, 2048, 3308, 4096])
        axes[0].axvspan(1490, 3307, alpha=.09, color="#15803d", label="finite proof window")
        axes[0].legend(fontsize=8, loc="lower right")
        axes[0].set_title("Three dyadic arcs and the finite certificate")
        axes[1].set_title("Last violations: 1487 and 1489")
        axes[1].set_ylim(.04855, .05028)
        for n in violations:
            axes[1].annotate(str(n), (n, c[n] / n - .5), xytext=(0, 11), textcoords="offset points", ha="center", fontsize=9)
        fig.suptitle("Conway's classical 1/20 bound: exact arithmetic + induction", fontsize=13)
        fig.tight_layout()
        dest = Path(__file__).resolve().parents[2] / "05-knowledge/results/conway_tail_barrier_20261013c.png"
        fig.savefig(dest, dpi=180)


if __name__ == "__main__":
    main()
