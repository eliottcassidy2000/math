"""Census of permanent unit-gap D flips with one arbitrary outside address."""

from argparse import ArgumentParser
from collections import Counter
from functools import lru_cache
from itertools import permutations
from math import gcd, isqrt
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "05-knowledge/results/chanpak_decimal_kernel_20261122"))
from fermat_receiver import (
    both_faces,
    literal_receiver,
    tagged_data,
    flipped,
    EDGE,
    D,
    join,
    TARGET,
    p,
)


@lru_cache(None)
def machine(order, first, last, cap_mask):
    tagged = {"t"}
    tagged.update(v for j, v in enumerate("ab") if cap_mask & (1 << j))

    def arrow(u, v):
        if u == "|" or v == "|":
            return True
        if u == "*":
            return v in tagged
        if v == "*":
            return u not in tagged
        return (
            (u in ("t", "u") and v == "a")
            or (u == "a" and v == "b")
            or (u == "b" and v in ("t", "u"))
        )

    def contact(v, wanted):
        return wanted is None or (v != "|" and (v in tagged) == wanted)

    width = len(order) + 1
    transitions = [[], []]
    accepting = []
    for state in range(3**width):
        digits = [state // 3**j % 3 for j in range(width)]
        valid = True
        for j, d in enumerate(digits):
            if d == 2:
                valid = False
            elif d == 0:
                if j == 0:
                    valid &= contact(order[0], first)
                elif j == width - 1:
                    valid &= contact(order[-1], last)
                else:
                    valid &= arrow(order[j - 1], order[j])
        if valid:
            accepting.append(state)
        for untagged, symbol in enumerate("tu"):
            moves = []
            for j, d in enumerate(digits):
                may_enter = arrow(order[j - 1], symbol) if j else contact(symbol, first)
                may_leave = (
                    arrow(symbol, order[j]) if j < width - 1 else contact(symbol, last)
                )
                if d or may_enter:
                    moves.append(state + 3**j * ((1 if may_leave else 2) - d))
            transitions[untagged].append(moves)
    return transitions, accepting


def count(word, markers, caps, first=None, last=None):
    total = 0
    for order in permutations(markers):
        transitions, accepting = machine(order, first, last, caps)
        mass = [0] * len(transitions[0])
        mass[0] = 1
        for bit in word:
            nxt = [0] * len(mass)
            for i, v in enumerate(mass):
                if v:
                    for j in transitions[bit][i]:
                        nxt[j] += v
            mass = nxt
        total += sum(mass[i] for i in accepting)
    return total


@lru_cache(None)
def recv(n, r, caps):
    u = r - 2 + caps.bit_count()
    word = tuple(int(v) for v in format(u, f"0{n}b"))

    def reading(markers, first=None, last=None):
        return count(word, markers, caps, first, last)

    data = (
        reading("ab"),
        reading("ab|"),
        reading("ab|", first=True),
        reading("ab|", last=False),
        reading("ab", last=True),
        reading("ab", first=False),
        reading("ab*"),
    )
    assert data[0] == 1 + (1 << n) and data[5] == r
    if n <= 5:
        rows, mark = literal_receiver(word)
        assert tagged_data(rows, (mark % (1 << n)) + (caps << n)) == data
    return data


def faces(data, target=TARGET, axis=None, metrics=None):
    h, q, c, t, e, r, k = data
    for side, (s, odd) in enumerate(((t, r), (c, e))):
        largest = max(0, target - k)
        width = (isqrt(largest * gcd(s, h)) // s).bit_length()
        if axis and width >= 15 and h < 1 << 48:
            command = [axis, s, h, k, k + odd, target, 10]
            lines = subprocess.check_output(
                list(map(str, command)), text=True
            ).splitlines()
            for line in lines:
                fields = line.split()
                if fields[:1] == ["address"]:
                    x = int(fields[1])
                    assert k + s * x + h * p(x) + (x & 1) * odd == target
                    yield side, x, True
                if fields[:1] == ["parity"] and metrics is not None:
                    metrics["native bucket candidates"] += int(
                        fields[fields.index("candidates") + 1]
                    )
            continue
        for j in (0, 1):
            M = target - k - odd * j
            if M < 0:
                continue
            g = gcd(s, h)
            if M % g:
                continue
            m = h // g
            bits = max(1, (isqrt(M * g) // s).bit_length() - 1)
            B = 1 << bits
            inverse = pow((B * s) // g, -1, m)
            for i in range(j, B, 2):
                rhs = M - s * i - h * p(i)
                if rhs < 0:
                    continue
                weight = B - 1 - int(f"{i:0{bits}b}"[::-1], 2)
                slope = B * s + h * weight
                lo = (rhs + slope + h - 1) // (slope + h)
                hi = rhs // slope
                res = (rhs // g) * inverse % m
                first = lo + (res - lo) % m
                for X in range(first, hi + 1, m):
                    a = B * X + i
                    yield side, a, s * a + h * p(a) == M


def sources():
    left = {}
    for o in range(1, 6):
        for line in (
            (ROOT / f"session212/cycle-levels/classes{o}.rows").read_text().splitlines()
        ):
            rows = tuple(map(int, line.split()))
            for tag in range(1 << o):
                left[tagged_data(rows, tag)] = rows, tag
    mids = []
    for m in range(1 << len(D)):
        a, b = tagged_data(flipped(D, *EDGE), m), tagged_data(D, m)
        if a[:4] == b[:4] and a != b:
            mids.append((m, a, b))
    return left, mids


def controls():
    charts = ((21, 19, 12, 8, 5, 2, 13), (5, 11, 17, 25, 0, 3, 7))
    for chart in charts:
        h, _, c, t, e, r, k = chart
        for planted in range(50):
            for target in (
                k + t * planted + r * (planted & 1) + h * p(planted),
                k + c * planted + e * (planted & 1) + h * p(planted) + 1,
            ):
                expected = {
                    (side, x)
                    for side, (s, d) in enumerate(((t, r), (c, e)))
                    for x in range((target - k) // s + 1)
                    if k + s * x + d * (x & 1) + h * p(x) == target
                }
                assert expected == {
                    (side, x) for side, x, ok in faces(chart, target) if ok
                }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--axis", help="optional compiled farey_axis for the largest suffix partitions"
    )
    args = parser.parse_args()
    controls()
    for caps in range(4):
        for z in (0, 1, 3, 6, 11):
            recv(4, z + 2 - caps.bit_count(), caps)
    left, mids = sources()
    print(len(left), len(mids), flush=True)
    signatures = {}
    for A, (arows, atag) in left.items():
        for m, D0, D1 in mids:
            n0, n1 = join(A, D0), join(A, D1)
            signatures[n0, n1] = (arows, atag, m)
    print("prefix pairs", len(signatures), flush=True)
    num = Counter()
    hit = []
    for n in range(24, 121):
        G = (1 << n) + 1
        active = {}
        for (v, w), names in signatures.items():
            c0 = G * (w[-1] - v[-1])
            slope = w[4] - v[4]
            if not slope:
                continue
            for sign in (2, -2):
                r, rem = divmod(sign - c0, slope)
                if rem:
                    continue
                for caps in range(4):
                    if 0 <= r - 2 + caps.bit_count() < (1 << n):
                        active.setdefault((r, caps), []).append((v, w, names))
        local = Counter()
        for (r, caps), items in active.items():
            rec = recv(n, r, caps)
            for v, w, names in items:
                a, b = join(v, rec), join(w, rec)
                assert a[:-1] == b[:-1] and abs(a[-1] - b[-1]) == 2
                data = min(a, b)
                local["charts"] += 1
                if data[-1] > TARGET:
                    continue
                local["live"] += 1
                _, q, c, t, _, _, k = data
                if n >= 52 and k + q + c + t <= TARGET:
                    local["positive-pair searches"] += 1
                    pairs = both_faces(data, TARGET, limit=30000)
                    assert pairs is not None, (n, r, caps, data)
                    for x, y in pairs:
                        print("TWO-FLANK HIT", n, r, caps, names, x, y, flush=True)
                        hit.append((n, r, caps, names, x, y))
                for side, x, ok in faces(data, axis=args.axis, metrics=local):
                    local["literal evaluations and native hits"] += 1
                    if ok:
                        print("HIT", n, r, caps, names, side, x, flush=True)
                        hit.append((n, r, caps, names, side, x))
        if active:
            print(n, len(active), dict(local), flush=True)
        local["receivers"] += len(active)
        num.update(local)
    print("TOTAL", dict(num), "hits", hit)


if __name__ == "__main__":
    main()
