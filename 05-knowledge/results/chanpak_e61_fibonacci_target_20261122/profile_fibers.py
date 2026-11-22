"""Count classes and full weak-slot profiles in the two E61 cyclic grammars."""
import argparse
import ast
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

PACKET = Path(__file__).resolve().parent
SOURCE = PACKET.parent / "chanpak_e61_compact_volchkov_20261122" / "check_strong_93.py"


def family_inputs():
    assignments = {}
    for item in ast.parse(SOURCE.read_text()).body:
        if isinstance(item, ast.Assign) and isinstance(item.targets[0], ast.Name):
            name = item.targets[0].id
            if name in {"MODULES", "B91"}:
                assignments[name] = ast.literal_eval(item.value)
    return {
        91: [(rows, h, odd << valuation) for h, odd, valuation, rows in assignments["B91"]],
        93: assignments["MODULES"],
    }


def converse_classes(rows):
    """Exact isomorphism test against arc reversal, restricted by out-degree."""
    n = len(rows)
    mask = (1 << n) - 1
    reverse = [m ^ mask ^ (1 << i) for i, m in enumerate(rows)]
    options = [[j for j in range(n) if reverse[j].bit_count() == m.bit_count()] for m in rows]
    order = sorted(range(n), key=lambda i: len(options[i]))
    image = {}

    def search(index, used):
        if index == n:
            return True
        i = order[index]
        for j in options[i]:
            if used >> j & 1:
                continue
            if any((rows[i] >> u ^ reverse[j] >> v) & 1 for u, v in image.items()):
                continue
            image[i] = j
            if search(index + 1, used | 1 << j):
                return True
            del image[i]
        return False

    return 1 if search(0, 0) else 2


def full_profile(rows):
    """Coefficients of (x)_j in the weak ordered-slot counting polynomial."""
    n = len(rows)
    ends = [[0] * n for _ in range(1 << n)]
    h = [1]
    covers = [[int(s == 0 and j == 0) for j in range(n + 1)] for s in range(1 << n)]
    for s in range(1, 1 << n):
        for i in range(n):
            if s >> i & 1:
                rest = s ^ 1 << i
                ends[s][i] = int(not rest) + sum(ends[rest][u] for u in range(n) if rows[u] >> i & 1)
        h.append(sum(ends[s]))
        anchor = s & -s
        remaining = s ^ anchor
        sub = remaining
        while True:
            chosen = sub | anchor
            rest = s ^ chosen
            for j in range(rest.bit_count() + 1):
                covers[s][j + 1] += h[chosen] * covers[rest][j]
            if not sub:
                break
            sub = (sub - 1) & remaining
    return tuple(covers[-1][1:])


def census(worker, targets):
    profiles = defaultdict(Counter)
    for n in sorted({s[0] for s in targets if s[0] < 10}):
        pairs = sorted({(h, q) for order, h, q in targets if order == n})
        with subprocess.Popen(["gentourng", "-cq", str(n)], stdout=subprocess.PIPE) as gen:
            args = [str(worker), str(n)] + [str(x) for pair in pairs for x in pair]
            data = subprocess.check_output(args, stdin=gen.stdout, text=True)
            gen.stdout.close()
            if gen.wait():
                raise RuntimeError("gentourng failed")
        for line in data.splitlines():
            if line.startswith("inputs ") or line.startswith("sample "):
                print(f"order {n}: {line}")
                continue
            h, q, count, *profile = map(int, line.split())
            profiles[n, h, q][tuple(profile)] = count
    return profiles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("worker", type=Path, help="compiled profile_fibers.cpp")
    families = family_inputs()
    targets = {(len(rows), h, q): rows for family in families.values() for rows, h, q in family}
    profiles = census(parser.parse_args().worker.resolve(), targets)
    choices, full_matches = {}, {}
    for key, rows in sorted(targets.items()):
        n, h, q = key
        conv = converse_classes(rows)
        profile = full_profile(rows)
        if profile[0] != h or 2 * (profile[0] + profile[1]) != q:
            raise ValueError("source profile mismatch")
        if n == 10:
            profiles[key][profile] = conv
        choices[key] = sum(profiles[key].values())
        full_matches[key] = profiles[key][profile]
        assert choices[key] >= full_matches[key] >= conv
        print(key, "H,Q choices", choices[key], "same full profile", full_matches[key], "converse choices", conv)
    for order, family in families.items():
        singletons = {91: 35, 93: 24}[order]
        mult = Counter((len(rows), h, q) for rows, h, q in family)
        permutations = math.factorial(len(family) + singletons) // math.factorial(singletons)
        permutations //= math.prod(math.factorial(m) for m in mult.values())
        for desc, table in (("H,Q", choices), ("full profile of originals", full_matches)):
            multiplier = math.prod(table[key] ** m for key, m in mult.items())
            print(order, desc, "multiplier", multiplier, "strong classes", permutations * multiplier)


if __name__ == "__main__":
    main()
