"""Reproduce two-component birth exclusions or the order-14 complete interval."""
import argparse
from collections import defaultdict
from itertools import product
from pathlib import Path
import subprocess
import tempfile

from fibers import census_rows


def targets(shell):
    birth, by_order = {}, defaultdict(set)
    for n in range(3, 12):
        for h, _, _, _ in census_rows(n):
            if h > 1 and h not in birth:
                birth[h] = n
                by_order[n].add(h)
    if shell == 'interval14':
        out = {h for h in range(3, 3423418, 2) if h not in (7, 21) and h not in birth}
        assert (len(out), min(out), max(out)) == (1667796, 80453, 3423417)
        return out
    out = set()
    if shell == 15:
        # Delete the first vertex of each path: M(12) <= 12 M(11).
        limit = 12 * max(birth)
        out.update(3 * b for b in range(3, limit + 1, 2)
                   if b not in (7, 21) and 3 * b not in birth)
    for n in range(max(3, shell - 11), shell // 2 + 1):
        for a, b in product(by_order[n], by_order[shell - n]):
            if a * b not in birth:
                out.add(a * b)
    expected = {14: (29253, 80517, 285285), 15: (584052, 80453, 3423417)}
    assert (len(out), min(out), max(out)) == expected[shell]
    return out


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('mode', choices=('14', '15', 'interval14'))
    cli.add_argument('--witness', type=int, help='emit and independently recount one literal target')
    args = cli.parse_args()
    shell = 15 if args.mode == 'interval14' else int(args.mode)
    missing = targets(args.mode if args.mode == 'interval14' else shell)
    if args.witness is not None and args.witness not in missing:
        cli.error('the requested count is outside the residual target set')
    print(f'{len(missing)} residual target heights', flush=True)
    with tempfile.TemporaryDirectory() as tmp:
        exe = Path(tmp) / 'birth_cut_cover'
        src = Path(__file__).with_suffix('.cpp')
        subprocess.run(['c++', '-O3', '-std=c++17', str(src), '-o', str(exe)], check=True)
        requested = sorted(missing) if args.witness is None else [args.witness]
        subprocess.run([str(exe), str(shell), str(args.witness or 0)],
                       input='\n'.join(map(str, requested)), text=True, check=True)


if __name__ == '__main__':
    main()
