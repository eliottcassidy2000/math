"""Emit certified Farey charts on stdout, or test a compiled address inverse."""
import argparse
import subprocess
from pathlib import Path

from fermat_receiver import A, D, EDGE, TARGET, address_height, flipped, join, receiver, tagged_data


def chart(n):
    return join(join(tagged_data(A, 12), tagged_data(flipped(D, *EDGE), 35)), receiver(n))


def check(binary):
    for n in (4, 12, 24, 32):
        data = chart(n)
        row = " ".join(map(str, (n, *data))) + "\n"
        for x0, y0 in ((0, 0), (0, 13), (8, 0), (1, 1), (7, 11), (47, 3)):
            target = address_height(data, x0, y0)
            h, q, c, t, e, r, k = data
            expected = set()
            for x in range(1+(target-k)//t):
                for y in range(1+(target-k-t*x)//(c+q*x)):
                    if address_height(data, x, y) == target:
                        expected.add((x, y))
            raw = subprocess.check_output([str(binary),str(n),str(target)],input=row,text=True)
            actual = [(int(a),int(b)) for tag,a,b in
                      (line.split() for line in raw.splitlines() if line.startswith("HIT "))]
            assert set(actual) == expected and len(actual) == len(set(actual)), (n,x0,y0,actual,expected)
        # Exercise the long quotient and tangent inversion with planted inner, balanced
        # and outer addresses otherwise too large for a full Cartesian control.
        period = 290 if n == 12 else 10
        for x, y in ((period*10_000_001+3, 13_000_003), (period*13_000_003+7, 10_000_001),
                     (31, (1<<37)+123), ((1<<39)+777, 3), (75003, 99991), (740791,736852)):
            target = address_height(data, x, y)
            out = subprocess.check_output([str(binary),str(n),str(target)],input=row,text=True)
            assert f"HIT {x} {y}\n" in out, (n,x,y)
    print("Full small inverse controls and long planted address checks passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path, help="compiled farey_hyperbola binary")
    args = parser.parse_args()
    if args.check:
        check(args.check.resolve())
    else:
        for n in range(4, 33, 4):
            print(n, *chart(n))
