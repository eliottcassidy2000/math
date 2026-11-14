"""Complete isomorphism-class check through 8, not a sampled tail."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent /
                      "session328_fermat_gluing_and_phase_repair_20261114"))
from fermat_profiles import rooted, canon, extend, strong, glue, ends


def main():
    expected = {1: (3, 17), 3: (6, 51), 4: (16, 149), 5: (59, 443),
                6: (313, 1617), 7: (1891, 6111), 8: (10746, 28503)}
    classes = {(0,)}
    covered = {}
    for n in range(1, 9):
        if n > 1:
            classes = {canon(extend(A, mask)) for A in classes
                       for mask in range(1 << (n - 1))}
        for A in sorted(classes):
            if not strong(A):
                continue
            for H, L, _, _, _, root in rooted(A):
                for j in (1, 2, 4):
                    covered[H + (L << j)] = A, root, j
        if n in expected and (len(covered), max(covered)) != expected[n]:
            raise RuntimeError((n, len(covered), max(covered), expected[n]))
        print(n, len(classes), len(covered), max(covered))
    for value in (35, 69, 28503):
        A, root, j = covered[value]
        if sum(ends(glue(A, root, j))[-1]) != value:
            raise RuntimeError(("literal witness", value))
    if len(classes) != 6880:
        raise RuntimeError("complete isomorphism-class count failed")


if __name__ == "__main__":
    main()
