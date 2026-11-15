"""Exact table for the pentagon origin and the old size-eight boundary example."""
from movable_census import F, append
from opposite_fans import reverse_leaf, SOURCE

C = (6, 12, 24, 17, 3)
P, Q = 1, 4

def point(a, b):
    r = C
    for _ in range(a):
        r = append(r, Q)
    for _ in range(b):
        r = reverse_leaf(r, P)
    return r

def main():
    table = [
        [F(point(a, b))[0] for b in range(3)]
        for a in range(3)
    ]
    if table != [[15, 33, 65], [33, 67, 131], [65, 131, 259]]:
        raise RuntimeError(table)
    # Infer the old 8-vertex source by moving q to the last label.
    old_order = (0, 1, 2, 3, 5, 6, 7, 4)
    r = point(3, 0)
    relabeled = tuple(
        sum(((r[i] >> j) & 1) << k for k, j in enumerate(old_order))
        for i in old_order
    )
    if relabeled != SOURCE:
        raise RuntimeError(relabeled)
    for row in table:
        print(*row)
    print("Size-eight source recovered as point (3,0).")

if __name__ == "__main__":
    main()
