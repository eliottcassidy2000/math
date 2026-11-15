"""Finite base checks for the all-height opposite-fan theorem."""

from movable_census import append, F

SOURCE = (118, 124, 248, 241, 224, 192, 128, 3)
P, Q = 1, 7

def converse(r):
    return tuple(x ^ ((1 << len(r)) - 1) ^ (1 << i) for i, x in enumerate(r))

def reverse_leaf(r, p):
    return converse(append(converse(r), p))

def grid(a, b):
    r = SOURCE
    for _ in range(a):
        r = append(r, Q)
    for _ in range(b):
        r = reverse_leaf(r, P)
    return r

def main():
    expected = [[129, 259, 515], [257, 515, 1027], [513, 1027, 2051]]
    # Use explicit assertion so optimized Python cannot skip the evidence.
    rows = [[F(grid(a,b))[0] for b in range(3)] for a in range(3)]
    if rows != expected:
        raise RuntimeError(rows)
    for a in range(3):
        r = grid(a, 0)
        if F(r)[1][Q] != -1 or F(converse(r))[1][P] != 1:
            raise RuntimeError((a, r))
    print("Verified source, all nine grid values, and both boundary signs.")
    for row in rows:
        print(*row)

if __name__ == "__main__":
    main()
