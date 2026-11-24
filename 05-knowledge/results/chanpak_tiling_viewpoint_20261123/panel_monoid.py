"""Integer panel matrices and the blue/black binary families."""
from itertools import product
from line_targets import partner_height


def transpose(a):
    return tuple(zip(*a))


def act(v, a):
    return tuple(sum(x*y for x, y in zip(v, col)) for col in transpose(a))


def mul(a, b):
    return tuple(act(v, b) for v in a)


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def power(a, n):
    result = ((1,0,0),(0,1,0),(0,0,1))
    while n:
        if n & 1:
            result = mul(result, a)
        a = mul(a, a)
        n //= 2
    return result


# Row coordinates for the first three singleton tower states. The direction
# in these matrices is left-to-right along the anchored source.
S = ((0,1,0),(0,0,1),(1,1,1))
A = ((-1,5,3),(0,8,5),(0,17,8))
B = ((-1,4,4),(-1,7,7),(-2,14,11))
D = ((4,14,-3),(8,25,-6),(15,48,-12))
Z = ((1,1,1),(1,0,2),(3,1,3))
K = ((1,3,5),(3,5,9),(5,9,17))
EMPTY = (-1,-1,1)


def boundary(z, other):
    return z[0]*other[0] - z[1]*other[1] + z[2]*other[2]


def main():
    panels = [(((0,3),), A),
              (((0,2),(0,3)), B),
              (((0,3),(1,3)), D)]
    for arcs, matrix in panels:
        piece = ((5, (4, arcs)),)
        measurements = tuple(tuple(partner_height(piece, (i,j)) for j in (1,2,3))
                             for i in (1,2,3))
        assert mul(matrix, K) == measurements
    assert act(EMPTY, Z) == (1,0,0)
    zero, one = mul(A, B), mul(A, D)
    p, q = (2,1,1), (1,2,0)
    for start in (p, q):
        for matrix, expected in ((zero, p), (one, q)):
            raw = act(start, matrix)
            mod = tuple(x % 3 for x in raw)
            assert mod in (expected, tuple(-x % 3 for x in expected))

    for k in (1, 3, 8):
        r = 112*(k+1)
        rows = [act(EMPTY, zero)]
        for _ in range(k):
            rows = [act(v, m) for v in rows for m in (zero, one)]
        assert len(set(rows)) == 2**k
        z = [act(v, Z) for v in rows]
        after = [act(act(v, power(S, r)), Z) for v in rows]
        blue = [boundary(a,b) for a,b in zip(after, z)]
        open_word = [boundary(a, z[0]) for a in after]
        assert len(set(blue)) == len(set(open_word)) == 2**k
        assert all(h>0 and h%2 for h in blue+open_word)
        print("bits",k,"order",128*(k+1),"distinct blue outputs",len(set(blue)),
              "distinct black outputs tested",len(set(open_word[1:])))
    # A small direct orientation replay of an unequal left/right panel pair.
    left, right = panels[0], panels[1]
    word = ((5, (4, left[0])), (5, (4, right[0])))
    v = act(act(EMPTY, left[1]), right[1])
    assert dot(v, (1,1,3)) == partner_height(word, (0,0,0))


if __name__ == "__main__":
    main()
