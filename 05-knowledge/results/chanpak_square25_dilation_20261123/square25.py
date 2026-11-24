"""Exact strong heights on the two sheets of the square-sum graph through 25."""
from collections import Counter
from functools import lru_cache
from math import isqrt

# (odd part of the base label, offset): (slope, contact mask, outgoing rows)
RAILS = {
    (1, -3): (128, 40, (2, 16, 11, 19, 69, 31, 47)),
    (1, -1): (1024, 104, (32, 9, 11, 81, 71, 30, 167, 63)),
    (3, -3): (768, 10, (350, 508, 504, 496, 480, 193, 384, 257, 32)),
    (3, -1): (6, 3, (0, 1, 3)),
    (5, -3): (10240, 42, (318, 108, 472, 176, 482, 452, 393, 259, 10)),
    (5, -1): (5, 5, (0, 1, 3)),
    (7, -3): (448, 38, (32, 5, 65, 7, 15, 30, 59)),
    (7, -1): (224, 160, (2, 4, 17, 71, 75, 31, 167, 63)),
    (9, -3): (9, 7, (2, 4, 1)),
    (9, -1): (9, 9, (2, 4, 1, 7)),
    (11, -3): (88, 9, (2, 4, 33, 7, 15, 27)),
    (11, -1): (176, 20, (2, 64, 3, 7, 15, 31, 61)),
    (13, -3): (1664, 37, (6, 36, 192, 7, 15, 29, 59, 123)),
    (13, -1): (26, 12, (2, 16, 3, 7, 13)),
    (15, -3): (15, 17, (2, 4, 1, 7, 15)),
    (15, -1): (15, 3, (0, 5, 9, 3)),
    (17, -3): (17, 20, (2, 8, 3, 5, 15)),
    (17, -1): (34, 9, (2, 4, 17, 7, 11)),
    (19, -3): (19, 18, (2, 8, 3, 5, 15)),
    (19, -1): (608, 136, (4, 65, 18, 7, 43, 143, 61, 95)),
    (21, -3): (336, 9, (2, 64, 67, 23, 39, 15, 57)),
    (21, -1): (336, 28, (8, 5, 33, 70, 15, 27, 55)),
    (23, -3): (368, 56, (16, 9, 35, 37, 14, 83, 31)),
    (23, -1): (46, 24, (2, 32, 3, 7, 15, 29)),
    (25, -3): (200, 34, (16, 5, 9, 19, 70, 31, 47)),
    (25, -1): (1600, 82, (6, 36, 80, 71, 75, 29, 163, 63)),
    (1, -5): (2048, 42, (10, 36, 41, 194, 15, 153, 55, 87)),
    (1, 1): (4, 1, (0, 1, 3)),
    (3, -5): (96, 3, (0, 1, 11, 35, 15, 23)),
    (3, 1): (6, 3, (2, 4, 1)),
    (5, -5): (160, 5, (0, 1, 3, 23, 71, 31, 47)),
    (5, 1): (20, 24, (2, 4, 17, 7, 11)),
    (7, -5): (7, 7, (0, 1, 3)),
    (7, 1): (112, 24, (2, 32, 3, 23, 39, 77, 31)),
    (9, -5): (18, 14, (2, 8, 3, 5)),
    (9, 1): (1152, 100, (6, 36, 136, 35, 15, 85, 31, 123)),
    (11, -5): (88, 48, (2, 4, 9, 67, 15, 31, 55)),
    (11, 1): (11, 12, (2, 8, 3, 5)),
    (13, -5): (416, 50, (32, 13, 25, 49, 67, 22, 47)),
    (13, 1): (208, 48, (2, 64, 3, 7, 47, 79, 157, 63)),
    (15, -5): (480, 3, (0, 5, 9, 3, 47, 143, 63, 95)),
    (15, 1): (60, 40, (16, 5, 17, 7, 42, 15)),
    (17, -5): (2176, 140, (144, 21, 41, 67, 140, 27, 55, 110)),
    (17, 1): (34, 48, (2, 8, 3, 37, 15, 23)),
    (19, -5): (38, 9, (0, 5, 17, 7, 11)),
    (19, 1): (2432, 80, (374, 460, 504, 497, 482, 450, 384, 1, 128)),
    (21, -5): (21, 11, (2, 8, 3, 5)),
    (21, 1): (168, 66, (16, 5, 17, 7, 74, 31, 47)),
    (23, -5): (23, 17, (2, 8, 3, 5, 15)),
    (23, 1): (368, 44, (32, 5, 17, 7, 43, 78, 31)),
    (25, -5): (25, 40, (2, 8, 3, 5, 15, 31)),
    (25, 1): (100, 10, (4, 5, 48, 7, 11, 27)),
}

# Strong initial witnesses before the one-return width becomes positive.
INITIAL = {
    13: (22, 12, 24, 17, 2),
    15: (6, 12, 24, 17, 3),
    29: (38, 60, 24, 49, 33, 4),
    31: (46, 20, 56, 50, 33, 2),
    45: (38, 28, 24, 49, 33, 6),
    61: (118, 76, 120, 113, 34, 66, 16),
    63: (246, 252, 248, 49, 224, 192, 136, 8),
    77: (94, 108, 120, 48, 98, 65, 8),
    93: (62, 92, 120, 112, 96, 66, 1),
    109: (78, 60, 120, 112, 97, 65, 2),
    111: (78, 28, 120, 112, 97, 67, 2),
    125: (42, 92, 89, 112, 33, 70, 17),
    127: (106, 60, 89, 48, 97, 68, 10),
    157: (78, 28, 56, 112, 97, 67, 6),
    175: (14, 28, 56, 112, 97, 67, 7),
    189: (22, 44, 88, 49, 98, 69, 11),
    205: (222, 60, 248, 240, 224, 193, 130, 2),
    221: (158, 188, 248, 240, 96, 193, 131, 16),
    223: (222, 188, 184, 240, 224, 65, 134, 32),
    255: (190, 156, 248, 112, 224, 194, 131, 8),
    303: (190, 28, 248, 240, 96, 194, 131, 18),
    317: (222, 124, 184, 240, 224, 193, 132, 2),
    333: (158, 188, 248, 112, 224, 65, 131, 40),
    335: (190, 220, 120, 240, 96, 194, 129, 20),
    365: (190, 220, 56, 112, 224, 194, 133, 12),
    381: (158, 124, 248, 176, 224, 65, 137, 34),
    399: (126, 188, 152, 240, 224, 196, 6, 65),
    445: (158, 124, 184, 240, 224, 193, 5, 66),
    511: (158, 124, 56, 240, 224, 193, 133, 6),
    637: (170, 92, 153, 112, 225, 70, 133, 42),
    765: (446, 380, 376, 496, 480, 192, 385, 6, 160),
    829: (286, 508, 376, 496, 224, 449, 129, 261, 80),
    1023: (414, 508, 248, 112, 480, 321, 385, 296, 12),
    1277: (126, 444, 504, 368, 224, 448, 386, 265, 17),
    1599: (350, 444, 56, 496, 224, 449, 390, 261, 20),
    2557: (62, 60, 472, 496, 224, 452, 387, 259, 19),
    5117: (638, 216, 658, 756, 224, 710, 388, 513, 703, 82),
    10237: (418, 984, 75, 625, 485, 390, 801, 76, 140, 437),
}


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def path_engine(rows):
    n = len(rows)
    full = (1 << n) - 1
    assert all(0 <= row <= full for row in rows)
    assert all(((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) == (i != j)
               for i in range(n) for j in range(n))
    incoming = [full ^ row ^ (1 << i) for i,row in enumerate(rows)]

    @lru_cache(None)
    def ending(subset, end):
        rest = subset ^ (1 << end)
        return sum(ending(rest, i) for i in bits(rest & incoming[end])) if rest else 1

    @lru_cache(None)
    def starting(subset, start):
        rest = subset ^ (1 << start)
        return sum(starting(rest, i) for i in bits(rest & rows[start])) if rest else 1

    @lru_cache(None)
    def total(subset):
        return sum(ending(subset, i) for i in bits(subset)) if subset else 1

    return full, ending, starting, total


def reaches(rows, contact):
    for _ in rows:
        for i, row in enumerate(rows):
            if row & contact:
                contact |= 1 << i
    return contact == (1 << len(rows)) - 1


def compile_rows(core, contact, width):
    n = len(core)
    assert width > 0
    module_mask = ((1 << width) - 1) << n
    full_core = (1 << n) - 1
    apex = n + width
    fixed = tuple(row | (((contact >> i) & 1) << apex) for i, row in enumerate(core))
    module = tuple((module_mask & ~((1 << (n+j+1)) - 1)) | full_core for j in range(width))
    return fixed + module + (module_mask | (full_core ^ contact),)


def profile(rows, contact):
    full, end, start, total = path_engine(rows)
    gain = correction = 0
    for subset in range(1, full+1):
        rest = full ^ subset
        left = sum(end(subset,i) for i in bits(subset & contact))
        gain += left * total(rest)
        correction += left * sum(start(rest,j) for j in bits(rest & contact))
    return gain, total(full) - correction


def witness(k, v, offset):
    assert k >= (6 if offset in (-5,1) else 2) and 1 <= v <= 25 and offset in (-5,-3,-1,1)
    center = v * (4 ** k)
    oddpart = v // (v & -v)
    gain, contact, core = RAILS[oddpart,offset]
    if center < 2 * gain:
        return INITIAL[center + offset]
    width = (center // gain).bit_length() - 1
    assert gain * (1 << width) == center
    return compile_rows(core, contact, width)


# On the two orientations of one marked arc the leading slope agrees.
FLIPS = (
    (9, 9, (4,1,2,7), (0,2)),
    (15, 17, (10,4,1,6,15), (0,3)),
    (736, 130, (4,5,24,131,107,15,47,119), (4,5)),
)


def square(x):
    return isqrt(x)**2 == x


def board(n, removed=()):
    verts = {(v,s) for v in range(1,n+1) for s in (-1,1)} - set(removed)
    return {(v,s): {(w,-s) for w in range(1,n+1) if (w,-s) in verts and square(v+w)}
            for v,s in verts}


def cycles(adj):
    # Root and orientation select each undirected cycle exactly once. Each
    # unused vertex needs two available neighbors among unused vertices and ends.
    root = min(adj, key=lambda v: (len(adj[v]),v))

    def visit(path, unused):
        if not unused:
            if path[-1] in adj[root] and path[1] < path[-1]:
                yield tuple(path)
            return
        for nxt in adj[path[-1]] & unused:
            rest = unused - {nxt}
            available = rest | {nxt,root}
            if all(len(adj[x] & available) >= 2 for x in rest):
                yield from visit(path+[nxt],rest)

    return list(visit([root],set(adj)-{root}))


def edges(cycle):
    return {frozenset((a,b)) for a,b in zip(cycle,cycle[1:]+cycle[:1])}


def verify():
    for (a,b),(gain,contact,rows) in RAILS.items():
        assert gain // (gain & -gain) == a and reaches(rows,contact)
        assert profile(rows,contact) == (gain,b)
        host = compile_rows(rows,contact,2)
        f,_,_,total = path_engine(host)
        assert total(f) == 4 * gain + b
    for gain, contact, rows, (u,v) in FLIPS:
        changed = list(rows)
        changed[u] ^= 1 << v
        changed[v] ^= 1 << u
        assert reaches(rows, contact) and reaches(changed, contact)
        assert {profile(rows, contact), profile(tuple(changed), contact)} == {
            (gain,-1), (gain,-3)}
    for h,rows in INITIAL.items():
        f,_,_,total = path_engine(rows)
        assert total(f) == h and all(reaches(rows,1<<i) for i in range(len(rows)))
    for k in range(2,8):
        for v in range(1,26):
            for b in (-3,-1):
                assert len(witness(k,v,b)) <= 2*k+14
    for v in range(1,26):
        for b in (-5,1):
            assert len(witness(6,v,b)) <= 26
    print("52 affine strong profiles and",len(INITIAL),"strong initial witnesses")
    for n,expected in [(23,1),(24,0),(25,2),(26,4)]:
        cs = cycles(board(n))
        assert len(cs) == expected
        print("labels per sheet",n,"cycles",len(cs))
        if n == 25:
            inventory = Counter(tuple(sorted(v for edge in edges(c) for v in range(1,26)
                                              if edge=={(v,1),(v,-1)})) for c in cs)
            assert inventory == {(2,18):1,(8,18):1}
            for c in cs:
                e = edges(c)
                for left,right in ((16,9),(13,12),(19,6)):
                    for sign in (-1,1):
                        assert {(left,sign),(right,-sign)} in e
    for u,v,U,V in ((16,9,13,12),(16,9,19,6)):
        assert all(square(x) for x in (u+v,U+V,u+64*V,v+64*U))
    for k in (6,10):
        heights = [4**k*v+b for v in range(1,26) for b in (-5,-3,-1,1)]
        assert len(set(heights)) == 100 and {x%24 for x in heights} == set(range(1,24,2))
        primary = [h for h in heights if h%8 in (5,7)]
        secondary = [h for h in heights if h%8 in (1,3)]
        assert all(not square(h+h2+4) for h in primary for h2 in secondary)
    deleted = board(25, [(2,1),(6,-1)])
    leaf_heights = sorted(4*v-2+s for (v,s),nb in deleted.items() if len(nb)==1)
    assert leaf_heights == [39,75,89]
    print("three leaves at first scale",leaf_heights)


if __name__ == "__main__":
    verify()
