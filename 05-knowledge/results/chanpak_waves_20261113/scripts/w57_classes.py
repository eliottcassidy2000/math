"""Wave 57.  Isomorphism classes of tournaments of order 1..7 (A000568: 1, 1, 2, 4, 12, 56, 456) by orbit marking.
Writes w57_reps.txt: one line 'n out_0 ... out_{n-1}' per class (out_i = out-neighbour bitmask)."""
import itertools
def classes(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]; m = len(pairs)
    idx = {pr: t for t, pr in enumerate(pairs)}; perms = list(itertools.permutations(range(n)))
    seen = bytearray(1 << m)
    def relabel(bits, p):
        b = 0
        for t, (i, j) in enumerate(pairs):
            x = bits >> t & 1; a, c = p[i], p[j]
            if a < c: b |= x << idx[(a, c)]
            else: b |= (1 - x) << idx[(c, a)]
        return b
    for bits in range(1 << m):
        if seen[bits]: continue
        for p in perms: seen[relabel(bits, p)] = 1
        out = [0] * n
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out
if __name__ == '__main__':
    counts = []
    with open('w57_reps.txt', 'w') as fo:
        for n in range(1, 8):
            k = 0
            for out in classes(n):
                fo.write(f"{n} " + " ".join(map(str, out)) + "\n"); k += 1
            counts.append(k)
    print("classes by order 1..7:", counts)
