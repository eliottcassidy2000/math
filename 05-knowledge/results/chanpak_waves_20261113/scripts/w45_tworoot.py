"""Two-root appends.  Letter i appends a vertex z beaten by every old vertex except root x_i, which z beats
(the rooted append of inverse_fences, charge a = 2r-(h-1) conserved at that root).  Mirror letter i' : z beats every
old vertex except x_i, which beats z.  Binary odometer test: H(w i) = 2 H(w) - 1 + a_i for all words w."""
import itertools, random, sys
from w45_hk import H_of
def append(n, out, x, mirror):
    o = out[:] + [0]; z = n
    for v in range(n):
        if (v == x) != mirror:   # z beats v
            o[z] |= 1 << v
        else:
            o[v] |= 1 << z
    return n + 1, o
def Hword(n, out, roots, word):
    for (i, mir) in word:
        n, out = append(n, out, roots[i], mir)
    return H_of(n, out)
def test(n, out, roots, letters, L=3):
    H0 = H_of(n, out)
    d = {}
    for l in letters:
        d[l] = Hword(n, out, roots, [l]) - 2 * H0
    for Lw in range(1, L):
        for w in itertools.product(letters, repeat=Lw):
            hw = Hword(n, out, roots, list(w))
            for l in letters:
                if Hword(n, out, roots, list(w) + [l]) - 2 * hw != d[l]:
                    return None
    return H0, d
def tours(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in range(1 << len(pairs)):
        out = [0] * k
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out
if __name__ == '__main__':
    found = 0; tried = 0
    letter_sets = [((0, False), (1, False)), ((0, True), (1, True)), ((0, False), (1, True)), ((0, False), (0, True))]
    for n in range(3, 6):
        seen = set()
        for out in tours(n):
            for x1, x2 in itertools.permutations(range(n), 2):
                for ls in letter_sets:
                    if ls == ((0, False), (0, True)) and x2 != x1 + 1: pass
                    tried += 1
                    r = test(n, out, (x1, x2), ls, L=2)
                    if r is None: continue
                    H0, d = r
                    vals = sorted(set(d.values()))
                    if len(vals) > 1:
                        r3 = test(n, out, (x1, x2), ls, L=3)
                        if r3 is None: continue
                        key = (n, H0, tuple(sorted(d.items())))
                        if key in seen: continue
                        seen.add(key); found += 1
                        if found <= 12:
                            print(f"n={n} H={H0} roots=({x1},{x2}) letters={ls} digits(H(wi)-2H(w))={d} out={out}")
        print(f"n={n} done, odometer seeds so far {found} (tests {tried})", flush=True)
