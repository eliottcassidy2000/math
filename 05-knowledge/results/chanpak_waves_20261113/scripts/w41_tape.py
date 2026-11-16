# Wave 41 (mac-mini-2026-11-16-dials). Run in a directory holding lit_base.pkl from w39_lit_setup.py (slowtree) / masks files from gentourng (arc_affine).
# two-apex tapes: spine sites 0..n-1 transitive (i -> j for i < j); letter z = 2a + b: A beats site iff a, B beats site iff b.
import subprocess, sys
HQ = sys.argv[1] if len(sys.argv) > 1 else './hq'
def tour(tape, AB):
    n = len(tape); A, Bv = n, n + 1
    out = [0] * (n + 2)
    for i in range(n):
        for j in range(i + 1, n): out[i] |= 1 << j
    for i, z in enumerate(tape):
        a, b = z >> 1, z & 1
        if a: out[A] |= 1 << i
        else: out[i] |= 1 << A
        if b: out[Bv] |= 1 << i
        else: out[i] |= 1 << Bv
    if AB: out[A] |= 1 << Bv
    else: out[Bv] |= 1 << A
    return out
def H(tapes_ab):
    inp = "".join(f"{len(t)} " + " ".join(map(str, t)) + "\n" for t in (tour(tp, ab) for tp, ab in tapes_ab))
    o = subprocess.run([HQ], input=inp.encode(), capture_output=True, check=True).stdout.decode().split("\n")
    return [int(l.split()[0]) for l in o if l.strip()]
def rc(s): return [3 - z for z in reversed(s)]
def bits(T, k): return [3 if (T >> (k - 1 - i)) & 1 else 0 for i in range(k)]
if __name__ == '__main__':
    k = 4
    for core in ([1, 2], [0, 1, 2, 3]):
        for U in range(8):
            u = [3] + bits(U, 3)
            rows = []
            tapes = []
            for T in range(2 ** k):
                s = bits(T, k)
                W = u + s + core + list(reversed(s)) + rc(u)
                tapes += [(W, True), (W, False)]
            hs = H(tapes)
            up = hs[0::2]; lo = hs[1::2]
            d1 = [lo[T + 1] - lo[T] for T in range(2 ** k - 1)]
            d2 = [d1[i + 1] - d1[i] for i in range(len(d1) - 1)]
            gaps = [up[T] - lo[T] for T in range(2 ** k)]
            print('core', core, 'u', u, 'lower', lo[:6], '... first diffs', d1[:8], 'second diffs nonzero', sum(1 for x in d2 if x), 'gaps', gaps[:5])
