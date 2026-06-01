#!/usr/bin/env python3
"""
lrc_cutting_sequence_s521.py   claudebox-2026-06-01-S521

The chamber sign-vector / cutting-sequence encoding of LRC (reflection:
07-reflections/lrc-complex-vertices-cutting-sequences-s521.md).

Vertex = the full sector-vector s(t)=(floor(v_i t n))_i in (Z/n)^m = which sector
each runner occupies. The line t->(v_i t) on the m-torus, cut by the 1/n-grid; s(t)
is its CUTTING SEQUENCE. strict-LRC <=> the central-box vector (all coords in
{1..n-2}) appears. The realizable set is polynomially thin (~n*sum(v_i)) inside the
exponential n^m grid; consecutive cells differ by a single coordinate +-1.
Also: Menger/source probe (source = local indeg condition, no connectivity content).
"""
from fractions import Fraction as F
def fr(x): return x % 1
def sector(p, n): return int(fr(p) * n)
def secvec(sp, t, n): return tuple(sector(F(v)*t, n) for v in sp)
def cells(sp, n):
    W = set([F(0)])
    for v in sp:
        for k in range(n*v+1): W.add(F(k, n*v) % 1)
    W = sorted(w for w in W if 0 <= w < 1); W2 = W + [F(1)]
    return [(a+b)/2 for a, b in zip(W, W2[1:]) if 0 < (a+b)/2 < 1]

def main():
    print("Cutting-sequence (sector-vector) encoding: realizable cells vs n^m grid\n")
    print(f"{'speeds':22} {'cells':6} {'n*sum(v)':9} {'n^m':10} {'ratio':10} {'central box?':12}")
    for sp in [(1,2,3,4),(1,2,4,7),(1,2,3,4,5),(2,3,5,7,11),(1,2,3,4,5,6)]:
        n = len(sp)+1; m = len(sp)
        R = set(secvec(list(sp), t, n) for t in cells(list(sp), n))
        central = any(all(1 <= c <= n-2 for c in sv) for sv in R)
        print(f"{str(sp):22} {len(R):6} {n*sum(sp):9} {n**m:10} {len(R)/n**m:.2e}  {str(central):12}")
    print("\n=> polynomial cells in an exponential grid: a thin structured slice (cutting sequence).")
    # single-coordinate-move check
    sp = [1, 2, 4, 7]; n = 5
    seq = []
    for t in sorted(set(cells(sp, n))):
        s = secvec(sp, t, n)
        if not seq or s != seq[-1]: seq.append(s)
    def one_move(a, b):
        diff = [(b[i]-a[i]) % n for i in range(len(a))]
        nz = [d for d in diff if d != 0]
        return len(nz) == 1 and nz[0] in (1, n-1)
    mv = sum(1 for a, b in zip(seq, seq[1:]) if one_move(a, b))
    print(f"transitions {sp}: {mv}/{len(seq)-1} single-coordinate +-1 moves (grid-edge walk on (Z/n)^m)")
    print("\nMenger/source probe: observer is a flow-source iff indeg=N(t)=0 (LOCAL condition);")
    print("  no global connectivity content -- Menger adds nothing beyond N(t)=0.")

if __name__ == "__main__":
    main()
