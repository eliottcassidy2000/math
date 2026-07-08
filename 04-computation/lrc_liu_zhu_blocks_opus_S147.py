"""
lrc_liu_zhu_blocks_opus_S147.py   (opus-2026-07-07-S147, HYP-5277 part 2)

THE UNIFORM LOWER-BOUND CONSTRUCTION for Liu-Zhu Conjecture 2.

CONJECTURE (from the S147 pattern mining): for M = {x,y,y-x,y+x}, x=2k+1, y=2m+1,
gcd(x,y)=1, N = 4(k+1)m+1, the optimal Motzkin avoiding set is (up to rotation)

    B = union_{t=0..k} [2m*t, 2m*t + m - 1]     ((k+1) blocks, length m, step 2m)

which lives in [0, 2(k+1)m) = [0, (N-1)/2) (the first half of Z_N), density (k+1)m/N.
Then A = { a^{-1} * b mod N : b in B } avoids M for a SPECIFIC rotation a = a(x,y).

TASK: (1) find, per instance, the rotation(s) a making a*A_something... concretely: find
all c with gcd(c,N)=1 such that C = {c*j mod N : j in B} avoids M; identify c = c(x,y).
(2) If c is a clean function of (x,y): VERIFY avoidance + density symbolically over a big
(x,y) grid -> the general lower bound mu >= (k+1)m/N for ALL Liu-Zhu instances.
(3) Also confirm the x=1 special case (should be the classic slab, m blocks... k=0 => 1
    block of length m, step 2m, in [0,2m) -- a single slab of length m).
"""
from fractions import Fraction as F
from math import gcd


def blocks(k, m, N):
    B = set()
    for t in range(k + 1):
        for s in range(m):
            B.add((2 * m * t + s) % N)
    return B


def avoids(A, M, N):
    Aset = A
    return all(((a + d) % N) not in Aset for a in Aset for d in M)


def find_rotations(k, m, x, y):
    N = 4 * (k + 1) * m + 1
    M = sorted({x, y, y - x, y + x})
    B = blocks(k, m, N)
    good = []
    for c in range(1, N):
        if gcd(c, N) != 1:
            continue
        C = set((c * j) % N for j in B)
        if len(C) == (k + 1) * m and avoids(C, M, N):
            good.append(c)
    return N, M, B, good


def main():
    print("=" * 100)
    print("UNIFORM LOWER-BOUND: A = rotation of (k+1) length-m blocks step 2m; which rotations c avoid M?")
    print("=" * 100)
    print(f"  {'(x,y)':>7} {'M':>16} {'N':>4}  {'#rot avoiding':>13}   {'rotations c (a few)':>24}"
          f"   candidate c(x,y)")
    all_ok = True
    rows = []
    for y in range(3, 26, 2):
        for x in range(1, y, 2):
            if gcd(x, y) != 1 or x + y > 26:
                continue
            k, m = (x - 1) // 2, (y - 1) // 2
            N, M, B, good = find_rotations(k, m, x, y)
            if not good:
                all_ok = False
                print(f"  {str((x,y)):>7} {str(M):>16} {N:>4}  {'NONE!':>13}   -- blocks ansatz FAILS here")
                continue
            # candidate closed forms for c: test x, y, y-x, y+x, 2m, x*inv, etc. mod N
            cands = {}
            for name, val in [("x", x), ("y", y), ("y-x", y - x), ("y+x", y + x),
                              ("2m", 2 * m), ("2k+2", 2 * k + 2), ("-x", (-x) % N),
                              ("x^{-1}", pow(x, -1, N) if gcd(x, N) == 1 else -1),
                              ("y^{-1}", pow(y, -1, N) if gcd(y, N) == 1 else -1)]:
                if val % N in good:
                    cands[name] = val % N
            rows.append((x, y, k, m, N, good, cands))
            print(f"  {str((x,y)):>7} {str(M):>16} {N:>4}  {len(good):>13}   {str(good[:5]):>24}"
                  f"   {cands}")
    print()
    print(f"  BLOCKS ANSATZ avoids M for some rotation on EVERY instance: {all_ok}")
    print()
    print("  Which named closed-form rotation works UNIFORMLY (in EVERY row)?")
    for name in ["x", "y", "y-x", "y+x", "2m", "2k+2", "-x", "x^{-1}", "y^{-1}"]:
        hits = sum(1 for (_, _, _, _, _, _, cands) in rows if name in cands)
        print(f"    c = {name:8s}: works in {hits}/{len(rows)} rows"
              f"{'   <-- UNIFORM' if hits == len(rows) else ''}")
    print()
    print("  If a uniform c is found, the lower bound mu >= (k+1)m/N holds for ALL Liu-Zhu")
    print("  instances by exhibiting A = c^{-1}*B (density (k+1)m/N, avoids M). QED lower half.")


if __name__ == "__main__":
    main()
