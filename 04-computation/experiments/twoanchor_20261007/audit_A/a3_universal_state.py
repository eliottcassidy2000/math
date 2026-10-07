#!/usr/bin/env python3
"""Audit A, item 3: THM-4601 (iii) universal residual state, independent code.

* clearing identity F_(2,2,2)(x) - 1 = 27 (F_(4,1,1)(y_3) - 1) as affine maps, and the three formulas of the proof;
* validity: for x = 1 mod 8 with 27 | x+1, y_3 = (x+1)/27 - 1 has U-letters (4,1,1) iff x = 1 mod 128
  (exhaustive over x mod 27*2^12), and y_4 has (2,2,1,1) iff x = 1 mod 128 (81 | x+1);
* random residual sources (own generator), J = 3 (j = 7, 8) and large J, K = 4, 5 and large K:
  at Terras time 2J after x the D = 3 chain is in (3, -26) with Y = T^2J(y_3) = T^2J(y_4) (D = 4 needs K >= 5);
* D = 5, 6: state (6, 1 - 3^6) once J >= 6, and the (J-dependent) states for J = 3, 4, 5.
"""
import random
from fractions import Fraction as Fr
from collections import defaultdict
from a2_barrier import make_source, T

NCHK = 0
def check(c, m):
    global NCHK
    if not c: raise AssertionError(m)
    NCHK += 1

def v2(n):
    return (n & -n).bit_length() - 1

def Fw(w, x):
    for a in w: x = (3*x + 1) / Fr(2**a)
    return x

def uword(z, L):
    out = []
    for _ in range(L):
        y = 3*z + 1; a = v2(y); out.append(a); z = y >> a
    return tuple(out), z

def part_identity():
    for x in (Fr(0), Fr(1), Fr(12345, 7), Fr(-3)):
        y3 = (x + 1)/27 - 1; y4 = (x + 1)/81 - 1
        check(Fw((2, 2, 2), x) - 1 == 27*(Fw((4, 1, 1), y3) - 1), "clearing identity")
        check(Fw((4,), y3) == (x - 17)/144 and Fw((4, 1), y3) == (x + 31)/96 and Fw((4, 1, 1), y3) == (x + 63)/64, "proof formulas")
        check(Fw((2, 2, 2), x) == (27*x + 37)/64, "F_222")
        check(Fw((2, 2, 1, 1), y4) == (x + 63)/64, "y_4 via (2,2,1,1)")
    # exhaustive validity over x mod 27*2^12 (x = 1 mod 8, x = -1 mod 27)
    M2 = 1 << 12
    cnt = {True: 0, False: 0}
    for xr in range(1, M2, 8):
        # CRT: x = xr mod 2^12, x = -1 mod 27
        x = xr + M2 * ((-1 - xr) * pow(M2, -1, 27) % 27)
        x += 27 * M2 * 1000                         # positive lift, letters below 2^12 decided
        check((x + 1) % 27 == 0 and x % M2 == xr, "CRT")
        y3 = (x + 1)//27 - 1
        w3, _ = uword(y3, 3)
        cond = (x % 128 == 1)
        check((w3 == (4, 1, 1)) == cond, "y_3 letters (4,1,1) iff x = 1 mod 128")
        cnt[cond] += 1
        if (x + 1) % 81 == 0:
            y4 = (x + 1)//81 - 1
            w4, _ = uword(y4, 4)
            check((w4 == (2, 2, 1, 1)) == cond, "y_4 letters (2,2,1,1) iff x = 1 mod 128")
        w0, _ = uword(x, 3)
        check((w0 == (2, 2, 2)) == (x % 128 == 1), "x letters (2,2,2) iff x = 1 mod 128 (j >= 7)")
    print(f"   validity exhaustive mod 27*2^12: (4,1,1) iff x=1 mod 128 ({cnt[True]} classes yes, {cnt[False]} no); same for y_4 (2,2,1,1)")

def chain_state_after(x, y, D, steps):
    u, v, k = x, y, D
    for _ in range(steps):
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
    return k, u, v

def part_sources():
    rnd = random.Random(31337)
    nsrc = 0
    Ks = [4, 5, 6, 7, 8, 13, 40, 200, 777, 1500]
    js = [7, 8, 9, 10, 11, 12, 13, 14, 20, 33, 64, 120, 301, 302]
    states56 = defaultdict(set)
    for K in Ks:
        for j in js:
            for _ in range(6):
                n, t, x = make_source(rnd, K, j, extra_bits=120)
                J = (j - 1)//2
                check(J >= 3, "J >= 3")
                y3 = (x + 1)//27 - 1
                k3, X3, Y3 = chain_state_after(x, y3, 3, 2*J)
                check(k3 == 3 and X3 - 1 == 27*(Y3 - 1), "D=3 state (3,-26)")
                check(X3 - 27*Y3 == -26, "e = -26")
                # words: x has 2^J, y_3 has (4,1,1) 2^(J-3)
                wx, _ = uword(x, J + 1); wy, _ = uword(y3, J + 1)
                check(wx[:J] == (2,)*J and wx[J] != 2, "x has exactly J letters 2")
                check(wy[:3] == (4, 1, 1) and wy[3:J] == (2,)*(J - 3) and wy[J] != 2, "y_3 word (4,1,1) 2^(J-3), then not 2")
                if K >= 5:
                    y4 = (x + 1)//81 - 1
                    k4, X4, Y4 = chain_state_after(x, y4, 4, 2*J)
                    check(k4 == 3 and X4 == X3 and Y4 == Y3, "D=4 lands in the same state and same Y")
                if K >= 7:
                    for D in (5, 6):
                        yD = (x + 1)//3**D - 1
                        kD, XD, YD = chain_state_after(x, yD, D, 2*J)
                        states56[(D, min(J, 7))].add((kD, XD - 3**kD*YD))
                nsrc += 1
    print(f"   {nsrc} random sources, K in {Ks}, j in {js} (J = 3..150): D=3 state (3,-26) and words ok; "
          f"D=4 same state and same Y (K >= 5)")
    for key in sorted(states56):
        print(f"   D={key[0]}, J={'>=7' if key[1] == 7 else key[1]}: end-of-run states {sorted(states56[key])}")
    for D in (5, 6):
        check(states56[(D, 7)] == {(6, 1 - 729)}, "D=5,6 universal (6, 1-3^6) for J >= 6 (sampled J >= 7)")

if __name__ == "__main__":
    part_identity()
    part_sources()
    # J = 6 exactly for D = 5, 6
    rnd = random.Random(99)
    for j in (13, 14):
        for K in (7, 9, 30, 300):
            n, t, x = make_source(rnd, K, j, 100)
            for D in (5, 6):
                yD = (x + 1)//3**D - 1
                k, X, Y = chain_state_after(x, yD, D, 2*((j - 1)//2))
                check((k, X - 3**k*Y) == (6, -728), "J = 6: (6, -728)")
    print(f"   J = 6 exactly: D = 5, 6 in (6, -728)")
    print(f"ALL (iii) CHECKS PASSED ({NCHK} assertions)")
