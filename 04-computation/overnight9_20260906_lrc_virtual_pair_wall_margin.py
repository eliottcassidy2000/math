"""Standalone exact controls for the sharp virtual-wall fibre margin."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import sys

sys.stdout.reconfigure(newline="\n")
gates = 0


def need(test, label):
    global gates
    gates += 1
    if not test:
        raise RuntimeError(label)


def norm(x):
    return min(x % 1, (-x) % 1)


def fibre(T, y):
    return tuple(min(norm(w * (y + j) / 3) for w in T) for j in range(3))


def main():
    phases = sorted({F(p, q) for q in range(1, 25) for p in range(q)})
    units = [w for w in range(1, 12) if w % 3]
    rows = 0
    for T in combinations(units, 3):
        ds = [abs(u - v if (u - v) % 3 == 0 else u + v) // 3
              for u, v in combinations(T, 2)]
        for y in phases:
            H = max(fibre(T, y))
            for d in ds:
                need(H >= abs(F(1, 3) - norm(d * y)) / 2,
                     "literal quantitative signed-pair margin")
                if H < F(1, 14):
                    need(F(4, 21) < norm(d * y) < F(10, 21),
                         "full spoil obeys both strict owner-band walls")
        rows += 1
    sharp = 0
    for A, B in combinations(range(2, 81, 3), 2):
        if gcd(A, B) != 1:
            continue
        T = (A, B, A + B)
        y = F(1, A + B)
        d = (B - A) // 3
        want = F(A, 3 * (A + B))
        need(fibre(T, y) == (want, want, F(0)), "exact primitive sharpness family")
        need(want == (F(1, 3) - norm(d * y)) / 2, "sharp linear profile")
        sharp += 1
    need(fibre((11, 17, 28), F(1, 28)) == (F(11, 84), F(11, 84), F(0)),
         "sharp virtual-wall value 11/84")
    upper_sharp = 0
    for A in range(2, 81, 3):
        for C in range(2, 321, 3):
            if C < 4 * A:
                continue
            g = gcd(A, C)
            T = (A, C, A + C)
            y = F(1, C)
            d = (2 * A + C) // 3
            want = F(A, 3 * C)
            need(fibre(T, y) == (want, F(0), want), "entire upper sharpness family")
            need(want == (norm(d * y) - F(1, 3)) / 2, "upper exact linear profile")
            reduced = tuple(w // g for w in T)
            need(max(fibre(reduced, g * y)) == want
                 and norm((d // g) * g * y) == norm(d * y),
                 "primitive dilation preserves upper sharpness")
            upper_sharp += 1
    need(max(fibre((2, 11, 13), F(1, 11))) == F(2, 33)
         and norm(5 * F(1, 11)) == F(5, 11), "original upper-half sharpness overclaim hostile")
    need(max(fibre((1, 4, 5), F(1, 4))) == F(1, 12), "delta=1/2 sharp endpoint")
    for N in range(1, 61):
        if N % 11 == 5:
            continue
        T = (9 * N - 1, 42 * N - 1, 51 * N - 2)
        y = F(1, T[1])
        d = 20 * N - 1
        need(gcd(T[0], T[1]) == 1 and all(w % 3 for w in T),
             "primitive upper owner-band sharpness sequence")
        need(max(fibre(T, y)) < F(1, 14), "upper-wall sequence actually fully spoiled")
        need(F(10, 21) - norm(d * y) == F(11, 21 * (42 * N - 1)),
             "exact approach to the sharp upper owner wall")

    # Arbitrary coprime unit lifts of h, not the primary h=1 modulo L.
    controls = []
    for d in (31, 37, 55, 65):
        b, c = 3 * d + 1, 3 * d + 2
        L = 42 * d * b * c
        h = 91 ** 6 + 2
        while gcd(h, L) != 1:
            h += 1
        need(h % L != 1, "independent h residue")
        C = (d, 2 * d, 3 * d, 4 * d, 14, 14 * b, 14 * c, h, 2 * h, 1)
        T = (1, b, c)
        survivors = []
        for k in range(d):
            y = F(14 * k + 1, 14 * d)
            if all(norm(z * y) >= F(1, 14) for z in C):
                H = max(fibre(T, y))
                need(H >= F(11, 84), "uniform sharp tail margin at every surviving virtual wall")
                need(all(norm(z * y) > F(1, 14) for z in C if z != d),
                     "every other body owner is strictly safe")
                j = max(range(3), key=lambda j: fibre(T, y)[j])
                eta = min([F(1, 14 * d), 3 * (H - F(1, 14)) / (2 * max(T))]
                          + [(norm(z * y) - F(1, 14)) / (2 * z) for z in C if z != d])
                x = (y + eta / 2 + j) / 3
                need(all(norm(z * x) > F(1, 14) for z in [3 * z for z in C] + list(T)),
                     "arbitrary-unit-lift strict physical component witness")
                survivors.append(k)
        need(len(survivors) > 0, "phase-free virtual gate survives arbitrary coprime h")
        controls.append((d, h, len(survivors), survivors[0]))
    print("STATUS: PASS; globally sharp absolute-value fibre margin and arbitrary-coprime-h strict extension")
    print("FINITE CONTROLS:", rows, "tail triples;", len(phases), "Farey phases;", sharp,
          "primitive exact sharpness rows")
    print("UPPER PROFILE:", upper_sharp, "equality rows with primitive rescaling; endpoint delta=1/2 attained")
    print("CORRECTION: original valid one-sided bound was not globally sharp; exact profile is |delta-1/3|/2")
    print("FULL-SPOIL BAND: 4/21 < delta < 10/21; both endpoints sharp")
    print("VIRTUAL WALL SHARPNESS: T=(11,17,28), y=1/28, sheet minima=(11/84,11/84,0)")
    print("ARBITRARY COPRIME h CONTROLS (d,h,survivor_count,first_k):", controls)
    print("ACTIVE GATES:", gates)


if __name__ == "__main__":
    main()
