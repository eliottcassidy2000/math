#!/usr/bin/env python3
"""S93 lane seeds: checks of the one-apex rules of PROOFS.txt section 6 against the complete census d_one_14.npz
(all words of length n <= 14, every O; records with r = +1 or -1).
R1  Theorem 6 (d): for strong words the unit pairs {apex, z_j} with r = +1 are exactly (1 0^(n-1), j = 1),
    (1^(n-1) 0, j = n), (110, j = 1), (100, j = 3); slopes 3 * 2^(n-1), 16, 11.
R2  y = 1^a 00, twins {z_j, z_(j+1)}, j < a: seed, slope 50 * 3^(a-2) - 3 (j = 1), 100 * 2^(j-2) * 3^(a-j-1) (j >= 2).
R3  y = 11 0^b, twins {z_j, z_(j+1)} in the 0-block: seed; slope of the first pair 30 * 3^(b-2) - 1.
R4  twin criterion: consecutive spine vertices with equal letters (a 0-block with prefix Y and suffix Z) form a unit
    pair iff val(Y) = 1 + [y_1 = 1] + [y_n = 0] + J(y) + I(Y) + I(Z)   (all words, n <= 12).
R5  every unit set has at most two vertices.
"""
import numpy as np

d = np.load("d_one_14.npz")
ok1 = ok2 = ok3 = ok4 = ok5 = True
n2 = n3 = n4 = 0
slopes3 = {}
for n in range(2, 15):
    R = d[f"rec_{n}"]
    seeds = {}
    for w, O, sg, M, H in R.tolist():
        if bin(O).count("1") > 2:
            ok5 = False
        if sg == 1:
            seeds[(w, O)] = M
    # R1
    got = set()
    for (w, O), M in seeds.items():
        y = format(w, f"0{n}b")
        if y[0] == "1" and y[-1] == "0" and (O >> n) & 1 and bin(O).count("1") == 2:
            j = [i for i in range(n) if (O >> i) & 1][0] + 1
            got.add((y, j, M))
    exp = {("1" + "0" * (n - 1), 1, 3 << (n - 1)), ("1" * (n - 1) + "0", n, 3 << (n - 1))}
    if n == 3:
        exp |= {("110", 1, 16), ("100", 3, 11)}
    if got != exp:
        ok1 = False
        print("R1 fails", n, got ^ exp)
    # R2
    a = n - 2
    if a >= 2:
        w = int("1" * a + "00", 2)
        for j in range(1, a):
            O = (1 << (j - 1)) | (1 << j)
            M = seeds.get((w, O))
            pred = 50 * 3 ** (a - 2) - 3 if j == 1 else 100 * 2 ** (j - 2) * 3 ** (a - j - 1)
            n2 += 1
            if M != pred:
                ok2 = False
                print("R2 fails", n, j, M, pred)
    # R3
    b = n - 2
    if b >= 2:
        w = int("11" + "0" * b, 2)
        for j in range(3, n):
            O = (1 << (j - 1)) | (1 << j)
            M = seeds.get((w, O))
            n3 += 1
            if M is None:
                ok3 = False
                print("R3 fails", n, j)
            slopes3[(b, j)] = M
        if seeds.get((w, 0b1100)) != 30 * 3 ** (b - 2) - 1:
            ok3 = False
    # R4
    if n <= 12:
        for w in range(1 << n):
            y = [(w >> (n - 1 - i)) & 1 for i in range(n)]
            J = sum(1 for k in range(n - 1) if y[k] == 0 and y[k + 1] == 1)
            def I(seq):
                return sum(1 << (k - i - 1) for i in range(len(seq)) for k in range(i + 1, len(seq)) if seq[i] == 1 and seq[k] == 0)
            for j in range(n - 1):                 # twins at positions j, j+1 (0-based), both letters 0
                if y[j] == 0 and y[j + 1] == 0:
                    Y, Z = y[:j], y[j + 2:]
                    p = sum(Y[i] << (len(Y) - 1 - i) for i in range(len(Y)))
                    C = (1 if y[0] == 1 else 0) + (1 if y[-1] == 0 else 0) + J + I(Y) + I(Z)
                    pred = (p == C + 1)
                    O = (1 << j) | (1 << (j + 1))
                    n4 += 1
                    if pred != ((w, O) in seeds):
                        ok4 = False
print(f"R1 classification of the unit pairs {{apex, spine vertex}} over strong words, n = 2..14: {ok1}")
print(f"R2 slopes of the twins of 1^a 00, a = 2..12 ({n2} pairs): {ok2}")
print(f"R3 twins of 11 0^b are seeds, b = 2..12 ({n3} pairs), first pair has slope 30 * 3^(b-2) - 1: {ok3}")
print("   slopes of the pair (z_j, z_(j+1)) of 11 0^b, b = 2..6:", {k: v for k, v in sorted(slopes3.items()) if k[0] <= 6})
print(f"R4 twin criterion val(Y) = 1 + [y_1 = 1] + [y_n = 0] + J + I(Y) + I(Z) for 0-blocks, all words n <= 12 ({n4} twin pairs): {ok4}")
print(f"R5 every unit set (both signs) over a one-apex base of length <= 14 has at most two vertices: {ok5}")
