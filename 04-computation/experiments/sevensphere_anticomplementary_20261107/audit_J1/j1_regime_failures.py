#!/usr/bin/env python3
# audit J1 (THM-4911 (6) title): diagnose the non-equicentric tournaments (n = 5, 7) with NO vertex in THM-4742's regime: which condition fails at each vertex (O'=O, G'=O', or v,s,w dependent).
import numpy as np
exec(open("j1_transfer_exact.py").read().split("tot_checked = tot_ok = 0")[0])
for n in (5, 7):
    for b, A in enumerate(classes(n)):
        M = A - A.T; x = M.sum(1); y = M @ x
        if (y == y[0]).all(): continue
        info = []
        for v in range(n):
            S = setup(A, v)
            if S["regime"]: break
            why = []
            if S["Op"] == S["O"]: why.append("O'=O")
            if S["Gp"] == S["Op"]: why.append("G'=O'")
            if not why: why.append("v,s,w dependent")
            info.append("/".join(why))
        else:
            arcs = [(i, j) for i in range(n) for j in range(n) if A[i, j]]
            print(f"n={n} class {b}: scores {A.sum(1).tolist()}; per-vertex failure {info}; arcs {arcs}")
