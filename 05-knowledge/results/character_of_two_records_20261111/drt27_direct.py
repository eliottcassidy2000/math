#!/usr/bin/env python3
"""Run H3: H mod 256 of doubly regular tournaments of order 27 by direct counting, H = sum over the start vertex v of
the number of Hamiltonian paths of T - v that start in N^+(v) (27 layered programmes on 26 vertices each).

  drt27_direct.py LAYERED control                 P_27 (must give 95)
  drt27_direct.py LAYERED classes MAXMAT STRIDE K the classes of drt_zoo.py g3 MAXMAT STRIDE (same order): the
                                                  self-converse ones and the first K without a converse symmetry

LAYERED = compiled hampaths_layered.c."""
import sys, subprocess
from paley_tools import paley_out
from fold import fold, paths_ending_in
from drt_zoo import gs_solutions, goethals_seidel, descendant, is_drt, classes, anti_involution, is_skew_hadamard

def H_mod256(layered, out):
    n = len(out); tot = 0
    for v in range(n):
        verts = [w for w in range(n) if w != v]; idx = {w: i for i, w in enumerate(verts)}
        lines = ["%d" % (n - 1)]
        for w in verts:
            mask = 0
            for x in verts:
                if x != w and (out[w] >> x) & 1: mask |= 1 << idx[x]
            lines.append("%d %d" % (mask, 1 if (out[v] >> w) & 1 else 0))
        res = subprocess.run([layered], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True).stdout
        tot += int(res.strip().split("\n")[-1].split(":")[1])
    return tot % 256

if __name__ == "__main__":
    layered, mode = sys.argv[1], sys.argv[2]
    if mode == "control":
        h = H_mod256(layered, paley_out(27))
        print("control P_27: H mod 256 = %d (%s)" % (h, "ok" if h == 95 else "MISMATCH")); sys.exit(0)
    maxmat, stride, K = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    reps = []; nsol = 0; nmat = 0
    for a, b, c, d in gs_solutions(7):
        nsol += 1
        if (nsol - 1) % stride or nmat >= maxmat: continue
        M = goethals_seidel(a, b, c, d); assert is_skew_hadamard(M); nmat += 1
        tours = [descendant(M, idx) for idx in range(28)]; assert all(is_drt(o) for o in tours)
        classes(tours, reps)
    print("classes of drt_zoo.py g3 %d %d: %d" % (maxmat, stride, len(reps))); sys.stdout.flush()
    done = 0; ok4 = ok8 = True
    order = [i for i, r in enumerate(reps)]
    info = []
    for i, (inv, rep, members) in enumerate(reps):
        psi = anti_involution(rep); info.append(psi)
    todo = [i for i in order if info[i] is not None] + [i for i in order if info[i] is None][:K]
    for i in todo:
        rep = reps[i][1]; psi = info[i]
        h = H_mod256(layered, rep)
        line = "  class %2d (%d members): self-converse: %s; H = %d (mod 256), %d (mod 4), %d (mod 8)" % (i, len(reps[i][2]), "yes" if psi else "no", h, h % 4, h % 8)
        if psi:
            Z, D1, plus = fold(rep, 27, psi); t = paths_ending_in(Z, D1)
            line += "; folding formula: %d (mod 4) (%s)" % ((1 + 2 * t) % 4, "agrees" if (1 + 2 * t) % 4 == h % 4 else "DISAGREES")
        print(line); sys.stdout.flush()
        ok4 = ok4 and h % 4 == 3; ok8 = ok8 and h % 8 == 7
    print("R9 (H3: H = 3 mod 4): %s; R10 (H3: H = 7 mod 8): %s" % ("HOLDS" if ok4 else "FAILS", "HOLDS" if ok8 else "FAILS"))
