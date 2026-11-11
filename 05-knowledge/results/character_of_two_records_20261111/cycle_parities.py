#!/usr/bin/env python3
"""Run I: parities of the numbers of directed cycles of every odd length on skew Hadamard switching classes.
Conjecture (HYP-9667 (bK)): if S^2 = -(N-1) I then c_(2k+1) = [N = 4 mod 8] C(N/2 - 1, k) (mod 2) for every k.

  cycle_parities.py CYC controls          P_7, P_11, P_19, P_23: c_3, Hamiltonian cycles, 1 + 2 c_odd = H (mod 4)
  cycle_parities.py CYC insample          P_3, P_7 (all members), P_11, CD_4 classes (random members)
  cycle_parities.py CYC local             the vertex-local form fails (CD_4, D_19)
  cycle_parities.py CYC i1 JSON           order 23: the 13 classes of JSON (drt23_classes.json of the Sierpinski lane)
                                          and the double of P_11
  cycle_parities.py CYC i1mod4 JSON       the same classes: odd cycle counts modulo 4 (exploratory)
  cycle_parities.py CYC short JSON        cycle counts of lengths 3..12 modulo 256, orders 19 and 23 (exploratory)
  cycle_parities.py CYC pairs NAMES       all single and double switchings of (source) + T, T in P_7, P_11, CD_4,
                                          P_19, D_19, P_23, DP_11 (the double of P_11)
  cycle_parities.py CYC i2 [members]      N = 24: random members of the classes of P_23^+ and D(P_11)^+
  cycle_parities.py CYC i3 [members]      N = 20: D_19, P_19 (with a source) and random members of their classes
  cycle_parities.py CYC i4 MAXMAT STRIDE K  order 27: P_27 and the classes of drt27_direct.py (same order)

CYC = compiled cycle_lengths.c.  A tournament T on n vertices is doubly regular iff (source) + T satisfies
S^2 = -(N-1) I with N = n + 1; (source) + T has the cycles of T."""
import sys, random, subprocess, json
from math import comb
from paley_tools import paley_out
from cayley_dickson import tournament
from drt_zoo import (gs_solutions, goethals_seidel, descendant, double, paley_matrix, is_drt, classes, anti_involution,
                     is_skew_hadamard, find_map)
from switch_class import plus, switch, is_sh

def cycles(exe, out):
    inp = "%d\n" % len(out) + "".join("%d\n" % o for o in out)
    res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
    c = {}
    for l in res:
        if l.strip(): a, b = l.split(); c[int(a)] = int(b)
    return c
def odd_vector(c, N):
    return ''.join(str(c[l] % 2) for l in range(3, N + 1, 2))
def predicted(N):
    """N = order of the skew Hadamard tournament (even); lengths 3, 5, ..., N - 1."""
    m = N // 2 - 1
    return ''.join(str((comb(m, k) % 2) if N % 8 == 4 else 0) for k in range(1, m + 1))
def code_to_out(code, n):
    out = [0] * n; it = iter(code)
    for i in range(n):
        for j in range(i + 1, n):
            if next(it) == '1': out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out
def check(name, exe, out, N_sh, results):
    """out: a tournament; N_sh: the order of the skew Hadamard tournament it represents (len(out) or len(out)+1)."""
    c = cycles(exe, out); v = odd_vector(c, len(out)); p = predicted(N_sh)
    v = v + '0' * (len(p) - len(v)) if len(v) < len(p) else v       # a member of order N has no cycle of length N + 1
    ok = (v == p); results.append(ok)
    even = ''.join(str(c[l] % 2) for l in range(4, len(out) + 1, 2))
    print("  %-28s odd lengths 3,5,...: %s  predicted %s  %s   (even lengths 4,6,...: %s)" % (name, v, p, "ok" if ok else "FAILS", even)); sys.stdout.flush()
    return c

if __name__ == "__main__":
    exe, mode = sys.argv[1], sys.argv[2]; results = []
    if mode == "controls":
        known = {7: (189, 24), 11: (95095, 5505), 19: (1172695746915, 34358763933), 23: (15760206976379349, 374127973957716)}
        for q, (H, C) in known.items():
            c = cycles(exe, paley_out(q)); c3 = q * (q * q - 1) // 24; codd = sum(c[l] for l in range(3, q + 1, 2))
            print("control P_%d: c_3 = %d (mod 256) against %d: %s; Hamiltonian cycles %d (mod 256) against %d: %s; 1 + 2 c_odd = %d (mod 4) against H: %d: %s" % (
                q, c[3], c3 % 256, "ok" if c[3] == c3 % 256 else "MISMATCH", c[q], C % 256, "ok" if c[q] == C % 256 else "MISMATCH",
                (1 + 2 * codd) % 4, H % 4, "ok" if (1 + 2 * codd) % 4 == H % 4 else "MISMATCH"))
    elif mode == "insample":
        random.seed(9)
        for name, out in [("P_3", paley_out(3)), ("P_7", paley_out(7)), ("P_11", paley_out(11)), ("CD_4", tournament(4)[0])]:
            T = plus(out); N = len(T); assert is_sh(T)
            Xs = list(range(1 << (N - 1))) if N <= 8 else [0] + [random.randrange(1 << (N - 1)) for _ in range(40 if N <= 12 else 12)]
            vs = {}
            for X in Xs:
                v = odd_vector(cycles(exe, switch(T, X)), N); vs[v] = vs.get(v, 0) + 1
            print("class of %s^+ (N = %d): %d members; parity vectors of the odd cycle counts: %s; predicted %s" % (name, N, len(Xs), sorted(vs.items()), predicted(N)))
    elif mode == "local":
        # the vertex-local form (the cycles through each single vertex) is false
        def through(out, v):
            n = len(out); perm = [v] + [w for w in range(n) if w != v]; pos = {w: i for i, w in enumerate(perm)}; o = [0] * n
            for x in range(n):
                for y in range(n):
                    if x != y and (out[x] >> y) & 1: o[pos[x]] |= 1 << pos[y]
            inp = "%d\n" % n + "".join("%d\n" % t for t in o)
            res = subprocess.run([exe, "root0"], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
            c = {int(l.split()[0]): int(l.split()[1]) for l in res if l.strip()}
            return odd_vector(c, n)
        a, b, c, d = next(gs_solutions(5)); G = descendant(goethals_seidel(a, b, c, d), 0)
        for name, out in [("P_7", paley_out(7)), ("P_11", paley_out(11)), ("CD_4", tournament(4)[0]), ("D_19", G)]:
            vs = {}
            for v in range(len(out)):
                t = through(out, v); vs[t] = vs.get(t, 0) + 1
            print("%s: parity vectors of the numbers of odd cycles through a vertex, by vertex: %s; global prediction %s" % (name, sorted(vs.items()), predicted(len(out) + 1)))
    elif mode == "i1":
        d = json.load(open(sys.argv[3]))
        print("I1: order 23 (N = 24): predicted " + predicted(24))
        for rec in d:
            out = code_to_out(rec['code'], 23); assert is_drt(out)
            check("%s (converse class %s)" % (rec['label'], rec['converse_class']), exe, out, 24, results)
        D = descendant(double(paley_matrix(11)), 0); assert is_drt(D); check("double of P_11", exe, D, 24, results)
        print("R12 (I1): %s (%d of %d)" % ("HOLDS" if all(results) else "FAILS", sum(results), len(results)))
    elif mode == "i1mod4":
        # exploratory: the odd cycle counts modulo 4 and the implied parity of the number of disjoint pairs
        d = json.load(open(sys.argv[3]))
        rows = [(rec['label'], code_to_out(rec['code'], 23), rec['H']) for rec in d] + [("double", descendant(double(paley_matrix(11)), 0), 15759752857737405)]
        for name, out, H in rows:
            c = cycles(exe, out); odd = [c[l] for l in range(3, 24, 2)]; a1 = sum(odd)
            print("%-7s c_l mod 4 (l = 3, 5, ..., 23): %s   c_odd mod 8 = %d   H mod 8 = %d, mod 16 = %d   parity of the number of disjoint pairs of odd cycles = %d" % (
                name, ''.join(str(x % 4) for x in odd), a1 % 8, H % 8, H % 16, ((H - 1 - 2 * a1) // 4) % 2))
    elif mode == "short":
        # the cycle counts of small length are the same for all doubly regular tournaments of one order
        d = json.load(open(sys.argv[3]))
        rows = [(rec['label'], code_to_out(rec['code'], 23)) for rec in d] + [("double", descendant(double(paley_matrix(11)), 0))]
        print("order 23: cycle counts modulo 256 for the lengths 3, 4, ..., 12")
        for name, out in rows:
            c = cycles(exe, out); print("  %-7s %s" % (name, ' '.join("%3d" % c[l] for l in range(3, 13))))
        a, b, c0, d0 = next(gs_solutions(5)); G = descendant(goethals_seidel(a, b, c0, d0), 0)
        print("order 19: the same")
        for name, out in [("P_19", paley_out(19)), ("D_19", G)]:
            c = cycles(exe, out); print("  %-7s %s" % (name, ' '.join("%3d" % c[l] for l in range(3, 13))))
    elif mode == "pairs":
        # all single and double switchings of (source) + T for the classes on which the law was seen to hold:
        # if H mod 4 and the cycle parities are quadratic functions of the switching set (observation QF), a class
        # on which all singles and pairs agree with the normalised member is constant.
        names = sys.argv[3:]
        a, b, c0, d0 = next(gs_solutions(5)); G = descendant(goethals_seidel(a, b, c0, d0), 0)
        pool = {"P_7": paley_out(7), "P_11": paley_out(11), "CD_4": tournament(4)[0], "P_19": paley_out(19), "D_19": G,
                "P_23": paley_out(23), "DP_11": descendant(double(paley_matrix(11)), 0)}
        for name in names:
            T = plus(pool[name]); N = len(T); base = odd_vector(cycles(exe, T), N); vs = {}
            for u in range(N):
                for v in range(u, N):
                    X = (1 << u) | (1 << v)
                    w = odd_vector(cycles(exe, switch(T, X)), N); vs[w] = vs.get(w, 0) + 1
            print("class of (source) + %s (N = %d): %d single and double switchings; odd parity vectors: %s; normalised member %s; predicted %s" % (
                name, N, N + N * (N - 1) // 2, sorted(vs.items()), base, predicted(N))); sys.stdout.flush()
    elif mode == "i2":
        k = int(sys.argv[3]) if len(sys.argv) > 3 else 10; random.seed(2412)
        print("I2: N = 24: predicted " + predicted(24))
        D = descendant(double(paley_matrix(11)), 0)
        for name, out in [("P_23", paley_out(23)), ("D(P_11)", D)]:
            T = plus(out); assert is_sh(T)
            for i in range(k): check("%s^+ member %d" % (name, i), exe, switch(T, random.randrange(1 << 23)), 24, results)
        print("R12 (I2): %s (%d of %d)" % ("HOLDS" if all(results) else "FAILS", sum(results), len(results)))
    elif mode == "i3":
        k = int(sys.argv[3]) if len(sys.argv) > 3 else 20; random.seed(2012)
        print("I3: N = 20: predicted " + predicted(20))
        a, b, c, d = next(gs_solutions(5)); G = descendant(goethals_seidel(a, b, c, d), 0); assert is_drt(G) and find_map(G, paley_out(19)) is None
        for name, out in [("D_19", G), ("P_19", paley_out(19))]:
            check(name, exe, out, 20, results); T = plus(out); assert is_sh(T)
            for i in range(k): check("%s^+ member %d" % (name, i), exe, switch(T, random.randrange(1 << 19)), 20, results)
        print("R12 (I3): %s (%d of %d)" % ("HOLDS" if all(results) else "FAILS", sum(results), len(results)))
    elif mode == "i4":
        maxmat, stride, K = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
        print("I4: order 27 (N = 28): predicted " + predicted(28))
        check("P_27", exe, paley_out(27), 28, results)
        reps = []; nsol = 0; nmat = 0
        for a, b, c, d in gs_solutions(7):
            nsol += 1
            if (nsol - 1) % stride or nmat >= maxmat: continue
            M = goethals_seidel(a, b, c, d); assert is_skew_hadamard(M); nmat += 1
            tours = [descendant(M, idx) for idx in range(28)]; classes(tours, reps)
        info = [anti_involution(r[1]) for r in reps]
        todo = [i for i in range(len(reps)) if info[i] is not None] + [i for i in range(len(reps)) if info[i] is None][:K]
        for i in todo:
            check("class %d (%s)" % (i, "self-converse" if info[i] else "no converse symmetry"), exe, reps[i][1], 28, results)
        print("R12 (I4): %s (%d of %d)" % ("HOLDS" if all(results) else "FAILS", sum(results), len(results)))
