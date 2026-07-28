#!/usr/bin/env python3
"""
ak_mode3_glue.py — klein-S691: MODE THREE — identification gluing.

The benchmark's INTUITIVE definition allows, at each gluing stage, either an
edge OR an IDENTIFICATION of corresponding vertices; the verifiable format
has no identification device (and its literal rule-1 is separately unsound —
see ak_exploit_cert.py). Mode three explores the intuitive game: grid
structures plus a merge partition. Merging lowers n (denominator n − t)
while m, r are unchanged — the 9→7 move of the owner's μ₃ fragment.

Implementation: AKInstance on the quotient — vertices = merge classes;
edge vectors act on class ids (an edge whose endpoints merge becomes a
zero/self vector and is dropped from the lattice but STILL COSTS in m,
so honest merges only help when the merged pair is NOT directly edged, or
when the edge's cost is worth the denominator drop).

VALIDATION CAVEAT: constructibility of a given merge pattern (that it
arises from the recursive gluing) is NOT checked here; this is an upper-
bound explorer for what identification could buy. Anything interesting
must be re-derived as an explicit recursive gluing before it is called a
certificate. Statuses: FINITE-EXACT for the quotient forcing runs;
HEURISTIC for constructibility.
"""
import random, math, time
from fractions import Fraction
from itertools import product
from ak_forcing_engine import AKInstance, force
from ak_strict_search import XFULL, POOL3, POOL5, X0, greedy_seeds


class GluedInstance:
    """Quotient of a grid instance by a merge partition (union-find)."""

    def __init__(self, dims, fs, T0, R, merges):
        self.base = AKInstance(XFULL, dims, fs, [], [])
        parent = {v: v for v in self.base.verts}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        for (a, b) in merges:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
        self.rep = {v: find(v) for v in self.base.verts}
        classes = sorted(set(self.rep.values()))
        self.cid = {c: i for i, c in enumerate(classes)}
        self.n = len(classes)
        self.dims, self.fs = dims, fs
        self.T0 = [self.cid[self.rep[tuple(v)]] for v in T0]
        self.R = [(self.cid[self.rep[tuple(v)]], tuple(xx)) for (v, xx) in R]
        self.merges = merges

    def m(self):
        return self.base.m()

    def r(self):
        return len(self.R)

    def t(self):
        return len(set(self.T0))

    def score(self):
        den = self.n - self.t()
        return Fraction(self.m() + self.r(), den) if den > 0 else None

    def generators(self):
        gens = [{v: (a, b)} for (v, (a, b)) in self.R]
        for i0, f in enumerate(self.fs):
            suffixes = list(product(*[range(1, d + 1)
                                      for d in self.dims[i0 + 1:]]))
            for pre, lab in f.items():
                lab = tuple(lab)
                if lab == (0, 0):
                    continue
                prep = pre[:i0] + (pre[i0] + 1,)
                for s in suffixes:
                    u = self.cid[self.rep[pre + s]]
                    w = self.cid[self.rep[prep + s]]
                    if u == w:
                        continue  # merged pair: edge collapses (cost stays)
                    g = {}
                    g[u] = lab
                    g[w] = (-lab[0], -lab[1])
                    gens.append(g)
        return gens


def force_glued(gi):
    """Forcing closure on the quotient (same algorithm as engine.force)."""
    # adapt: build a fake AKInstance-like shim
    class Shim:
        pass
    sh = Shim()
    sh.n = gi.n
    sh.verts = list(range(gi.n))
    sh.vid = {v: v for v in range(gi.n)}
    sh.T0 = gi.T0
    sh.generators = lambda mode: gi.generators()
    return force(sh, "any")


def quick_experiments():
    print("== mode-3 experiments: does identification buy sub-13/7? ==")
    rng = random.Random(7)
    best = (Fraction(10), None)
    # μ₃-flavored [3,3]: three branches of three; merge two leaves of branch 1
    for trial in range(600):
        dims = rng.choice([[3, 3], [2, 2, 2], [2, 2, 2, 2], [4, 3]])
        k = len(dims)
        verts = list(product(*[range(1, d + 1) for d in dims]))
        fs = []
        for i0 in range(k):
            dom = list(product(*([range(1, d + 1) for d in dims[:i0]]
                                 + [range(1, dims[i0])])))
            f = {}
            for pre in dom:
                if rng.random() < 0.7:
                    f[pre] = rng.choice(POOL5)
            fs.append(f)
        nm = rng.choice([1, 1, 2, 3])
        merges = []
        for _ in range(nm):
            a, b = rng.sample(verts, 2)
            merges.append((a, b))
        T0 = [rng.choice(verts)] if rng.random() < 0.6 else []
        # greedy seeds on the quotient: emulate via repeated closure
        gi = GluedInstance(dims, fs, T0, [], merges)
        R = []
        for _ in range(16):
            gi = GluedInstance(dims, fs, T0, R, merges)
            ok, T, _ = force_glued(gi)
            if ok:
                break
            unf = [c for c in range(gi.n) if c not in T]
            # place a seed at first unfired class (map back to any base vertex)
            inv = {}
            for v, rp in gi.rep.items():
                inv.setdefault(gi.cid[rp], v)
            R.append((inv[unf[0]], rng.choice(POOL3)))
        else:
            continue
        if not ok:
            continue
        # prune
        changed = True
        while changed:
            changed = False
            for i in range(len(R) - 1, -1, -1):
                gi2 = GluedInstance(dims, fs, T0, R[:i] + R[i + 1:], merges)
                ok2, _, _ = force_glued(gi2)
                if ok2:
                    R = R[:i] + R[i + 1:]
                    changed = True
        gi = GluedInstance(dims, fs, T0, R, merges)
        ok, _, _ = force_glued(gi)
        s = gi.score()
        if ok and s is not None and s < best[0]:
            best = (s, (dims, fs, T0, R, merges))
            print(f"  trial {trial}: score {s} = {float(s):.4f}  "
                  f"m={gi.m()} r={gi.r()} n={gi.n} t={gi.t()} "
                  f"merges={merges} dims={dims}")
    print("best mode-3 found:", best[0])
    if best[1]:
        dims, fs, T0, R, merges = best[1]
        print(" dims=", dims, "\n fs=", fs, "\n T0=", T0,
              "\n R=", R, "\n merges=", merges)


if __name__ == "__main__":
    quick_experiments()
