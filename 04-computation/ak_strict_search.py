#!/usr/bin/env python3
"""
ak_strict_search.py — klein-S691 (2026-07-28)

Systematic search for STRICT-mode forcing certificates with score < 2 in
the arithmetic-Kakeya benchmark game (see ak_forcing_engine.py).

GAME RECAP (strict reading): vertices of a d_1×…×d_k grid; axis-i edges
join (p, e_i, s)–(p, e_i+1, s) with a label f_i(p, e_i) ∈ X SHARED across
all suffixes s (cost = number of edges anyway); seeds x·δ_v cost 1;
initial-T wildcards cost via the denominator (n−t).  A vertex v fires
iff the generated lattice contains a vector supported on T ∪ {v} whose
v-column is a nonzero multiple of (1,−1).

STRUCTURE THEORY (hand analysis, klein-S691, recorded for searchers):
  * FLOW PICTURE: lattice vectors = Z²-species flows; an edge labelled ℓ
    transports species ℓ; interior (non-T) vertices need exact species
    cancellation, so straight transport is monochromatic ("telescopes"),
    turning requires ≥3 pairwise-independent incident species satisfying
    a Z-relation (e.g. w = x+y junctions), and T-vertices are universal
    junctions (wildcards).  Firing v = delivering two independent species
    into v (v is where a net (a,−a) — the FORBIDDEN edge species — is
    absorbed; note (1,0)−(0,1)=(1,−1) is exactly the excluded label
    direction, so a fired vertex is an "illegal junction").
  * SLOT BOUND: every fired vertex needs ≥2 incident generators (edges
    at v / seeds at v) with independent directions ⟹ 2m + r ≥ 2(n−t).
    This allows score → 1; the blocking is global (ordering/frontier).
  * FRONTIER CURSE (empirical): on paths, ladders, rail grids, cycles,
    every natural design lands on score exactly 2 — at a straight
    frontier each new vertex has one backward direction and forward
    flows cannot cancel back without a junction species that itself
    exits forward.  CONJECTURE (path case): strict k=1 admits no score
    < 2.  This searcher exists to break or support that conjecture in
    wider classes (junction fabrics, hypercubes, wildcard skeletons).

Everything here is FINITE-EXACT; results print as
    BEST class=<..> score=<Fraction> struct=<repr>
and any score < 2 is re-verified with a fresh closure before reporting.
"""

import random, json, sys, time
from fractions import Fraction
from itertools import product
from ak_forcing_engine import AKInstance, force

X0 = (0, 0)
x, y, w, u, p2 = (1, 0), (0, 1), (1, 1), (1, 2), (2, 1)
POOL3 = [x, y, w]
POOL5 = [x, y, w, u, p2]
XFULL = [X0] + POOL5

TWO = Fraction(2)


def try_force(inst):
    ok, T, _ = force(inst, "strict")
    return ok, T


def greedy_seeds(dims, fs, T0, seed_pool, max_seeds=24):
    """Greedy seed completion: repeatedly add the seed that fires the most
    new vertices (lookahead 1); returns instance or None."""
    R = []
    inst = AKInstance(XFULL, dims, fs, T0, R)
    ok, T = try_force(inst)
    while not ok and len(R) < max_seeds:
        best = None
        unfired = [v for v in inst.verts if inst.vid[v] not in T]
        # candidate seeds only at unfired vertices
        for v in unfired:
            for lab in seed_pool:
                cand = AKInstance(XFULL, dims, fs, T0, R + [(v, lab)])
                ok2, T2 = try_force(cand)
                gain = len(T2) - len(T)
                key = (1 if ok2 else 0, gain)
                if best is None or key > best[0]:
                    best = (key, (v, lab), T2, ok2)
        if best is None or best[0][1] == 0 and not best[3]:
            # no seed helps at all -> add 2 seeds at first unfired vertex
            v = unfired[0]
            R += [(v, x), (v, y)]
        else:
            R.append(best[1])
        inst = AKInstance(XFULL, dims, fs, T0, R)
        ok, T = try_force(inst)
    if not ok:
        return None
    # prune pass: drop seeds that are unnecessary (in reverse add order)
    changed = True
    while changed:
        changed = False
        for i in range(len(R) - 1, -1, -1):
            cand = AKInstance(XFULL, dims, fs, T0, R[:i] + R[i + 1:])
            ok2, _ = try_force(cand)
            if ok2:
                R = R[:i] + R[i + 1:]
                changed = True
    return AKInstance(XFULL, dims, fs, T0, R)


def score_of(inst):
    s = inst.score()
    return s if s is not None else Fraction(10 ** 9)


BEST = {}


def record(cls, inst, note=""):
    s = score_of(inst)
    cur = BEST.get(cls)
    if cur is None or s < cur[0]:
        # re-verify from scratch before recording
        chk = AKInstance(inst.X, inst.dims, inst.fs, inst.T0, inst.R)
        ok, _ = try_force(chk)
        assert ok, "record of non-forcing instance"
        BEST[cls] = (s, inst, note)
        flag = "  *** < 2 ***" if s < TWO else ""
        print(f"  BEST[{cls}] score={s} ({float(s):.4f}) "
              f"m={inst.m()} r={inst.r()} n={inst.n} t={inst.t()}{flag}")
        if s < TWO:
            print("    dims=", inst.dims)
            print("    fs  =", inst.fs)
            print("    T0  =", inst.T0)
            print("    R   =", inst.R)


# ---------------------------------------------------------------- class A
def class_A_paths(nmax=6, pool=POOL3):
    """Exhaustive k=1 paths with labels from pool; T0 in {∅, {v}, pairs};
    greedy seeds.  Conjecture: nothing < 2."""
    print(f"== class A: paths n<= {nmax}, pool={pool}")
    for n in range(2, nmax + 1):
        verts = [(i,) for i in range(1, n + 1)]
        t0_opts = [[]] + [[v] for v in verts]
        if n <= 5:
            t0_opts += [[verts[i], verts[j]] for i in range(n) for j in range(i + 1, n)]
        for labs in product(pool, repeat=n - 1):
            f1 = {(i,): labs[i - 1] for i in range(1, n)}
            for T0 in t0_opts:
                inst = greedy_seeds([n], [f1], T0, POOL3)
                if inst is not None:
                    record("A-path", inst)
    print("   class A done;",
          "min =", BEST.get("A-path", (None,))[0])


# ---------------------------------------------------------------- class B
def class_B_ladders(Lmax=6, iters=400, rng=None):
    """[L,2] ladders: axis-1 labels free per column-step (SHARED top+bottom),
    axis-2 rung labels free per column.  Random patterns + greedy seeds."""
    rng = rng or random.Random(20260728)
    print("== class B: ladders [L,2] random patterns")
    for L in range(3, Lmax + 1):
        for it in range(iters):
            f1 = {}
            for j in range(1, L):
                if rng.random() < 0.95:
                    f1[(j,)] = rng.choice(POOL5)
            f2 = {}
            for j in range(1, L + 1):
                if rng.random() < 0.7:
                    f2[(j, 1)] = rng.choice(POOL5)
            t0_opts = [[], [(1, 1)], [(1, 1), (L, 1)], [(1, 1), (1, 2)]]
            T0 = rng.choice(t0_opts)
            inst = greedy_seeds([L, 2], [f1, f2], T0, POOL3)
            if inst is not None:
                record("B-ladder", inst)
    print("   class B done;", "min =", BEST.get("B-ladder", (None,))[0])


# ---------------------------------------------------------------- class C
def class_C_strips(iters=300, rng=None):
    """[L,3] strips: 4-direction middle-column junction fabrics."""
    rng = rng or random.Random(1675)
    print("== class C: strips [L,3] junction fabrics")
    for L in (3, 4, 5):
        for it in range(iters):
            f1 = {}
            for j in range(1, L):
                if rng.random() < 0.95:
                    f1[(j,)] = rng.choice(POOL5)
            f2 = {}
            for j in range(1, L + 1):
                for e in (1, 2):
                    if rng.random() < 0.65:
                        f2[(j, e)] = rng.choice(POOL5)
            T0 = rng.choice([[], [(1, 2)], [(1, 1)], [(1, 1), (L, 3)]])
            inst = greedy_seeds([L, 3], [f1, f2], T0, POOL3)
            if inst is not None:
                record("C-strip", inst)
    print("   class C done;", "min =", BEST.get("C-strip", (None,))[0])


# ---------------------------------------------------------------- class D
def class_D_cubes(iters=250, rng=None):
    """Hypercubes [2]*k and [3,3,3]: axis-shared labels, sparse entries.
    This is the constructible 'iterated doubling' shape of the KT proofs."""
    rng = rng or random.Random(423)
    print("== class D: hypercubes/deep products")
    shapes = [[2, 2, 2], [2, 2, 2, 2], [3, 3, 3], [2, 3, 2], [3, 2, 2]]
    for dims in shapes:
        k = len(dims)
        for it in range(iters):
            fs = []
            for i0 in range(k):
                dom = list(product(*([range(1, d + 1) for d in dims[:i0]]
                                     + [range(1, dims[i0])])))
                f = {}
                for pre in dom:
                    if rng.random() < rng.choice([0.35, 0.7, 1.0]):
                        f[pre] = rng.choice(POOL5)
                fs.append(f)
            T0 = []
            if rng.random() < 0.5:
                verts = list(product(*[range(1, d + 1) for d in dims]))
                T0 = [rng.choice(verts)]
                if rng.random() < 0.4:
                    T0.append(rng.choice(verts))
                T0 = list(set(T0))
            inst = greedy_seeds(dims, fs, T0, POOL3)
            if inst is not None:
                record(f"D-{'x'.join(map(str,dims))}", inst)
    print("   class D done")


# ---------------------------------------------------------------- hillclimb
def hillclimb(base_cls, steps=400, rng=None):
    """Local moves on the current best of a class: relabel/toggle one entry,
    re-run greedy seeds, keep improvements."""
    if base_cls not in BEST:
        return
    rng = rng or random.Random(89)
    print(f"== hillclimb on {base_cls}")
    s, inst, _ = BEST[base_cls]
    for step in range(steps):
        fs = [dict(f) for f in inst.fs]
        i0 = rng.randrange(len(fs))
        dom = list(product(*([range(1, d + 1) for d in inst.dims[:i0]]
                             + [range(1, inst.dims[i0])])))
        if not dom:
            continue
        pre = rng.choice(dom)
        cur = fs[i0].get(pre, X0)
        choices = [c for c in [X0] + POOL5 if c != cur]
        nv = rng.choice(choices)
        if nv == X0:
            fs[i0].pop(pre, None)
        else:
            fs[i0][pre] = nv
        cand = greedy_seeds(inst.dims, fs, inst.T0, POOL3)
        if cand is not None and score_of(cand) <= score_of(inst):
            if score_of(cand) < score_of(inst):
                record(base_cls, cand, "hillclimb")
            inst = cand
    print("   hillclimb done")


if __name__ == "__main__":
    t0 = time.time()
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "A"):
        class_A_paths(nmax=6)
    if which in ("all", "B"):
        class_B_ladders()
    if which in ("all", "C"):
        class_C_strips()
    if which in ("all", "D"):
        class_D_cubes()
    for cls in list(BEST):
        hillclimb(cls, steps=200)
    print("\n==== FINAL BESTS ====")
    for cls, (s, inst, note) in sorted(BEST.items()):
        print(f"{cls:12s} score={s} ({float(s):.4f}) "
              f"m={inst.m()} r={inst.r()} n={inst.n} t={inst.t()} {note}")
        if s < TWO:
            print("  dims=", inst.dims, " fs=", inst.fs,
                  " T0=", inst.T0, " R=", inst.R)
    print(f"[{time.time()-t0:.1f}s]")
