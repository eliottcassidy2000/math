#!/usr/bin/env python3
"""
ak_mode3_v2.py — klein-S691: LEGAL mode-three (the intuitive gluing game).

The intuitive benchmark definition, faithfully: at stage i, consecutive
copies H_{e_i}, H_{e_i+1} are glued; FOR EACH corresponding vertex pair
(i.e. each suffix s, given the prefix p): choose
    'M'  — identify the two vertices (merge),
    a label v ∈ X — connect by an edge labelled v (per-slot label!), or
    None — no connection.
So unlike the verifiable format, labels vary per suffix AND merges exist;
m counts actual edges; n counts vertices AFTER merging.

v1 (ak_mode3_glue.py) allowed arbitrary merges — ILLEGAL (its 11/6 used
cross-suffix and non-adjacent merges). Here merges are restricted to
grid-edge slots ((p, e_i, s) ~ (p, e_i+1, s)), which is exactly what the
recursion permits. NOTE: a faithful stage-wise recursion would also demand
the copies at stage i be *isomorphic as glued objects below stage i*; the
per-prefix freedom used here matches the verifiable format's f_i(prefix
dependence) and is if anything SMALLER than the fully general recursion
with distinct H-copies; flag any record for a re-derivation pass.

Statuses: forcing runs FINITE-EXACT; constructibility LEGAL-BY-SLOT
(stage-wise re-derivation still owed before calling anything a
benchmark-grade certificate).
"""
import random, math, time, sys
from fractions import Fraction
from itertools import product
from ak_forcing_engine import force
from ak_strict_search import XFULL, POOL3, POOL5, X0

MERGE = "M"


class Mode3Instance:
    """dims + slots: slots[i0][(prefix, suffix)] ∈ {label, MERGE}; absent = none."""

    def __init__(self, dims, slots, T0, R):
        self.dims, self.slots = dims, slots
        self.baseverts = list(product(*[range(1, d + 1) for d in dims]))
        parent = {v: v for v in self.baseverts}

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        self.edges = []  # (u, v, label) on base vertices
        k = len(dims)
        for i0 in range(k):
            for (pre, s), val in slots[i0].items():
                u = pre + s
                w = pre[:i0] + (pre[i0] + 1,) + s
                if val == MERGE:
                    ra, rb = find(u), find(w)
                    if ra != rb:
                        parent[ra] = rb
                else:
                    self.edges.append((u, w, tuple(val)))
        self.rep = {v: find(v) for v in self.baseverts}
        classes = sorted(set(self.rep.values()))
        self.cid = {c: i for i, c in enumerate(classes)}
        self.n = len(classes)
        self.verts = list(range(self.n))
        self.vid = {v: v for v in self.verts}
        self.T0base = [tuple(v) for v in T0]
        self.Rbase = [(tuple(v), tuple(xx)) for (v, xx) in R]
        self.T0 = sorted(set(self.cid[self.rep[tuple(v)]] for v in T0))
        self.R = [(self.cid[self.rep[tuple(v)]], tuple(xx)) for (v, xx) in R]

    def m(self):
        # edges whose endpoints ended up merged still count (they were chosen
        # as edges); count all label slots
        return len(self.edges)

    def r(self):
        return len(self.R)

    def t(self):
        return len(self.T0)

    def score(self):
        den = self.n - self.t()
        return Fraction(self.m() + self.r(), den) if den > 0 else None

    def generators(self, mode="any"):
        gens = [{v: lab} for (v, lab) in self.R]
        for (u, w, lab) in self.edges:
            cu, cw = self.cid[self.rep[u]], self.cid[self.rep[w]]
            if cu == cw:
                continue
            gens.append({cu: lab, cw: (-lab[0], -lab[1])})
        return gens


def try_force(mi):
    ok, T, _ = force(mi, "any")
    return ok, T


def greedy(dims, slots, T0, rng, max_seeds=18):
    R = []
    for _ in range(max_seeds):
        mi = Mode3Instance(dims, slots, T0, R)
        ok, T = try_force(mi)
        if ok:
            break
        inv = {}
        for v, rp in mi.rep.items():
            inv.setdefault(mi.cid[rp], v)
        unf = [c for c in range(mi.n) if c not in T]
        # greedy-1: best (vertex,label) by fired gain
        best = None
        for c in unf[:8]:
            for lab in POOL3:
                mi2 = Mode3Instance(dims, slots, T0, R + [(inv[c], lab)])
                ok2, T2 = try_force(mi2)
                key = (1 if ok2 else 0, len(T2) - len(T))
                if best is None or key > best[0]:
                    best = (key, (inv[c], lab))
        R.append(best[1])
    else:
        return None
    # prune
    changed = True
    while changed:
        changed = False
        for i in range(len(R) - 1, -1, -1):
            mi2 = Mode3Instance(dims, slots, T0, R[:i] + R[i + 1:])
            ok2, _ = try_force(mi2)
            if ok2:
                R = R[:i] + R[i + 1:]
                changed = True
    return Mode3Instance(dims, slots, T0, R)


def anneal(dims, steps, seed, log=print):
    rng = random.Random(seed)
    k = len(dims)
    slotkeys = []
    for i0 in range(k):
        pres = list(product(*([range(1, d + 1) for d in dims[:i0]]
                              + [range(1, dims[i0])])))
        sufs = list(product(*[range(1, d + 1) for d in dims[i0 + 1:]]))
        slotkeys.append([(p, s) for p in pres for s in sufs])
    verts = list(product(*[range(1, d + 1) for d in dims]))

    def rand_slots():
        slots = [dict() for _ in range(k)]
        for i0 in range(k):
            for key in slotkeys[i0]:
                import os
                z = rng.random()
                if z < 0.45:
                    slots[i0][key] = rng.choice(POOL5)
                elif z < 0.58 and not os.environ.get('AK_NO_MERGE'):
                    slots[i0][key] = MERGE
        return slots

    slots, T0 = rand_slots(), []
    cur = greedy(dims, slots, T0, rng)
    cur_s = cur.score() if cur else Fraction(10)
    best, best_s = cur, cur_s
    for step in range(steps):
        temp = 0.3 * (1 - step / steps) + 0.02
        ns = [dict(d) for d in slots]
        nT0 = list(T0)
        mv = rng.random()
        if mv < 0.75:
            i0 = rng.randrange(k)
            if not slotkeys[i0]:
                continue
            key = rng.choice(slotkeys[i0])
            cval = ns[i0].get(key)
            import os
            opts = ([None] + POOL5) if os.environ.get('AK_NO_MERGE') else [None, MERGE] + POOL5
            nv = rng.choice([o for o in opts if o != cval])
            if nv is None:
                ns[i0].pop(key, None)
            else:
                ns[i0][key] = nv
        elif mv < 0.88 and len(nT0) < 2:
            v = rng.choice(verts)
            if v not in nT0:
                nT0.append(v)
        elif nT0:
            nT0.pop(rng.randrange(len(nT0)))
        cand = greedy(dims, ns, nT0, rng)
        if cand is None:
            continue
        s = cand.score()
        if s is None:
            continue
        if s <= cur_s or rng.random() < math.exp(float(cur_s - s) / temp):
            slots, T0, cur, cur_s = ns, nT0, cand, s
            if s < best_s:
                chk = Mode3Instance(dims, ns, nT0, cand.Rbase)
                ok, _ = try_force(chk)
                assert ok
                best, best_s = cand, s
                log(f"  [{dims} seed{seed} step{step}] NEW BEST {s} "
                    f"= {float(s):.4f} m={cand.m()} r={cand.r()} "
                    f"n={cand.n} t={cand.t()}")
                if s <= Fraction(11, 6):
                    log(f"   *** rung <= 11/6: slots={ns}")
                    log(f"   T0={nT0} Rbase={cand.Rbase}")
    return best, best_s


if __name__ == "__main__":
    t0 = time.time()
    overall = (Fraction(10), None, None)
    for i, (dims, steps) in enumerate([([3, 3], 1500), ([2, 2, 2], 1200),
                                       ([3, 3], 1500), ([2, 2, 2, 2], 600),
                                       ([3, 3, 3], 350), ([4, 3], 900)]):
        print(f"== mode3-legal anneal {dims} steps={steps}")
        b, s = anneal(dims, steps, 300 + i)
        print(f"   -> best {s} ({float(s):.4f})  [{time.time()-t0:.0f}s]")
        if b is not None and s < overall[0]:
            overall = (s, b, dims)
    s, b, dims = overall
    print("\n==== MODE-3 LEGAL OVERALL BEST ====")
    if b:
        print(f"score={s} = {float(s):.4f}  dims={dims} "
              f"m={b.m()} r={b.r()} n={b.n} t={b.t()}")
        print("slots=", b.slots)
        print("T0(base)=", b.T0, " R(classids)=", b.R)
