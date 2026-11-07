#!/usr/bin/env python3
"""Independent audit of THM-4892 (S85), part (E): the face poset of the permutohedron, the topologies J_X
(a sieve on a face F covers iff it contains every edge of F and every 2-face of F whose type is in X),
the presheaf O(B_1|...|B_k) = prod_j Tour(B_j), and its sheafifications.

Checks
  1. J_X is a Grothendieck topology: exhaustive on every face for n = 3 and on every face of dimension
     <= 2 for n = 4; randomised on the 3-dimensional face (the whole permutohedron) for n = 4.
  2. The order between the four topologies, as sets of covering sieves:  J_{sq,hex} is contained in J_{sq}
     and J_{hex}, which are contained in J_{empty}.  With the standard convention (J' is FINER than J iff
     every J-covering sieve is J'-covering; SGA4 II.1, Johnstone C2.1) J_{empty} is the finest and
     J_{sq,hex} the coarsest.  O is a sheaf for the coarsest one only.
  3. O is a presheaf, separated for every X, a sheaf for X = {sq, hex} on every face (n <= 5), and not a
     sheaf for the other X at the top face (n = 4, 5).
  4. The plus construction from the colimit over ALL covering sieves (n = 3, every face; Lemma 0).
  5. Global sections of a(O) = O^+ at the top face: counted from the matching-family relation (opposite
     edges of X-faces carry the same pair orientation); compared with THM-4891's class counts
     (n = 4: 2^24 square-coherent, 2^12 hexagon-coherent, 2^6 tournaments, 2^36 all orientations).
  6. Sink spectra n = 4: hexagon-coherent (exhaustive, 4,096) and square-coherent (exhaustive, 2^24, numpy);
     tournaments {1,3,5}.  THM-4891 states {0..9,12} (7 for exactly 48) and {0..6}.
  7. The n = 5 witness tables of THM-4891 (7 and 21 sinks under each partial coherence) re-counted.
Standard library + numpy.  Imports nothing from the session's scripts.

Run:  python aud_perm_site.py > aud_perm_site.out
"""
import ast
import itertools
import json
import os
import random
import re
import sys
import time
from collections import Counter

T0 = time.time()


# ----------------------------------------------------------------------------------------------------------
# faces = ordered set partitions
# ----------------------------------------------------------------------------------------------------------
def ordered_partitions(items):
    items = tuple(items)
    if not items:
        yield ()
        return
    n = len(items)
    for r in range(1, n + 1):
        for first in itertools.combinations(items, r):
            rest = tuple(x for x in items if x not in first)
            for tail in ordered_partitions(rest):
                yield (tuple(first),) + tail


def refines(Fp, F):
    """Fp <= F  iff  Fp refines F: consecutive blocks of Fp group into the blocks of F, in order."""
    i = 0
    for B in F:
        acc = set()
        while acc != set(B):
            if i >= len(Fp):
                return False
            acc |= set(Fp[i])
            if not acc <= set(B):
                return False
            i += 1
    return i == len(Fp)


def ftype(F):
    sizes = sorted(len(B) for B in F)
    big = [s for s in sizes if s > 1]
    if not big:
        return 'vertex'
    if big == [2]:
        return 'edge'
    if big == [3]:
        return 'hex'
    if big == [2, 2]:
        return 'sq'
    return 'dim>=3'


class FacePoset:
    def __init__(self, n):
        self.n = n
        self.faces = list(ordered_partitions(range(n)))
        self.below = {F: [G for G in self.faces if refines(G, F)] for F in self.faces}
        self.top = (tuple(range(n)),)

    def covers(self, X, F, R):
        for G in self.below[F]:
            t = ftype(G)
            if (t == 'edge' or t in X) and G not in R:
                return False
        return True

    def downsets(self, F):
        """all sieves (down-sets of the faces below F), by antichain enumeration."""
        P = self.below[F]
        le = {(a, b): refines(a, b) for a in P for b in P}
        anti = []

        def rec(i, chosen):
            anti.append(tuple(chosen))
            for j in range(i, len(P)):
                if all(not le[(P[j], c)] and not le[(c, P[j])] for c in chosen):
                    chosen.append(P[j])
                    rec(j + 1, chosen)
                    chosen.pop()
        rec(0, [])
        out = []
        for ac in anti:
            s = set()
            for a in ac:
                s.update(G for G in P if le[(G, a)])
            out.append(frozenset(s))
        return out

    def restrict(self, R, Fp):
        return frozenset(G for G in R if refines(G, Fp))


XS = {'{sq,hex}': ('sq', 'hex'), '{sq}': ('sq',), '{hex}': ('hex',), '{}': ()}


# ----------------------------------------------------------------------------------------------------------
# the presheaf O
# ----------------------------------------------------------------------------------------------------------
def O_elems(F):
    pairs = [p for B in F for p in itertools.combinations(sorted(B), 2)]
    out = []
    for bits in range(1 << len(pairs)):
        out.append(frozenset((a, b) if not (bits >> k & 1) else (b, a) for k, (a, b) in enumerate(pairs)))
    return out


def O_res(F, x, Fp):
    blk = {}
    for i, B in enumerate(Fp):
        for v in B:
            blk[v] = i
    return frozenset(a for a in x if blk[a[0]] == blk[a[1]])


# ----------------------------------------------------------------------------------------------------------
def check_axioms(fp, F, sieves, X, rng=None, sample=None):
    v = Counter()
    covset = {R for R in sieves if fp.covers(X, F, R)}
    if frozenset(fp.below[F]) not in covset:
        v['maximal'] += 1
    for R in covset:
        for Fp in fp.below[F]:
            if not fp.covers(X, Fp, fp.restrict(R, Fp)):
                v['stability'] += 1
    Rs = list(covset)
    S2 = sieves
    if sample:
        Rs = rng.sample(Rs, min(sample, len(Rs)))
        S2 = rng.sample(sieves, min(sample, len(sieves)))
    for R in Rs:
        for R2 in S2:
            if all(fp.covers(X, G, fp.restrict(R2, G)) for G in R) and R2 not in covset:
                v['transitivity'] += 1
    inter = frozenset.intersection(*covset)
    if inter not in covset:
        v['no smallest covering sieve'] += 1
    return v, covset


def random_sieves(fp, F, rng, k):
    """random down-sets of the faces below F, mixed with sieves near R_min so that coverings occur."""
    P = fp.below[F]
    out = set()
    edges_x = [G for G in P if ftype(G) in ('edge', 'sq', 'hex')]
    for _ in range(k):
        gens = [G for G in P if rng.random() < rng.choice([0.05, 0.2, 0.5, 0.8])]
        if rng.random() < 0.6:
            gens += [G for G in edges_x if rng.random() < 0.95]
        s = set()
        for g in gens:
            s.update(G for G in P if refines(G, g))
        out.add(frozenset(s))
    out.add(frozenset(P))
    # the candidate smallest covering sieves (generated by the edges and the X-faces), for every X, and the
    # same sieves with one generator removed (so that the 'smallest covering sieve' test is not vacuous)
    for X in XS.values():
        gens = [G for G in P if ftype(G) == 'edge' or ftype(G) in X]
        for drop in [None] + gens[:: max(1, len(gens) // 12)]:
            s = set()
            for g in gens:
                if g != drop:
                    s.update(G for G in P if refines(G, g))
            out.add(frozenset(s))
    return list(out)


def matching_families_faces(fp, R):
    maxR = [G for G in R if not any(G != H and refines(G, H) for H in R)]
    maxR.sort(key=lambda g: len(g))
    Rs = sorted(R)
    res = []
    assign = {}

    def meet_faces(A, B):
        return [G for G in R if refines(G, A) and refines(G, B)]

    pair_lb = {(i, j): meet_faces(maxR[i], maxR[j]) for i in range(len(maxR)) for j in range(i)}

    def rec(i):
        if i == len(maxR):
            fam = []
            for U in Rs:
                vals = {O_res(M, assign[M], U) for M in maxR if refines(U, M)}
                if len(vals) != 1:
                    return
                fam.append((U, vals.pop()))
            res.append(tuple(fam))
            return
        M = maxR[i]
        for x in O_elems(M):
            ok = True
            for j in range(i):
                M2 = maxR[j]
                for U in pair_lb[(i, j)]:
                    if O_res(M, x, U) != O_res(M2, assign[M2], U):
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                assign[M] = x
                rec(i + 1)
                del assign[M]
    rec(0)
    return res


def colimit_classes(fp, F, covset):
    """|O^+(F)| from the colimit over all covering sieves (equivalence generated by restriction)."""
    cov = list(covset)
    idx = {R: i for i, R in enumerate(cov)}
    fams = [matching_families_faces(fp, R) for R in cov]
    parent = {}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for i in range(len(cov)):
        for f in fams[i]:
            parent[(i, f)] = (i, f)
    for i, R in enumerate(cov):
        maxR = [G for G in R if not any(G != H and refines(G, H) for H in R)]
        for M in maxR:
            j = idx.get(R - {M})
            if j is None:
                continue
            for f in fams[i]:
                g = tuple(p for p in f if p[0] != M)
                a, b = find((i, f)), find((j, g))
                if a != b:
                    parent[a] = b
    classes = len({find(k) for k in parent})
    rmin = frozenset.intersection(*cov)
    return classes, len(fams[idx[rmin]])


# ----------------------------------------------------------------------------------------------------------
# matching families on R_min for O, counted by the forced equalities (opposite edges of X-faces carry the
# same orientation of the same pair).  Every O-section of an X-face restricts injectively to its edges, so
# a matching family on R_min(F) is the same thing as an orientation of the edges of F that is coherent on
# the X-faces; we count classes of edges under the relation and also return the class map (for sinks).
# ----------------------------------------------------------------------------------------------------------
def edge_classes(fp, F, X):
    P = fp.below[F]
    edges = [G for G in P if ftype(G) == 'edge']
    eidx = {e: i for i, e in enumerate(edges)}
    parent = list(range(len(edges)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for G in P:
        if ftype(G) in X:
            es = [e for e in P if ftype(e) == 'edge' and refines(e, G)]
            bypair = {}
            for e in es:
                pair = next(tuple(sorted(B)) for B in e if len(B) == 2)
                bypair.setdefault(pair, []).append(e)
            for pair, lst in bypair.items():
                for e in lst[1:]:
                    a, b = find(eidx[lst[0]]), find(eidx[e])
                    if a != b:
                        parent[a] = b
    cls = {}
    for e in edges:
        cls.setdefault(find(eidx[e]), []).append(e)
    return edges, [find(eidx[e]) for e in edges], len(cls)


def part_axioms_and_order():
    print('=' * 100)
    print('1-4. Axioms, the order between the four topologies, O, and the colimit plus construction')
    rng = random.Random(3)
    for n in (3, 4):
        fp = FacePoset(n)
        tot = Counter()
        order = Counter()
        sheafness = Counter()
        lemma0 = Counter()
        for F in fp.faces:
            dim = n - len(F)
            exhaustive = (n == 3) or dim <= 2
            sieves = fp.downsets(F) if exhaustive else random_sieves(fp, F, rng, 4000)
            covs = {}
            for name, X in XS.items():
                v, covset = check_axioms(fp, F, sieves, X, rng, None if exhaustive else 300)
                for k, c in v.items():
                    tot[f'{name} {k}'] += c
                covs[name] = covset
            # order between topologies on this face (as sets of covering sieves)
            order['J{sq,hex} <= J{sq}'] += covs['{sq,hex}'] <= covs['{sq}']
            order['J{sq,hex} <= J{hex}'] += covs['{sq,hex}'] <= covs['{hex}']
            order['J{sq} <= J{}'] += covs['{sq}'] <= covs['{}']
            order['J{hex} <= J{}'] += covs['{hex}'] <= covs['{}']
            order['faces'] += 1
            for a_, b_ in (('{sq,hex}', '{}'), ('{sq}', '{}'), ('{hex}', '{}'), ('{sq,hex}', '{sq}'),
                           ('{sq,hex}', '{hex}')):
                if covs[a_] < covs[b_]:
                    order[f'strict somewhere: J{a_} < J{b_}'] = 1
            # colimit plus construction (n = 3 every face; n = 4 faces of dim <= 2)
            if exhaustive:
                for name, X in XS.items():
                    cl, rm = colimit_classes(fp, F, covs[name])
                    lemma0['Lemma 0 (colimit = R_min families): ' + str(cl == rm)] += 1
                    nO = len(O_elems(F))
                    sheafness[f'{name}: O(F) -> O^+(F) bijective on {ftype(F) if dim <= 2 else "top"}: ' +
                              str(nO == cl)] += 1
        print(f'n={n}: faces {len(fp.faces)}; axiom violations: {dict(tot) or "none"}')
        print(f'      order of topologies (counts of faces where the inclusion holds / total): {dict(order)}')
        print(f'      {dict(lemma0)}')
        for k in sorted(sheafness):
            print(f'      {k}: {sheafness[k]}')


def part_sections():
    print('=' * 100)
    print('3/5. O separated for every X; O a sheaf exactly for X = {sq,hex}; global sections of a(O) = O^+')
    for n in (3, 4, 5):
        fp = FacePoset(n)
        rep = Counter()
        for F in fp.faces:
            P = fp.below[F]
            edges = [G for G in P if ftype(G) == 'edge']
            # separated: an O-section is determined by its restrictions to the edges
            secs = O_elems(F)
            sig = {tuple(O_res(F, x, e) for e in edges) for x in secs}
            rep['O separated (sections determined by edges): ' + str(len(sig) == len(secs))] += 1
            for name, X in XS.items():
                _, _, ncls = edge_classes(fp, F, X)
                rep[f'{name}: |O^+(F)| = 2^#classes = |O(F)|: ' + str(2 ** ncls == len(secs))] += 1
        top = fp.top
        cnt = {name: edge_classes(fp, top, X)[2] for name, X in XS.items()}
        print(f'n={n}: ' + '; '.join(f'{k}: {v}' for k, v in sorted(rep.items())))
        print(f'      top face: number of free edge classes (global sections = 2^classes): {cnt}')


def sinks_spectrum(fp, X):
    n = fp.n
    top = fp.top
    edges, cls, ncls = edge_classes(fp, top, X)
    roots = sorted(set(cls))
    ridx = {r: i for i, r in enumerate(roots)}
    # canonical orientation bit per class: bit 1 means "smaller element first" for the class representative
    # pair; all edges of a class swap the same pair, so the bit is well defined.
    eclass = {}
    for e, r in zip(edges, cls):
        eclass[e] = ridx[r]
    verts = [G for G in fp.faces if len(G) == n]
    # for each vertex, the 3 (n-1) incident edges and the required bit for the vertex to be a sink
    req = []
    for v in verts:
        arr = [B[0] for B in v]
        conds = []
        for i in range(n - 1):
            a, b = arr[i], arr[i + 1]
            e = tuple(tuple([x]) for x in arr[:i]) + (tuple(sorted((a, b))),) + tuple(tuple([x]) for x in arr[i + 2:])
            c = eclass[e]
            conds.append((c, 1 if a < b else 0))     # sink needs: a (= pi_i) first
        req.append(conds)
    return ncls, req


def part_spectra():
    import numpy as np
    print('=' * 100)
    print('6. Sink spectra at n = 4 (orientations coherent on the X-faces = global sections of a(O))')
    fp = FacePoset(4)
    for name in ('{hex}', '{sq}', '{sq,hex}'):
        ncls, req = sinks_spectrum(fp, XS[name])
        N = 1 << ncls
        cnt = np.zeros(N, dtype=np.int32)
        ar = np.arange(N, dtype=np.uint32)
        for conds in req:
            ok = np.ones(N, dtype=bool)
            for c, bit in conds:
                ok &= (((ar >> np.uint32(c)) & np.uint32(1)) == np.uint32(bit))
            cnt += ok
        spec = Counter(cnt.tolist())
        print(f'   X = {name}: {ncls} classes, {N} orientations; sink spectrum {sorted(spec)}; '
              f'orientations with 7 sinks: {spec.get(7, 0)}')


def part_witnesses():
    print('=' * 100)
    print("7. THM-4891's n = 5 witness tables re-counted (sinks under the stated relaxation)")
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(os.path.join(here, '..', '..', 'one_apex_collision_balance_20261106', 'permutohedron',
                                         'witnesses.out'))
    lines = open(path, encoding='utf-8').read().splitlines()
    for ln in lines:
        m = re.match(r'table (\S+) \(pair\|(.*?) -> first element\): (\{.*\})\s*$', ln)
        if not m:
            continue
        tag, kind, js = m.group(1), m.group(2), m.group(3)
        table = json.loads(js)
        n = int(tag[1])
        sinks = []
        for pi in itertools.permutations(range(n)):
            ok = True
            for i in range(n - 1):
                a, b = pi[i], pi[i + 1]
                pair = ''.join(map(str, sorted((a, b))))
                if 'set before' in kind:
                    before = ''.join(map(str, sorted(pi[:i])))
                    key = f'{pair}|{{{before}}}'
                else:
                    others = ''.join(str(x) for x in pi if x not in (a, b))
                    key = f'{pair}|{others}'
                if table[key] != a:
                    ok = False
                    break
            if ok:
                sinks.append(''.join(map(str, pi)))
        print(f'   {tag:9s} ({"square/set-context" if "set before" in kind else "hexagon/relative-order"} table, '
              f'{len(table)} classes): {len(sinks)} sinks')


if __name__ == '__main__':
    part_axioms_and_order()
    part_sections()
    part_spectra()
    part_witnesses()
    print(f'total time {time.time() - T0:.1f}s')
