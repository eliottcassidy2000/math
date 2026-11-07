#!/usr/bin/env python3
"""Independent audit of THM-4892 (S85), parts (A) and (B): Grothendieck-topology axioms and the plus
construction computed FROM THE DEFINITION, i.e. as the colimit over ALL covering sieves of the sets of
matching families (not only over the smallest covering sieve).

Sites (poset of subsets of V = {0..n-1}, morphisms = inclusions; sieves = down-sets):
  J2 : a sieve R on S covers iff R contains every subset of S with at most two elements;
  Jss: (for a fixed tournament T) R covers S iff R contains every strong component of T[S].

Everything is written from scratch: tournaments are out-neighbour bitmasks, strong components come from a
bitmask transitive closure, odd cycles from a depth-first search rooted at the least vertex, sieves from an
antichain enumeration of the Boolean lattice (Dedekind numbers 2,3,6,20,168,7581 are asserted).
Standard library only.  Imports nothing from the session's scripts.

Run:  python aud_sites_colimit.py > aud_sites_colimit.out
"""
import itertools
import random
import sys
import time
from collections import Counter

T0 = time.time()


# ----------------------------------------------------------------------------------------------------------
# bit utilities
# ----------------------------------------------------------------------------------------------------------
def popcount(x):
    return bin(x).count('1')


def elems_of(mask):
    out, i = [], 0
    while mask >> i:
        if mask >> i & 1:
            out.append(i)
        i += 1
    return out


def submasks(mask):
    res, s = [], mask
    while True:
        res.append(s)
        if s == 0:
            break
        s = (s - 1) & mask
    return res


def is_sub(a, b):
    return a & ~b == 0


# ----------------------------------------------------------------------------------------------------------
# down-sets of the Boolean lattice B_k (all sieves on a k-element object)
# ----------------------------------------------------------------------------------------------------------
_DS = {}
DEDEKIND = {0: 2, 1: 3, 2: 6, 3: 20, 4: 168, 5: 7581}


def downsets_local(k):
    if k in _DS:
        return _DS[k]
    N = 1 << k
    anti = []

    def comparable(a, b):
        return (a & b) == a or (a & b) == b

    def rec(i, chosen):
        anti.append(tuple(chosen))
        for j in range(i, N):
            if all(not comparable(j, c) for c in chosen):
                chosen.append(j)
                rec(j + 1, chosen)
                chosen.pop()
    rec(0, [])
    ds = []
    for ac in anti:
        cl = set()
        for a in ac:
            cl.update(submasks(a))
        ds.append(frozenset(cl))
    assert len(set(ds)) == len(ds) == DEDEKIND[k], (k, len(ds))
    _DS[k] = ds
    return ds


def sieves_on(S):
    el = elems_of(S)
    k = len(el)

    def glob(b):
        g = 0
        for i in range(k):
            if b >> i & 1:
                g |= 1 << el[i]
        return g
    return [frozenset(glob(b) for b in ds) for ds in downsets_local(k)]


def restrict_sieve(R, Sp):
    return frozenset(U for U in R if is_sub(U, Sp))


# ----------------------------------------------------------------------------------------------------------
# tournaments
# ----------------------------------------------------------------------------------------------------------
def tour_from_bits(n, bits):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> k & 1:
                out[j] |= 1 << i
            else:
                out[i] |= 1 << j
            k += 1
    return out


def random_tour(n, rng):
    return tour_from_bits(n, rng.getrandbits(n * (n - 1) // 2) if n > 1 else 0)


def components(out, S):
    """strong components of T[S] as bitmasks, via transitive closure of reachability inside S."""
    el = elems_of(S)
    reach = {v: (out[v] & S) | (1 << v) for v in el}
    changed = True
    while changed:
        changed = False
        for v in el:
            r = reach[v]
            acc = r
            for w in elems_of(r):
                acc |= reach[w]
            if acc != r:
                reach[v] = acc
                changed = True
    comps, seen = [], 0
    for v in el:
        if seen >> v & 1:
            continue
        c = 0
        for w in el:
            if (reach[v] >> w & 1) and (reach[w] >> v & 1):
                c |= 1 << w
        comps.append(c)
        seen |= c
    return comps


def odd_cycles(out, n):
    """directed odd cycles (length >= 3) as (vertex mask, tuple starting at least vertex); each once."""
    cyc = []
    for s in range(n):
        def dfs(path, used):
            last = path[-1]
            if len(path) >= 3 and len(path) % 2 == 1 and (out[last] >> s & 1):
                cyc.append((used, tuple(path)))
            for w in range(s + 1, n):
                if not (used >> w & 1) and (out[last] >> w & 1):
                    path.append(w)
                    dfs(path, used | (1 << w))
                    path.pop()
        dfs([s], 1 << s)
    return cyc


def ham_count(out, S):
    el = elems_of(S)
    k = len(el)
    if k == 0:
        return 1
    idx = {v: i for i, v in enumerate(el)}
    dp = [[0] * k for _ in range(1 << k)]
    for i in range(k):
        dp[1 << i][i] = 1
    for M in range(1 << k):
        for i in range(k):
            x = dp[M][i]
            if x == 0:
                continue
            for w in elems_of(out[el[i]] & S):
                j = idx[w]
                if not (M >> j & 1):
                    dp[M | (1 << j)][j] += x
    return sum(dp[(1 << k) - 1])


# ----------------------------------------------------------------------------------------------------------
# sites
# ----------------------------------------------------------------------------------------------------------
class Site:
    def __init__(self, n, covers, name):
        self.n, self.covers, self.name = n, covers, name
        self.objs = list(range(1 << n))
        self.sieves = {S: sieves_on(S) for S in self.objs}
        self.cov = {S: [R for R in self.sieves[S] if covers(S, R)] for S in self.objs}

    def check_axioms(self, full_transitivity=True, sample=400, rng=None):
        v = Counter()
        for S in self.objs:
            subs = submasks(S)
            maxs = frozenset(subs)
            if not self.covers(S, maxs):
                v['maximal sieve does not cover'] += 1
            covset = set(self.cov[S])
            # upward closure (consequence of the axioms; used by the Hasse-edge colimit below)
            for R in self.cov[S]:
                for R2 in self.sieves[S]:
                    if R <= R2 and R2 not in covset:
                        v['not upward closed'] += 1
            # stability under pullback
            for R in self.cov[S]:
                for Sp in subs:
                    if not self.covers(Sp, restrict_sieve(R, Sp)):
                        v['stability fails'] += 1
            # transitivity (local character)
            if full_transitivity:
                table = {}
                for R2 in self.sieves[S]:
                    table[R2] = {U: self.covers(U, restrict_sieve(R2, U)) for U in subs}
                for R in self.cov[S]:
                    for R2 in self.sieves[S]:
                        if all(table[R2][U] for U in R) and R2 not in covset:
                            v['transitivity fails'] += 1
            else:
                rng = rng or random.Random(0)
                sv = self.sieves[S]
                for R in rng.sample(self.cov[S], min(len(self.cov[S]), 60)):
                    for R2 in rng.sample(sv, min(len(sv), sample)):
                        if all(self.covers(U, restrict_sieve(R2, U)) for U in R) and R2 not in covset:
                            v['transitivity fails'] += 1
            # smallest covering sieve
            inter = frozenset.intersection(*self.cov[S]) if self.cov[S] else None
            if inter is None or inter not in covset:
                v['no smallest covering sieve'] += 1
        return v

    def rmin(self, S):
        return frozenset.intersection(*self.cov[S])

    def irreducible(self, S):
        return len(self.cov[S]) == 1

    def empty_covered(self, S):
        return frozenset() in self.cov[S]


# ----------------------------------------------------------------------------------------------------------
# presheaves:  P = (elems: dict S -> list, res: function (S, x, Sp) -> element of P(Sp))
# ----------------------------------------------------------------------------------------------------------
def check_functorial(objs, P):
    elems, res = P
    bad = 0
    for S in objs:
        for x in elems[S]:
            if res(S, x, S) != x:
                bad += 1
            for Sp in submasks(S):
                y = res(S, x, Sp)
                if y not in elems[Sp]:
                    bad += 1
                    continue
                for Spp in submasks(Sp):
                    if res(Sp, y, Spp) != res(S, x, Spp):
                        bad += 1
    return bad


def matching_families(R, P):
    """all matching families on the sieve R (a down-set), as full families ((U, x_U) for every U in R),
    found by backtracking on the maximal elements and then VERIFIED on every element of R."""
    elems, res = P
    maxR = [U for U in R if not any(U != W and is_sub(U, W) for W in R)]
    maxR.sort(key=lambda m: (-popcount(m), m))
    Rs = sorted(R)
    out = []
    assign = {}

    def rec(i):
        if i == len(maxR):
            fam = []
            for U in Rs:
                vals = {res(M, assign[M], U) for M in maxR if is_sub(U, M)}
                if len(vals) != 1:
                    return
                fam.append((U, vals.pop()))
            out.append(tuple(fam))
            return
        M = maxR[i]
        for x in elems[M]:
            ok = True
            for j in range(i):
                M2 = maxR[j]
                I = M & M2
                if res(M, x, I) != res(M2, assign[M2], I):
                    ok = False
                    break
            if ok:
                assign[M] = x
                rec(i + 1)
                del assign[M]
    rec(0)
    return out


class UF:
    def __init__(self):
        self.p = {}

    def add(self, a):
        self.p.setdefault(a, a)

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def plus_colimit(site, P, check_welldef=True):
    """F^+(S) = colim over ALL covering sieves R of S of Match(R, F), the transition maps being restriction
    of families along R' subset R.  The equivalence is generated by single steps R -> R minus {maximal
    element} that stay covering (enough, since the covering sieves are upward closed - checked above).
    Returns (P_plus, info) with info[S] = dict(classes, rmin_count, lemma0_ok)."""
    elems, res = P
    store = {}
    info = {}
    for S in site.objs:
        cov = site.cov[S]
        Ridx = {R: i for i, R in enumerate(cov)}
        fams = [matching_families(R, P) for R in cov]
        uf = UF()
        for i in range(len(cov)):
            for f in fams[i]:
                uf.add((i, f))
        for i, R in enumerate(cov):
            maxR = [U for U in R if not any(U != W and is_sub(U, W) for W in R)]
            for M in maxR:
                R2 = R - {M}
                j = Ridx.get(R2)
                if j is None:
                    continue
                for f in fams[i]:
                    g = tuple(p for p in f if p[0] != M)
                    assert (j, g) in uf.p, 'restricted family is not a matching family?'
                    uf.union((i, f), (j, g))
        roots = {}
        node_class = {}
        for node in uf.p:
            r = uf.find(node)
            if r not in roots:
                roots[r] = len(roots)
            node_class[node] = roots[r]
        rmin = site.rmin(S)
        im = Ridx[rmin]
        per_class = Counter(node_class[(im, f)] for f in fams[im])
        lemma0 = (len(fams[im]) == len(roots)) and all(per_class[c] == 1 for c in range(len(roots)))
        maxidx = Ridx[frozenset(submasks(S))]
        rmin_of_class = {node_class[(im, f)]: f for f in fams[im]}
        store[S] = dict(cov=cov, Ridx=Ridx, fams=fams, node_class=node_class, nclass=len(roots),
                        maxidx=maxidx, rmin_of_class=rmin_of_class)
        info[S] = dict(classes=len(roots), rmin_count=len(fams[im]), lemma0_ok=lemma0)
    # restriction maps of F^+
    rest = {}
    welldef_bad = 0
    for S in site.objs:
        st = store[S]
        reps = {}
        for node, c in st['node_class'].items():
            reps.setdefault(c, []).append(node)
        for Sp in submasks(S):
            stp = store[Sp]
            for c, nodes in reps.items():
                targets = set()
                use = nodes if check_welldef else nodes[:1]
                for (i, f) in use:
                    R = st['cov'][i]
                    Rp = restrict_sieve(R, Sp)
                    j = stp['Ridx'][Rp]
                    g = tuple(p for p in f if is_sub(p[0], Sp))
                    targets.add(stp['node_class'][(j, g)])
                if len(targets) != 1:
                    welldef_bad += 1
                rest[(S, c, Sp)] = targets.pop()
    plus_elems = {S: list(range(store[S]['nclass'])) for S in site.objs}

    def plus_res(S, c, Sp):
        return rest[(S, c, Sp)]

    # unit eta: F -> F^+
    def eta(S, x):
        st = store[S]
        f = tuple((U, res(S, x, U)) for U in sorted(submasks(S)))
        return st['node_class'][(st['maxidx'], f)]
    def decode(S, c):
        return store[S]['rmin_of_class'][c]
    return (plus_elems, plus_res), info, eta, welldef_bad, decode


# ----------------------------------------------------------------------------------------------------------
# concrete presheaves
# ----------------------------------------------------------------------------------------------------------
def constant(A, objs):
    return ({S: list(A) for S in objs}, lambda S, x, Sp: x)


def linear_orders(objs):
    el = {S: [tuple(p) for p in itertools.permutations(elems_of(S))] for S in objs}
    return (el, lambda S, x, Sp: tuple(v for v in x if Sp >> v & 1))


def tournaments_presheaf(objs):
    el = {}
    for S in objs:
        e = elems_of(S)
        pairs = list(itertools.combinations(e, 2))
        L = []
        for bits in range(1 << len(pairs)):
            L.append(frozenset((a, b) if not (bits >> k & 1) else (b, a) for k, (a, b) in enumerate(pairs)))
        el[S] = L
    return (el, lambda S, x, Sp: frozenset(a for a in x if (Sp >> a[0] & 1) and (Sp >> a[1] & 1)))


def set_partitions(items):
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        for i in range(len(p)):
            yield tuple(sorted(p[:i] + (tuple(sorted((first,) + p[i])),) + p[i + 1:]))
        yield tuple(sorted(p + ((first,),)))


def partitions_presheaf(objs):
    el = {S: sorted(set(set_partitions(elems_of(S)))) for S in objs}

    def res(S, x, Sp):
        return tuple(sorted(tuple(v for v in b if Sp >> v & 1) for b in x if any(Sp >> v & 1 for v in b)))
    return (el, res)


def product(P1, P2, objs):
    e1, r1 = P1
    e2, r2 = P2
    return ({S: [(a, b) for a in e1[S] for b in e2[S]] for S in objs},
            lambda S, x, Sp: (r1(S, x[0], Sp), r2(S, x[1], Sp)))


def forget_below(P, Z, k, objs):
    """F(S) = P(S) x Z when |S| >= k, and P(S) when |S| < k (the Z-label is forgotten on small sets).
    A presheaf; for |Z| >= 2 it is NOT separated for J2 once k = 3 (the label is invisible on pairs)."""
    e, r = P
    el = {S: ([(x, z) for x in e[S] for z in Z] if popcount(S) >= k else [(x, None) for x in e[S]]) for S in objs}

    def res(S, x, Sp):
        return (r(S, x[0], Sp), x[1] if popcount(Sp) >= k else None)
    return (el, res)


def packings_presheaf(out, n, objs):
    cyc = odd_cycles(out, n)
    el = {}
    for S in objs:
        inside = [c for c in cyc if is_sub(c[0], S)]
        packs = []

        def rec(i, used, chosen):
            packs.append(tuple(chosen))
            for j in range(i, len(inside)):
                if inside[j][0] & used == 0:
                    chosen.append(inside[j])
                    rec(j + 1, used | inside[j][0], chosen)
                    chosen.pop()
        rec(0, 0, [])
        L = []
        for p in packs:
            for cols in itertools.product((0, 1), repeat=len(p)):
                L.append(tuple(sorted((c[0], c[1], col) for c, col in zip(p, cols))))
        el[S] = L

    def res(S, x, Sp):
        return tuple(t for t in x if is_sub(t[0], Sp))
    return (el, res)


def swap(x):
    return tuple(sorted((a, b, 1 - c) for a, b, c in x))


def q_presheaf(out, n, objs):
    el, res = packings_presheaf(out, n, objs)

    def canon(x):
        return min(x, swap(x))
    qel = {S: sorted(set(canon(x) for x in el[S])) for S in objs}
    return (qel, lambda S, x, Sp: canon(res(S, x, Sp)))


def random_quotient_presheaf(objs, Y, Z, rng, k_ident):
    """X(S) = Y^S x Z (functions S -> Y plus a global label), divided by a random equivalence relation
    that is closed under restriction (so the quotient is a presheaf).  Typically NOT separated."""
    X = {}
    for S in objs:
        e = elems_of(S)
        X[S] = [(tuple(zip(e, vals)), z) for vals in itertools.product(Y, repeat=len(e)) for z in Z]

    def xres(S, x, Sp):
        return (tuple(p for p in x[0] if Sp >> p[0] & 1), x[1])
    uf = {S: UF() for S in objs}
    for S in objs:
        for x in X[S]:
            uf[S].add(x)
    for _ in range(k_ident):
        S = rng.choice(objs)
        a, b = rng.choice(X[S]), rng.choice(X[S])
        uf[S].union(a, b)
    changed = True
    while changed:
        changed = False
        for S in objs:
            for x in X[S]:
                r = uf[S].find(x)
                if r == x:
                    continue
                for v in elems_of(S):
                    Sp = S & ~(1 << v)
                    a, b = xres(S, x, Sp), xres(S, r, Sp)
                    if uf[Sp].find(a) != uf[Sp].find(b):
                        uf[Sp].union(a, b)
                        changed = True
    el = {S: sorted(set(uf[S].find(x) for x in X[S])) for S in objs}
    return (el, lambda S, x, Sp: uf[Sp].find(xres(S, x, Sp)))


# ----------------------------------------------------------------------------------------------------------
# generic sheaf-theoretic checks for one presheaf on one site
# ----------------------------------------------------------------------------------------------------------
def analyse(site, P, label):
    elems, res = P
    rep = Counter()
    assert check_functorial(site.objs, P) == 0, label + ': not a presheaf'
    Pp, info, eta, wd, decode = plus_colimit(site, P)
    rep['F+ restriction well defined'] += (wd == 0)
    rep['F+ restriction ill defined'] += (wd != 0)
    assert check_functorial(site.objs, Pp) == 0
    Ppp, info2, eta2, wd2, _ = plus_colimit(site, Pp)
    lemma0 = all(info[S]['lemma0_ok'] for S in site.objs) and all(info2[S]['lemma0_ok'] for S in site.objs)
    sep = sheaf = True
    psep = psheaf = True
    per_S = {}
    for S in site.objs:
        imgs = [eta(S, x) for x in elems[S]]
        inj = len(set(imgs)) == len(imgs)
        sur = set(imgs) == set(Pp[0][S])
        imgs2 = [eta2(S, c) for c in Pp[0][S]]
        inj2 = len(set(imgs2)) == len(imgs2)
        sur2 = set(imgs2) == set(Ppp[0][S])
        sep &= inj
        sheaf &= inj and sur
        psep &= inj2
        psheaf &= inj2 and sur2
        per_S[S] = dict(F=len(elems[S]), Fp=info[S]['classes'], Fpp=info2[S]['classes'], inj=inj, sur=sur,
                        img=set(imgs), inj2=inj2, sur2=sur2)
    return dict(lemma0=lemma0, separated=sep, sheaf=sheaf, plus_separated=psep, plus_sheaf=psheaf,
                per_S=per_S, Pp=Pp, Ppp=Ppp, eta=eta, eta2=eta2, info=info, welldef=(wd == 0 and wd2 == 0),
                decode=decode)


# ----------------------------------------------------------------------------------------------------------
# PART A: the pair site J2
# ----------------------------------------------------------------------------------------------------------
def covers_J2(S, R):
    return all(U in R for U in submasks(S) if popcount(U) <= 2)


def part_A(nmax=4):
    print('=' * 100)
    print('PART A: pair site J2 - axioms and the plus construction from the colimit definition')
    for n in range(0, nmax + 1):
        site = Site(n, covers_J2, 'J2')
        ax = site.check_axioms()
        ncov_top = len(site.cov[(1 << n) - 1])
        irr = [S for S in site.objs if site.irreducible(S)]
        emp = [S for S in site.objs if site.empty_covered(S)]
        print(f'n={n}: axiom violations {dict(ax) or "none"}; covering sieves on V: {ncov_top}; '
              f'irreducible objects = subsets of size <= 2: {sorted(irr) == sorted(S for S in site.objs if popcount(S) <= 2)}; '
              f'objects covered by the empty sieve: {len(emp)}')
        objs = site.objs
        # linear orders
        L = linear_orders(objs)
        a = analyse(site, L, 'L')
        ok_counts = all(a['per_S'][S]['Fp'] == 2 ** (popcount(S) * (popcount(S) - 1) // 2) for S in objs)
        sheaf_iff_small = all((a['per_S'][S]['inj'] and a['per_S'][S]['sur']) == (popcount(S) <= 2) for S in objs)
        # decode each class of L^+(S) as a tournament through its unique R_min = R_2(S) family (Lemma 0);
        # (i) eta_L(order) = comparability tournament of the order, image = transitive tournaments (|S|! of them)
        # (ii) A3: a class lies in the image iff its restriction to every triple does
        Pp = a['Pp']
        img_ok = True
        triple_ok = True
        trans_count = {}

        def arcs_of(f):
            return frozenset((x[0], x[1]) for U, x in f if popcount(U) == 2)

        def transitive(arcs, e):
            return not any((x, y) in arcs and (y, z) in arcs and (z, x) in arcs
                           for x, y, z in itertools.permutations(e, 3))
        image = {S: a['per_S'][S]['img'] for S in objs}
        for S in objs:
            e = elems_of(S)
            for order in L[0][S]:
                comp = frozenset((order[i], order[j]) for i in range(len(order)) for j in range(i + 1, len(order)))
                if arcs_of(a['decode'](S, a['eta'](S, order))) != comp:
                    img_ok = False
            trans = {c for c in Pp[0][S] if transitive(arcs_of(a['decode'](S, c)), e)}
            trans_count[popcount(S)] = len(trans)
            if trans != image[S] or len(trans) != len(L[0][S]):
                img_ok = False
            for c in Pp[0][S]:
                local = all(Pp[1](S, c, U) in image[U] for U in submasks(S) if popcount(U) == 3)
                if local != (c in image[S]):
                    triple_ok = False
        print(f'   L: Lemma 0 (colimit over all covering sieves = R_min families): {a["lemma0"]}; '
              f'|L^+(S)| = 2^C(|S|,2) for all S: {ok_counts}; L separated: {a["separated"]}; '
              f'L sheaf exactly at |S| <= 2: {sheaf_iff_small}; L^+ is a sheaf: {a["plus_sheaf"]}; '
              f'image of L = transitive tournaments (counts {trans_count}): {img_ok}; transitivity local on '
              f'triples: {triple_ok}')
        # (A1) one step suffices for every presheaf: test a range of presheaves incl. non-separated ones
        tests = [('Tour', tournaments_presheaf(objs)), ('const |A|=2', constant((0, 1), objs)),
                 ('const |A|=3', constant((0, 1, 2), objs)), ('partitions', partitions_presheaf(objs)),
                 ('L x Z2', product(L, constant((0, 1), objs), objs)),
                 ('L x Z2 (Z2 forgotten on |S|<3)', forget_below(L, (0, 1), 3, objs)),
                 ('pt x Z3 (forgotten on |S|<3)', forget_below(constant((0,), objs), (0, 1, 2), 3, objs))]
        rng = random.Random(100 + n)
        if n <= 3:
            for t in range(6):
                tests.append((f'random quotient #{t}', random_quotient_presheaf(objs, (0, 1), (0, 1), rng, 3 + t)))
        for name, P in tests:
            b = analyse(site, P, name)
            nonsep = not b['separated']
            print(f'   {name:32s}: Lemma 0 {b["lemma0"]}; F separated {b["separated"]}; F sheaf {b["sheaf"]}; '
                  f'F^+ separated {b["plus_separated"]}; F^+ sheaf (A1) {b["plus_sheaf"]}; restriction maps well '
                  f'defined {b["welldef"]}' + ('  [non-separated input]' if nonsep else ''))


# ----------------------------------------------------------------------------------------------------------
# PART B: the strong-component site Jss
# ----------------------------------------------------------------------------------------------------------
def make_covers_Jss(out):
    cache = {}

    def comps(S):
        if S not in cache:
            cache[S] = components(out, S)
        return cache[S]

    def covers(S, R):
        return all(C in R for C in comps(S))
    return covers, comps


def fibre_count(P, comps, S):
    """|F(C_1) x_F(empty) ... x_F(empty) F(C_c)|"""
    elems, res = P
    cs = comps(S)
    if not cs:
        return 1
    tot = 0
    for e in elems[0]:
        prod = 1
        for C in cs:
            prod *= sum(1 for x in elems[C] if res(C, x, 0) == e)
        tot += prod
    return tot


def second_step_criterion(P, comps, S):
    """F^+(S) -> F^{++}(S) = prod F(C_i) fails to be onto iff c >= 2, all F(C_i) nonempty and the images
    of F(C_i) -> F(empty) are not all one common singleton."""
    elems, res = P
    cs = comps(S)
    if len(cs) < 2 or any(len(elems[C]) == 0 for C in cs):
        return False
    imgs = [frozenset(res(C, x, 0) for x in elems[C]) for C in cs]
    return not (all(len(I) == 1 for I in imgs) and len(set(imgs)) == 1)


def part_B(nmax=4, n5_samples=4):
    print('=' * 100)
    print('PART B: strong-component site Jss - axioms and the plus construction from the colimit definition')
    rep = Counter()
    t_start = time.time()
    cases = []
    for n in range(1, nmax + 1):
        for bits in range(1 << (n * (n - 1) // 2)):
            cases.append((n, tour_from_bits(n, bits), True))
    rng = random.Random(7)
    for _ in range(n5_samples):
        cases.append((5, random_tour(5, rng), False))
    for n, out, full in cases:
        covers, comps = make_covers_Jss(out)
        site = Site(n, covers, 'Jss')
        ax = site.check_axioms(full_transitivity=full, rng=rng)
        rep['axiom violations'] += sum(ax.values())
        rep['tournaments'] += 1
        # points / irreducible objects
        for S in site.objs:
            strong = (S != 0 and len(comps(S)) == 1)
            rep['irreducible <=> nonempty strong: ' + str(site.irreducible(S) == strong)] += 1
            rep['empty covered by empty sieve only at S=0: ' + str(site.empty_covered(S) == (S == 0))] += 1
            rep['R_min = down-closure of the components: ' + str(site.rmin(S) == frozenset(
                U for C in comps(S) for U in submasks(C)))] += 1
        objs = site.objs
        presheaves = [('const |A|=2', constant((0, 1), objs)), ('const |A|=3', constant((0, 1, 2), objs)),
                      ('Pk2', packings_presheaf(out, n, objs)), ('Q', q_presheaf(out, n, objs)),
                      ('Tour', tournaments_presheaf(objs)) if n <= 4 else None,
                      ('L', linear_orders(objs)),
                      ('Pk2 x Z2', product(packings_presheaf(out, n, objs), constant((0, 1), objs), objs))]
        if n <= 3:
            for t in range(2):
                presheaves.append((f'random quotient #{t}',
                                   random_quotient_presheaf(objs, (0, 1), (0, 1), rng, 2 + 2 * t)))
        for item in presheaves:
            if item is None:
                continue
            name, P = item
            b = analyse(site, P, name)
            elems, res = P
            rep['Lemma 0 holds: ' + str(b['lemma0'])] += 1
            rep['input presheaf separated: ' + str(b['separated'])] += 1
            rep['restriction of F^+ well defined: ' + str(b['welldef'])] += 1
            for S in objs:
                ps = b['per_S'][S]
                cs = comps(S)
                # (B2) plus formula
                rep['(B2) |F^+(S)| = fibre product: ' + str(ps['Fp'] == fibre_count(P, comps, S))] += 1
                # sheafification in general: a(F)(S) = prod F(C_i)
                prodF = 1
                for C in cs:
                    prodF *= len(elems[C])
                rep['a(F)(S)=F^{++}(S) has size prod |F(C_i)|: ' + str(ps['Fpp'] == prodF)] += 1
                # second step changes F^+(S) iff criterion
                rep['(B4) second step acts exactly per criterion: ' +
                    str((ps['Fp'] != ps['Fpp'] or not ps['sur2']) == second_step_criterion(P, comps, S))] += 1
            # (B3) sheaf criterion
            mult = (len(elems[0]) == 1) and all(
                len(elems[S]) == fibre_count(P, comps, S) and
                len(set(tuple(res(S, x, C) for C in comps(S)) for x in elems[S])) == len(elems[S])
                for S in objs)
            rep['(B3) sheaf <=> F(empty)=pt and multiplicative: ' + str(b['sheaf'] == mult)] += 1
            # (B4) F(empty) = pt => F^+ is a sheaf
            if len(elems[0]) == 1:
                rep['(B4) F(empty)=pt => F^+ sheaf: ' + str(b['plus_sheaf'])] += 1
            rep['F^+ separated (Grothendieck): ' + str(b['plus_separated'])] += 1
            rep['F^{++} computed; F^+ sheaf when F separated: ' +
                str((not b['separated']) or b['plus_sheaf'])] += 1
            # constant presheaf formulas
            if name.startswith('const'):
                k = len(elems[0])
                okc = all(b['per_S'][S]['Fp'] == (k if S else 1) and b['per_S'][S]['Fpp'] == k ** len(comps(S))
                          for S in objs)
                rep[f'(B4) constant presheaf A^+ = A, A^++ = A^c ({name}): ' + str(okc)] += 1
            if name in ('Pk2',):
                rep['(C1) Pk2 is a Jss-sheaf: ' + str(b['sheaf'])] += 1
                H = ham_count(out, (1 << n) - 1)
                rep['(C1) |Pk2(V)| = H(T): ' + str(len(elems[(1 << n) - 1]) == H)] += 1
            if name == 'Q':
                rep['(D2) Q^+ = a(Q) (one step): ' + str(b['plus_sheaf'])] += 1
    for k in sorted(rep):
        print(f'   {k}: {rep[k]}')
    print(f'   [{len(cases)} labelled tournaments: all with n <= {nmax}, plus {n5_samples} random n = 5; '
          f'{time.time() - t_start:.1f}s]')


# ----------------------------------------------------------------------------------------------------------
# PART B': constant presheaf, one and two plus steps, on ALL subsets of random tournaments n <= 6,
# using the smallest covering sieve (Lemma 0, verified by the colimit in PART B for n <= 5).
# ----------------------------------------------------------------------------------------------------------
def plus_rmin(objs, rmin, P):
    elems, res = P
    fams = {S: matching_families(rmin[S], P) for S in objs}
    idx = {S: {f: i for i, f in enumerate(fams[S])} for S in objs}

    def pres(S, i, Sp):
        f = fams[S][i]
        g = tuple(p for p in f if p[0] in rmin[Sp])
        return idx[Sp][g]
    return ({S: list(range(len(fams[S]))) for S in objs}, pres)


def part_Bprime(trials=60, seed=11):
    print('=' * 100)
    print("PART B': constant presheaf A (|A| = 2, 3; A(empty) = A), one and two plus steps on all subsets")
    rng = random.Random(seed)
    rep = Counter()
    for t in range(trials):
        n = rng.randint(3, 6)
        out = random_tour(n, rng)
        objs = list(range(1 << n))
        comps = {S: components(out, S) for S in objs}
        rmin = {S: frozenset(U for C in comps[S] for U in submasks(C)) for S in objs}
        for A in ((0, 1), (0, 1, 2)):
            P = constant(A, objs)
            P1 = plus_rmin(objs, rmin, P)
            P2 = plus_rmin(objs, rmin, P1)
            for S in objs:
                c = len(comps[S])
                n1, n2 = len(P1[0][S]), len(P2[0][S])
                rep['A^+(S) = A (S nonempty), pt (S empty): ' + str(n1 == (len(A) if S else 1))] += 1
                rep['A^++(S) = A^c(T[S]): ' + str(n2 == len(A) ** c)] += 1
                if c >= 2:
                    rep['subsets with c >= 2 (both steps needed)'] += 1
                    rep['   ... and A^+(S) != A^++(S): ' + str(n1 != n2)] += 1
                rep['subset-checks'] += 1
    for k in sorted(rep):
        print(f'   {k}: {rep[k]}')


if __name__ == '__main__':
    part_A(4)
    part_B(4, n5_samples=3)
    part_Bprime()
    print(f'total time {time.time() - T0:.1f}s')
