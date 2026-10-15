"""Exact, source-bound transport of fixed gold comparison certificates.

All-order proofs are in the companion note. The exhaustive control universe is
naturally labelled posets of size 2..4, not a census proving the conjecture.
Run python -B [-O] -X utf8 this_file.py; stdout is deterministic JSON.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import comb
import json


def integer(x, low=0):
    if type(x) is not int or x < low:
        raise ValueError("exact integer outside domain")
    return x


@dataclass(frozen=True)
class Poset:
    pred: tuple

    @property
    def n(self):
        return len(self.pred)


def poset(n, edges=()):
    integer(n)
    p = [0] * n
    for e in edges:
        if type(e) is not tuple or len(e) != 2:
            raise ValueError("edge pair")
        a, b = e
        integer(a); integer(b)
        if a >= n or b >= n or a == b:
            raise ValueError("edge endpoint")
        p[b] |= 1 << a
    for k in range(n):
        for j in range(n):
            if p[j] >> k & 1:
                p[j] |= p[k]
    if any(p[j] >> j & 1 for j in range(n)):
        raise ValueError("cycle")
    return Poset(tuple(p))


def audit_poset(p):
    if type(p) is not Poset or type(p.pred) is not tuple:
        raise ValueError("poset type")
    for x in p.pred:
        integer(x)
        if x >= 1 << p.n:
            raise ValueError("predecessor range")
    q = poset(p.n, tuple((i, j) for j in range(p.n)
                        for i in range(p.n) if p.pred[j] >> i & 1))
    if p != q:
        raise ValueError("not transitive closure")
    return p


def pair(p, e):
    if type(e) is not tuple or len(e) != 2:
        raise ValueError("comparison pair")
    a, b = e
    integer(a); integer(b)
    if a >= p.n or b >= p.n or a == b:
        raise ValueError("comparison source")
    return e


def edges(p):
    return tuple((i, j) for j in range(p.n)
                 for i in range(p.n) if p.pred[j] >> i & 1)


@lru_cache(None)
def _count(pred):
    full = (1 << len(pred)) - 1
    @lru_cache(None)
    def go(done):
        if done == full:
            return 1
        return sum(go(done | 1 << i) for i, need in enumerate(pred)
                   if not done >> i & 1 and need & done == need)
    return go(0)


def count(p, conditions=()):
    audit_poset(p)
    if type(conditions) is not tuple:
        raise ValueError("conditions tuple")
    for c in conditions:
        pair(p, c)
    try:
        q = poset(p.n, edges(p) + conditions)
    except ValueError:
        return 0
    return _count(q.pred)


def extensions(p):
    audit_poset(p)
    if p.n > 9:
        raise ValueError("literal base enumeration cap 9")
    def go(word, done):
        if len(word) == p.n:
            yield word
        for i, need in enumerate(p.pred):
            if not done >> i & 1 and need & done == need:
                yield from go(word + (i,), done | 1 << i)
    return tuple(go((), 0))


def cells(p, first, second):
    audit_poset(p); pair(p, first); pair(p, second)
    return tuple(count(p, (first if i == 0 else first[::-1],
                           second if j == 0 else second[::-1]))
                 for i in range(2) for j in range(2))


def deficits(c):
    return tuple(c[2*(1-i)] + c[2*(1-i)+1] - c[2*i+j]
                 for i in range(2) for j in range(2))


@dataclass(frozen=True)
class Gold:
    source: Poset
    first: tuple
    second: tuple
    counts: tuple


def compile_gold(p, first, second):
    c = cells(p, first, second)
    if min(deficits(c)) < 0:
        raise ValueError("not a gold certificate")
    return Gold(p, first, second, c)


def audit_gold(g):
    if type(g) is not Gold or type(g.counts) is not tuple or len(g.counts) != 4:
        raise ValueError("certificate type")
    for v in g.counts:
        integer(v)
    if g != compile_gold(g.source, g.first, g.second):
        raise ValueError("unauthenticated counts")
    return g


def robustness(g):
    audit_gold(g)
    c = g.counts
    return min(Fraction(c[2*(1-i)] + c[2*(1-i)+1], c[2*i+j])
               for i in range(2) for j in range(2) if c[2*i+j])


def weighted_cells(g, weights):
    audit_gold(g)
    words = extensions(g.source)
    if type(weights) is not tuple or len(weights) != len(words):
        raise ValueError("one weight per marked base extension")
    for w in weights:
        integer(w, 1)
    out = [0] * 4
    for w, word in zip(weights, words):
        loc = {a: i for i, a in enumerate(word)}
        i = int(loc[g.first[0]] > loc[g.first[1]])
        j = int(loc[g.second[0]] > loc[g.second[1]])
        out[2*i+j] += w
    return tuple(out)


def ideal(p, D):
    audit_poset(p)
    if type(D) is not tuple or len(set(D)) != len(D):
        raise ValueError("ideal tuple")
    for i in D:
        integer(i)
        if i >= p.n:
            raise ValueError("ideal source")
    if D != tuple(sorted(D)):
        raise ValueError("canonical ideal order")
    mask = sum(1 << i for i in D)
    if any(p.pred[i] & mask != p.pred[i] for i in D):
        raise ValueError("not an order ideal")
    return D


@dataclass(frozen=True)
class InsertionPlan:
    certificate: Gold
    lower_ideal: tuple
    slopes: tuple
    intercepts: tuple
    minimum_buffer: object
    obstruction: object


def insertion_plan(g, D):
    audit_gold(g); ideal(g.source, D)
    n = g.source.n
    ell = tuple(n - max((word.index(a) + 1 for a in D), default=0) + 1
                for word in extensions(g.source))
    delta, gamma = deficits(g.counts), deficits(weighted_cells(g, ell))
    bad = next((i for i in range(4) if delta[i] == 0 and gamma[i] < 0), None)
    minimum = None if bad is not None else max(
        [0] + [(-gamma[i] + delta[i] - 1)//delta[i]
               for i in range(4) if delta[i] > 0])
    return InsertionPlan(g, D, delta, gamma, minimum, bad)


def audit_plan(p):
    if type(p) is not InsertionPlan:
        raise ValueError("plan type")
    for row in (p.slopes, p.intercepts):
        if type(row) is not tuple or len(row) != 4 or any(type(x) is not int for x in row):
            raise ValueError("exact deficit tuple")
    if p.minimum_buffer is not None:
        integer(p.minimum_buffer)
    if p.obstruction is not None:
        integer(p.obstruction)
    if p != insertion_plan(p.certificate, p.lower_ideal):
        raise ValueError("plan source mismatch")
    return p


def accepts(p, m):
    audit_plan(p); integer(m)
    return p.minimum_buffer is not None and m >= p.minimum_buffer


def insertion_source(p, m):
    audit_plan(p); integer(m)
    q = p.certificate.source
    z = q.n + m
    es = edges(q) + tuple((a, z) for a in p.lower_ideal)
    es += tuple((a, q.n) for a in range(q.n)) if m else ()
    es += tuple((q.n+i, q.n+i+1) for i in range(m-1))
    return poset(z+1, es)


def core_of(p, n):
    audit_poset(p); integer(n)
    if n > p.n:
        raise ValueError("core size")
    ideal(p, tuple(range(n)))
    return Poset(tuple(x for x in p.pred[:n]))


def buffered_source(p, n, m):
    core_of(p, n); integer(m)
    es = edges(p)
    if m:
        es += tuple((i, p.n) for i in range(n))
        es += tuple((p.n+i, p.n+i+1) for i in range(m-1))
    return poset(p.n+m, es)


def fibre_ratio(n, r, m):
    integer(n); integer(r); integer(m)
    return Fraction(comb(m+n+r, r), comb(m+r, r))


def context_buffer(g, p):
    audit_gold(g)
    if core_of(p, g.source.n) != g.source:
        raise ValueError("core source differs")
    n, r, R = g.source.n, p.n-g.source.n, robustness(g)
    if r == 0:
        return 0
    if R <= 1:
        raise ValueError("strict certificate required for a uniform buffer bound")
    lo, hi = -1, 0
    while fibre_ratio(n, r, hi) > R:
        hi = 2*hi+1
    while hi-lo > 1:
        mid = (hi+lo)//2
        if fibre_ratio(n, r, mid) <= R:
            hi = mid
        else:
            lo = mid
    return hi


def context_schedule(g, p):
    """Exact all-m classification, using a finite source extension enumeration.

    Returns a certified sufficient tail, the complete exceptional invalid set,
    and fibre coefficients in the binomial basis binom(m+h,h). No expanded
    buffered source is enumerated here. Literal source cap is nine vertices.
    """
    bound = context_buffer(g,p)
    n, r = g.source.n, p.n-g.source.n
    words = extensions(g.source)
    index = {word:i for i,word in enumerate(words)}
    hist = [[0]*(r+1) for _ in words]
    for word in extensions(p):
        restricted = tuple(a for a in word if a < n)
        h = p.n-1-max(i for i,a in enumerate(word) if a < n)
        hist[index[restricted]][h] += 1
    hist = tuple(tuple(row) for row in hist)
    invalid = []
    for m in range(bound):
        weights = tuple(sum(c*comb(m+h,h) for h,c in enumerate(row)) for row in hist)
        if min(deficits(weighted_cells(g,weights))) < 0:
            invalid.append(m)
    return bound, tuple(invalid), hist


def natural_posets(n):
    es = tuple(combinations(range(n), 2))
    return tuple(sorted({poset(n, tuple(e for i,e in enumerate(es) if mask >> i & 1)).pred
                         for mask in range(1 << len(es))}))


def verify():
    checks = 0
    def check(v):
        nonlocal checks
        checks += 1
        if not v:
            raise AssertionError(checks)
    def hostile(f):
        try:
            f()
        except ValueError:
            check(True)
        else:
            check(False)
    a3 = poset(3)
    g = compile_gold(a3, (0,1), (0,2))
    ip = insertion_plan(g, (0,))
    check(g.counts == (2,1,1,2) and robustness(g) == Fraction(3,2))
    check(ip.minimum_buffer == 2)
    check(ip.slopes == (1,2,2,1) and ip.intercepts == (-2,2,6,6))
    for m in range(9):
        c = cells(insertion_source(ip,m), g.first, g.second)
        check(deficits(c) == tuple(m*d+b for d,b in zip(ip.slopes,ip.intercepts)))
        check((min(deficits(c)) >= 0) == accepts(ip,m))
    weak = compile_gold(poset(2), (0,1), (0,1))
    wp = insertion_plan(weak,(0,))
    check(wp.minimum_buffer is None and wp.intercepts[0] == -1)
    for m in range(17):
        check(min(deficits(cells(insertion_source(wp,m),weak.first,weak.second))) == -1)
    # Weighted sharpness and exact finite laws, independently assigning all 3^6 weights.
    weighted = 0
    for ws in product(range(1,4), repeat=6):
        c = weighted_cells(g,ws)
        check(sum(c) == sum(ws))
        if Fraction(max(ws),min(ws)) <= robustness(g):
            check(min(deficits(c)) >= 0)
        weighted += 1
    ws = (2,2,1,1,1,1)  # first two extensions have 0 first
    check(min(deficits(weighted_cells(g,ws))) < 0)
    # Complete transport identities, NOT a conjecture census.
    stats = dict(posets=0, gold_certificates=0, insertion_plans=0, literal_buffers=0)
    for n in range(2,5):
        ps = natural_posets(n)
        pairs = tuple(combinations(range(n),2))
        for pred in ps:
            p = Poset(pred); stats['posets'] += 1
            check(count(p) == len(extensions(p)))
            ideals = []
            for mask in range(1 << n):
                D = tuple(i for i in range(n) if mask >> i & 1)
                try:
                    ideal(p,D); ideals.append(D)
                except ValueError:
                    pass
            for e,f in product(pairs, repeat=2):
                try:
                    cert = compile_gold(p,e,f)
                except ValueError:
                    continue
                stats['gold_certificates'] += 1
                for D in ideals:
                    plan = insertion_plan(cert,D); stats['insertion_plans'] += 1
                    trials = {0,1,2,3}
                    if plan.minimum_buffer is not None:
                        trials |= {plan.minimum_buffer, max(0,plan.minimum_buffer-1)}
                    for m in sorted(trials):
                        c = cells(insertion_source(plan,m),e,f)
                        check(deficits(c) == tuple(m*d+b for d,b in zip(plan.slopes,plan.intercepts)))
                        check((min(deficits(c)) >= 0) == accepts(plan,m))
                        stats['literal_buffers'] += 1
    # Arbitrary upper context with three different predecessor ideals; Q is not a module.
    ctx = poset(6, ((0,3),(3,4),(4,5),(1,4),(2,5)))
    mstar = context_buffer(g,ctx)
    bound, invalid, histogram = context_schedule(g,ctx)
    check(bound == mstar)
    check(all(row[-1] == 1 for row in histogram))
    check(fibre_ratio(3,3,mstar) <= robustness(g))
    check(mstar == 0 or fibre_ratio(3,3,mstar-1) > robustness(g))
    context_rows = []
    check(invalid == tuple(range(7)))
    for m in (0,1,2,3,6,7,mstar):
        p = buffered_source(ctx,3,m)
        actual = tuple(count(p, tuple(zip(word,word[1:]))) for word in extensions(a3))
        lower, upper = comb(m+3,3), comb(m+6,3)  # e(R)=1
        check(all(lower <= w <= upper for w in actual))
        check(sum(actual) == count(p))
        check(actual == tuple(sum(c*comb(m+h,h) for h,c in enumerate(row)) for row in histogram))
        c = cells(p,g.first,g.second)
        check(c == weighted_cells(g,actual))
        check(deficits(c)[0] == (m+1)*(m*m-4*m-18)//6)
        check((min(deficits(c)) < 0) == (m in invalid))
        if m == mstar:
            check(min(deficits(c)) >= 0)
        context_rows.append(dict(buffer=m,weights=actual,deficits=deficits(c)))
    # Type/source hostiles before cached counts or dataclass equality.
    hostile(lambda: poset(True))
    hostile(lambda: count(Poset((False,0))))
    hostile(lambda: compile_gold(a3,(False,1),(0,2)))
    hostile(lambda: compile_gold(a3,(0,3),(0,2)))
    hostile(lambda: audit_gold(replace(g, counts=(2.0,1,1,2))))
    hostile(lambda: audit_gold(replace(g, counts=(1,2,1,2))))
    hostile(lambda: insertion_plan(g,(True,)))
    hostile(lambda: insertion_plan(compile_gold(poset(3,((0,1),)),(0,2),(1,2)),(1,)))
    hostile(lambda: accepts(ip,True))
    hostile(lambda: accepts(ip,-1))
    hostile(lambda: audit_plan(replace(ip,minimum_buffer=2.0)))
    hostile(lambda: audit_plan(replace(ip,intercepts=(0,2,6,6))))
    hostile(lambda: context_buffer(g,poset(4,((3,0),))))
    hostile(lambda: context_buffer(weak,poset(3,((0,2),))))
    hostile(lambda: weighted_cells(g,(1,)*5))
    hostile(lambda: weighted_cells(g,(True,)*6))
    return dict(status="FINITE-EXACT; general transport theorems proved in companion note",
                checks=checks, universe=stats, weighted_assignments=weighted,
                strict_example=dict(counts=g.counts,slopes=ip.slopes,intercepts=ip.intercepts,
                                    exact_threshold=ip.minimum_buffer,box_threshold=3),
                weak_hostile=dict(slopes=wp.slopes,intercepts=wp.intercepts,threshold=None),
                nonmodule_context=dict(uniform_bound=mstar,invalid_buffers=invalid,
                                       exact_tail_start=max(invalid,default=-1)+1,
                                       binomial_fibres=histogram,rows=context_rows),
                open="No proof of the gold partition conjecture or of the original unbuffered source.")


if __name__ == '__main__':
    print(json.dumps(verify(),indent=2,sort_keys=True))
