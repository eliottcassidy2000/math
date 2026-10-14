"""Exact finite information, action-marking and cocycle-boundary controls.

No statistical model is used to assert an unknown integer reaches ROOT.
All public receipts validate exact types; all checks remain active under -O.
"""
from collections import defaultdict
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def rational(x):
    need(type(x) in (int, F), "exact rational required")
    return F(x)


def statistic(labels, size):
    need(type(labels) is tuple and len(labels) == size
         and all(type(x) is int and x >= 0 for x in labels), "exact statistic labels")


@dataclass(frozen=True)
class Fisher:
    full: F
    compressed: F
    loss: F
    score: tuple
    conditional_score: tuple


def fisher(probabilities, derivatives, labels):
    """One fixed parameter: positive finite law, normalized tangent, fixed statistic."""
    need(type(probabilities) is tuple and len(probabilities) > 0, "nonempty finite law")
    need(type(derivatives) is tuple and len(derivatives) == len(probabilities), "same tangent source")
    p, dp = tuple(map(rational, probabilities)), tuple(map(rational, derivatives))
    need(all(x > 0 for x in p) and sum(p) == 1 and sum(dp) == 0, "positive normalized law and tangent")
    statistic(labels, len(p))
    mass, change = defaultdict(F), defaultdict(F)
    for x, d, t in zip(p, dp, labels):
        mass[t] += x
        change[t] += d
    score = tuple(d/x for x, d in zip(p, dp))
    mean = tuple(change[t]/mass[t] for t in labels)
    full = sum((x*s*s for x, s in zip(p, score)), F(0))
    compressed = sum((change[t]**2/mass[t] for t in mass), F(0))
    loss = sum((x*(s-m)**2 for x,s,m in zip(p,score,mean)), F(0))
    need(full-compressed == loss, "conditional-score identity")
    return Fisher(full, compressed, loss, score, mean)


def verdict_tangent(probabilities, verdicts):
    p = tuple(map(rational, probabilities))
    need(type(verdicts) is tuple and len(verdicts) == len(p)
         and all(type(v) is bool for v in verdicts), "literal Boolean verdicts")
    need(all(x > 0 for x in p) and sum(p) == 1, "positive normalized law")
    mean = sum((x*int(v) for x,v in zip(p,verdicts)), F(0))
    return tuple(x*(int(v)-mean) for x,v in zip(p,verdicts))


def descends(verdicts, labels):
    need(type(verdicts) is tuple and all(type(v) is bool for v in verdicts), "Boolean verdicts")
    statistic(labels, len(verdicts))
    fibres = defaultdict(set)
    for v,t in zip(verdicts,labels):
        fibres[t].add(v)
    return all(len(f) == 1 for f in fibres.values())


def partitions(n):
    if n == 0:
        return ((),)
    rows = [(0,)]
    for _ in range(n-1):
        rows = [row+(x,) for row in rows for x in range(max(row)+2)]
    return tuple(rows)


def word_check(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word), "positive exact valuation word")


def carrier(word):
    word_check(word)
    p,q,b = 1,1,0
    for a in word:
        p,q,b = 3*p,q*2**a,3*b+q
    return p,q,b


def native_residue(word):
    p,q,b = carrier(word)
    return ((q-b)*pow(p,-1,2*q)) % (2*q), 2*q


def replay(source, word, require_root=False):
    need(type(source) is int and source >= 1 and source % 2 == 1, "positive odd exact source")
    need(type(require_root) is bool, "Boolean endpoint requirement")
    word_check(word)
    x = source
    for a in word:
        need(x != 1, "ROOT padding is forbidden")
        t = 3*x+1
        actual = (t & -t).bit_length()-1
        need(a == actual, "wrong original-source valuation")
        x = t >> a
    need(not require_root or x == 1, "supplied word does not reach ROOT")
    return x


@dataclass(frozen=True)
class RootReceipt:
    source: int
    word: tuple


def audit_root(receipt):
    need(type(receipt) is RootReceipt, "source-owned ROOT receipt")
    return replay(receipt.source, receipt.word, True) == 1


def prepend(source, prefix, child):
    audit_root(child)
    endpoint = replay(source, prefix)
    need(endpoint == child.source, "child certificate source mismatch")
    out = RootReceipt(source, prefix+child.word)
    audit_root(out)
    return out


def perfect_matchings(vertices):
    if not vertices:
        return ((),)
    a = vertices[0]
    return tuple(tuple(sorted(((a,b),)+tail)) for b in vertices[1:]
                 for tail in perfect_matchings(tuple(v for v in vertices if v not in (a,b))))


MATCHINGS = tuple(sorted(perfect_matchings(tuple(range(6)))))
MATCHING_INDEX = {m:i for i,m in enumerate(MATCHINGS)}
PENTADS = tuple(c for c in combinations(range(15),5)
                if len({edge for i in c for edge in MATCHINGS[i]}) == 15)
PENTAD_INDEX = {p:i for i,p in enumerate(PENTADS)}


def permutation(p):
    need(type(p) is tuple and len(p) == 6 and all(type(x) is int for x in p)
         and set(p) == set(range(6)), "permutation of six marked points")


def compose(p,q):
    permutation(p)
    permutation(q)
    return tuple(p[q[i]] for i in range(6))


def outer(p):
    permutation(p)
    on_matchings = tuple(MATCHING_INDEX[tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in m))]
                        for m in MATCHINGS)
    return tuple(PENTAD_INDEX[tuple(sorted(on_matchings[j] for j in pentad))] for pentad in PENTADS)


@dataclass(frozen=True)
class ActionReceipt:
    point_action: tuple
    pentad_action: tuple


def audit_action(receipt):
    need(type(receipt) is ActionReceipt, "marked two-action receipt")
    permutation(receipt.point_action)
    permutation(receipt.pentad_action)
    need(receipt.pentad_action == outer(receipt.point_action), "wrong representation transport")
    return True


def cycle_type(p):
    permutation(p)
    remaining = set(range(6))
    lengths = []
    while remaining:
        start = min(remaining)
        x, length = start,0
        while x in remaining:
            remaining.remove(x)
            length += 1
            x = p[x]
        lengths.append(length)
    return tuple(sorted(lengths,reverse=True))


def doubling_cycles(modulus):
    need(type(modulus) is int and modulus > 1 and modulus % 2 == 1, "odd exact modulus")
    from math import gcd
    remaining = {a for a in range(modulus) if gcd(a,modulus) == 1}
    cycles = []
    while remaining:
        a = min(remaining)
        row = []
        while a in remaining:
            row.append(a)
            remaining.remove(a)
            a = 2*a % modulus
        cycles.append(tuple(row))
    return tuple(cycles)


def cycle_cocycle(charges, modulus, killed=None):
    """n-cycle i->i+1; normalized exponent at vertex0, or at killed vertex."""
    need(type(charges) is tuple and len(charges) > 0
         and all(type(q) is int for q in charges), "exact charges on a nonempty directed cycle")
    need(type(modulus) is int and modulus >= 1, "positive exact character modulus")
    n = len(charges)
    need(killed is None or type(killed) is int and 0 <= killed < n, "typed killed terminal")
    if killed is None and sum(charges) % modulus:
        return None
    anchor = 0 if killed is None else killed
    values = [None]*n
    values[anchor] = 0
    for k in range(1,n):
        i = (anchor-k) % n
        values[i] = (charges[i]+values[(i+1)%n]) % modulus
    need(all(i == killed or (values[i]-charges[i]-values[(i+1)%n]) % modulus == 0
             for i in range(n)), "literal retained-edge cocycle")
    return tuple(values)


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)
    def reject(fn,*args):
        try:
            fn(*args)
        except ValueError:
            check(True,"hostile rejected")
        else:
            check(False,"invalid packet accepted")

    # All Boolean verifiers and all deterministic partitions on 1..5 points.
    verdict_controls = 0
    for n in range(1,6):
        p = tuple(F(i+1,n*(n+1)//2) for i in range(n))
        for t in partitions(n):
            for v in product((False,True),repeat=n):
                f = fisher(p,verdict_tangent(p,v),t)
                check((f.loss == 0) == descends(v,t), "verdict-specific Fisher equality iff verifier descends")
                verdict_controls += 1

    # Truncated geometric valuation family: all order information is ancillary.
    words = tuple(product((1,2,3),repeat=2))
    z = F(1,2)
    normalizer = 1+z+z*z
    law = tuple(z**(sum(w)-2)/normalizer**2 for w in words)
    tangent = tuple(p*((sum(w)-2)/z-2*(1+2*z)/normalizer) for p,w in zip(law,words))
    costs = tuple(sum(w) for w in words)
    geometric = fisher(law,tangent,costs)
    check(geometric.full == geometric.compressed == F(208,49), "exact geometric-model sufficiency")
    v = tuple(w == (1,3) for w in words)
    verdict = fisher(law,verdict_tangent(law,v),costs)
    check(verdict.loss == F(8,147) and not descends(v,costs), "same compression loses actual source verdict")
    check(native_residue((1,3)) == (19,32) and native_residue((3,1)) == (29,32), "same slope different source cells")
    check(carrier((1,3)) == (9,16,5) and carrier((3,1)) == (9,16,11), "affine carry retained")
    for w in words:
        try:
            replay(19,w)
            actual = True
        except ValueError:
            actual = False
        check(actual == (w == (1,3)), "literal original-source verifier")
    root = prepend(19,(1,3),RootReceipt(11,(1,2,3,4)))
    check(root == RootReceipt(19,(1,3,1,2,3,4)) and audit_root(root), "actual strict source-owned completion")
    reject(prepend,19,(1,3),RootReceipt(17,(2,3,4)))
    reject(replay,19,(3,1))
    reject(audit_root,RootReceipt(1,(2,)))
    reject(audit_root,RootReceipt(True,()))

    # A sufficient score need not retain an ancillary Boolean certificate.
    p = (F(1,8),F(1,8),F(3,8),F(3,8))
    dp = (F(-1,2),F(-1,2),F(1,2),F(1,2))
    f = fisher(p,dp,(0,0,1,1))
    check(f.loss == 0 and not descends((False,True,False,True),(0,0,1,1)), "ancillary verdict hostile")

    # Actual six-point / six-pentad action, not a cosmetic relabelling.
    check(len(MATCHINGS) == 15 and len(PENTADS) == 6, "complete matching-factorization universe")
    perms = tuple(permutations(range(6)))
    images = {p:outer(p) for p in perms}
    check(len(set(images.values())) == 720, "outer map is a bijection")
    generators = []
    for i in range(5):
        g = list(range(6)); g[i],g[i+1] = g[i+1],g[i]
        generators.append(tuple(g))
    for p in perms:
        check(audit_action(ActionReceipt(p,images[p])), "typed action transport")
        for g in generators:
            check(images[compose(p,g)] == compose(images[p],images[g]), "homomorphism on generating edges")
    swap = generators[0]
    check(cycle_type(swap) == (2,1,1,1,1) and cycle_type(images[swap]) == (2,2,2), "pair swap becomes matching swap")
    check(sum(swap[i] == i for i in range(6)) == 4
          and all(images[swap][i] != i for i in range(6)), "no pointwise intertwining bijection")
    p = tuple(F(i+1,720*721//2) for i in range(720))
    mean = sum((p[i]*int(perms[i][0] == 0) for i in range(720)),F(0))
    dp = tuple(p[i]*(int(perms[i][0] == 0)-mean) for i in range(720))
    labels_outer = tuple(perms.index(images[g]) for g in perms)
    check(fisher(p,dp,labels_outer).loss == 0, "invertible action recoding retains all Fisher information")
    reject(audit_action,ActionReceipt(swap,swap))
    reject(outer,(True,0,2,3,4,5))

    # THM-4631: cycle index discards the phase coordinate, even with one orbit.
    orbit_reports = []
    for modulus in (3,7,9,11,17,23):
        cyc = doubling_cycles(modulus)
        units = tuple(sorted(a for row in cyc for a in row))
        t = tuple(next(i for i,row in enumerate(cyc) if a in row) for a in units)
        p = (F(1,len(units)),)*len(units)
        dp = tuple(F(1,len(units)) if a == cyc[0][0] else
                   F(-1,len(units)) if a == cyc[0][1] else F(0) for a in units)
        f = fisher(p,dp,t)
        check(f.compressed == 0 and f.full == F(2,len(units)), "orbit label loses phase score")
        phase = tuple(next(j for row in cyc for j,x in enumerate(row) if a == x) for a in units)
        marked = tuple(t[i]*len(units)+phase[i] for i in range(len(units)))
        check(fisher(p,dp,marked).loss == 0, "cycle label plus phase is lossless")
        orbit_reports.append((modulus,len(cyc),len(cyc[0]),str(f.full)))
    check(doubling_cycles(7) == ((1,2,4),(3,6,5)), "Paley residue split with marked phases")

    # Exact holonomy, including the retained trivial Collatz cycle.
    cycles_checked = 0
    for n in range(1,8):
        for q in product((0,1),repeat=n):
            for m in range(1,8):
                c = cycle_cocycle(q,m)
                check((c is not None) == (sum(q)%m == 0), "cycle holonomy iff")
                cycles_checked += 1
    check(cycle_cocycle((1,0),2) is None, "1->2->1 forbids nonzero odd-sign cocycle")
    check(cycle_cocycle((1,0),2,killed=0) == (0,0), "ROOT-killed normalization succeeds")
    check(cycle_cocycle((1,0),1) == (0,0), "trivial character positive control")

    # Exact approach-ray example: rational continuous bounded f, not an orbit theorem.
    rays = []
    for slope in (F(0),F(1),F(2),F(-1)):
        expected = slope*slope/(1+slope*slope)
        for denominator in range(1,21):
            y=F(1,denominator); x=slope*y
            check(x*x/(x*x+y*y) == expected, "same boundary point different exact approach value")
        rays.append((str(slope),str(expected)))

    # Certainty about a state is not an assertion about its transition law.
    check(fisher((F(1),),(F(0),),(0,)).full == 0, "constant Dirac family has zero Fisher")
    root_map=(0,0); loop_map=(0,1)
    check(root_map[1] == 0 and loop_map[1] == 1, "same certain start, different terminal event")
    reject(fisher,(F(1),),(F(0),),(False,))
    reject(fisher,(F(0),F(1)),(F(0),F(0)),(0,1))
    reject(fisher,(0.5,0.5),(0,0),(0,1))
    reject(verdict_tangent,(F(1),),(1,))
    reject(cycle_cocycle,(1,0),2,True)

    report={"checks":checks,"status":"PROVED finite observer identities; FINITE-EXACT controls; no universal ROOT inference",
        "all_verdict_partition_controls":verdict_controls,
        "collatz_model":{"geometric_fisher":str(geometric.full),"cost_fisher_loss":str(geometric.loss),
            "source_verdict_full":str(verdict.full),"source_verdict_compressed":str(verdict.compressed),
            "source_verdict_loss":str(verdict.loss),"source_cells":[[19,32],[29,32]],
            "strict_root_source":root.source,"strict_root_word":root.word},
        "outer_s6":{"matchings":len(MATCHINGS),"pentads":len(PENTADS),"actions":len(images),
            "point_transposition":swap,"pentad_action":images[swap],"pentad_cycle_type":cycle_type(images[swap])},
        "doubling_cycle_controls":orbit_reports,"holonomy_controls":cycles_checked,
        "boundary_ray_values":rays}
    print(json.dumps(report,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
