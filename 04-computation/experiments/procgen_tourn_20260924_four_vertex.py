#!/usr/bin/env python3
"""procgen_tourn_20260924, part A: four-vertex tournaments built from the Collatz pair structure.

Question (owner): are the "3" and the "+-1" of 3n+-1 the two 4-vertex tournaments that swap under
reversal of all 6 arcs, (1,1,1,3) = source over a 3-cycle and (0,2,2,2) = 3-cycle over a sink?

Sections (every check raises on failure):
  A1  the 64 labelled 4-tournaments: 4 classes, H, c3, converse action, H = 1 + 2 c3 (OCF at n = 4).
  A2  the AM-fair pairing quadruple Q_b(s_o) = {s_o, s_e = s_o + b, t_o = T_b(s_o), t_e = s_e/2}:
      degeneracy, order patterns, classes of map+order tournaments (forward / inverse-tree arcs,
      ascending / descending order), the reflection rho(x) = sigma - x, the exact identity
      converse = negation o time-reversal, H = 2 + sgn(b s_o).
  A3  exhaustive rule census: 512 local rules (map arcs forward/backward; each of the 4 non-map pairs
      oriented by ascending order, descending order, or a fixed label direction), evaluated on the
      four (sheet, side) patterns.  Which rules give EXACTLY the converse pair across the sheets?
  A4  the four up/down choices on a pair (THM-4470 family + Mahler), other natural quadruples
      (step decomposition n -> 3n -> 3n+b -> (3n+b)/2, the inverse-tree fork, length-4 orbit
      windows, b = 0 on Q, Mahler's ceil/floor(3n/2)), and the involution dictionary per family.
Run: python3 04-computation/experiments/procgen_tourn_20260924_four_vertex.py   (a few seconds, < 50 MB)
"""
import itertools
import sys
from fractions import Fraction as Fr

TT, SRC, SNK, STR = "TT(0,1,2,3)", "SRC+C3(1,1,1,3)", "C3+SNK(0,2,2,2)", "STRONG(1,1,2,2)"
CLASS = {(0, 1, 2, 3): TT, (1, 1, 1, 3): SRC, (0, 2, 2, 2): SNK, (1, 1, 2, 2): STR}
CONV = {TT: TT, STR: STR, SRC: SNK, SNK: SRC}
HVAL = {TT: 1, SRC: 3, SNK: 3, STR: 5}


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


# ---------------------------------------------------------------- tournament utilities
def is_tournament(V, arcs):
    V = list(V)
    if len(set(V)) != len(V):
        return False
    seen = set()
    for (u, v) in arcs:
        if u == v or u not in V or v not in V:
            return False
        key = frozenset((u, v))
        if key in seen:
            return False
        seen.add(key)
    return len(seen) == len(V) * (len(V) - 1) // 2


def scores(V, arcs):
    s = {v: 0 for v in V}
    for (u, _) in arcs:
        s[u] += 1
    return s


def cls(V, arcs):
    check(is_tournament(V, arcs), "not a tournament: %r %r" % (V, sorted(arcs)))
    return CLASS[tuple(sorted(scores(V, arcs).values()))]


def ham(V, arcs):
    A = set(arcs)
    return sum(1 for p in itertools.permutations(V) if all((p[k], p[k + 1]) in A for k in range(len(p) - 1)))


def c3(V, arcs):
    A = set(arcs)
    cnt = 0
    for a, b, c in itertools.combinations(V, 3):
        if ((a, b) in A and (b, c) in A and (c, a) in A) or ((b, a) in A and (c, b) in A and (a, c) in A):
            cnt += 1
    return cnt


def conv(arcs):
    return frozenset((v, u) for (u, v) in arcs)


def relabel(arcs, f):
    return frozenset((f[u], f[v]) for (u, v) in arcs)


# ---------------------------------------------------------------- A1
def section_A1():
    print("=" * 100)
    print("A1  the 64 labelled tournaments on 4 vertices")
    V = [0, 1, 2, 3]
    pairs = list(itertools.combinations(V, 2))
    count = {}
    for bits in range(64):
        arcs = frozenset((u, v) if (bits >> k) & 1 else (v, u) for k, (u, v) in enumerate(pairs))
        c = cls(V, arcs)
        h, t3 = ham(V, arcs), c3(V, arcs)
        check(h == HVAL[c], "H per class")
        check(h == 1 + 2 * t3, "OCF at n=4: H = 1 + 2 c3")
        check(cls(V, conv(arcs)) == CONV[c], "converse action on classes")
        count[c] = count.get(c, 0) + 1
    check(count == {TT: 24, SRC: 8, SNK: 8, STR: 24}, "labelled class counts")
    for c in (TT, SRC, SNK, STR):
        print("  %-17s labelled %2d  H = %d  c3 = %d  converse -> %s" % (c, count[c], HVAL[c], (HVAL[c] - 1) // 2, CONV[c]))
    print("  converse swaps SRC+C3 <-> C3+SNK and fixes TT, STRONG; H = 1 + 2 c3 on all 64: ok")
    print("  a single arc flip changes two scores by +-1, so the labelled distance between an SRC+C3 and a")
    print("  C3+SNK tournament is >= 2 (from (0,2,2,2) one flip reaches only (1,1,2,2) or (0,1,2,3)):")
    dmin = 99
    for b1 in range(64):
        a1 = frozenset((u, v) if (b1 >> k) & 1 else (v, u) for k, (u, v) in enumerate(pairs))
        if cls(V, a1) != SNK:
            continue
        for b2 in range(64):
            a2 = frozenset((u, v) if (b2 >> k) & 1 else (v, u) for k, (u, v) in enumerate(pairs))
            if cls(V, a2) == SRC:
                dmin = min(dmin, bin(b1 ^ b2).count("1"))
    check(dmin == 2, "min labelled distance SNK-SRC is 2")
    print("  min labelled arc-distance between the two diamonds = %d: ok" % dmin)


# ---------------------------------------------------------------- the pairing quadruple
LAB = ("so", "se", "to", "te")
NONMAP = {"A": ("so", "se"), "B": ("to", "te"), "C": ("so", "te"), "D": ("se", "to")}
MAPF = (("so", "to"), ("se", "te"))


def Tb(x, b):
    """Shortcut map on Z (b odd) or on Q with odd denominators (any b)."""
    x = Fr(x)
    num = x.numerator  # parity of a rational with odd denominator = parity of its numerator
    check(x.denominator % 2 == 1, "odd denominator")
    return x / 2 if num % 2 == 0 else (3 * x + b) / 2


def quad(b, so, beta=None):
    """Pair {s_o, s_e} of sheet b: s_e = s_o + beta (beta = b for the AM-fair pairing)."""
    if beta is None:
        beta = b
    se = so + beta
    val = {"so": Fr(so), "se": Fr(se), "to": Tb(so, b), "te": Tb(se, b)}
    return val


def distinct(val):
    return len(set(val.values())) == 4


def build(val, mu="fwd", order="asc", rule=None):
    """Map arcs (forward s->t or backward t->s) + non-map pairs oriented by rule.
    rule: dict pair-name -> 'asc' | 'desc' | 'L0' | 'L1'; default: all = order."""
    arcs = set()
    for (s, t) in MAPF:
        arcs.add((s, t) if mu == "fwd" else (t, s))
    for name, (u, v) in NONMAP.items():
        r = rule[name] if rule else order
        if r == "L0":
            arcs.add((u, v))
        elif r == "L1":
            arcs.add((v, u))
        else:
            up = val[u] < val[v]
            if r == "desc":
                up = not up
            arcs.add((u, v) if up else (v, u))
    return frozenset(arcs)


def pattern(val):
    return tuple(sorted(LAB, key=lambda l: val[l]))


RHO = {"so": "se", "se": "so", "to": "te", "te": "to"}


def section_A2():
    print("=" * 100)
    print("A2  the AM-fair pairing quadruple Q_b = {s_o, s_e = s_o + b, t_o = T_b(s_o), t_e = s_e/2}")
    print("    (b = +1: pairs {2i-1, 2i}; b = -1: pairs {2i, 2i+1}; THM-4470 section 1)")
    R = 2001
    pats = {}
    nondeg = 0
    for b in (1, -1):
        for so in range(-R, R + 1, 2):
            val = quad(b, so)
            se, to, te = val["se"], val["to"], val["te"]
            check(val["to"] + val["te"] == val["so"] + val["se"], "AM-fairness t_o + t_e = s_o + s_e")
            check((val["so"] - se) * (to - te) == -b * so, "inversion identity (s_o-s_e)(t_o-t_e) = -b s_o")
            if not distinct(val):
                check(abs(so) == 1, "degenerate only at |s_o| = 1 (s_o=%d, b=%d)" % (so, b))
                continue
            check(abs(so) != 1, "|s_o| = 1 is degenerate")
            nondeg += 1
            u = b * so
            side = "contracting side (b s_o >= 3)" if u >= 3 else "expanding side (b s_o <= -3)"
            pats.setdefault((b, "pos" if so > 0 else "neg"), set()).add(pattern(val))
            Ffa, Ffd = build(val, "fwd", "asc"), build(val, "fwd", "desc")
            Fba, Fbd = build(val, "bwd", "asc"), build(val, "bwd", "desc")
            V = list(LAB)
            want_f = SNK if u >= 3 else TT
            want_b = SRC if u >= 3 else TT
            check(cls(V, Ffa) == want_f and cls(V, Ffd) == want_f, "Theorem A(a) forward class, %s" % side)
            check(cls(V, Fba) == want_b and cls(V, Fbd) == want_b, "Theorem A(b) backward class")
            # rho is an isomorphism (M, O) -> (M, O^op) and (F+)^op -> F-  (labelled equalities)
            check(relabel(Ffa, RHO) == Ffd, "rho maps forward/asc onto forward/desc")
            check(relabel(conv(Ffa), RHO) == Fba, "rho maps (F+)^op onto F- (time reversal = converse)")
            # negation: Q' = Q_{-b}(-s_o) carries the same labels, reversed order
            valn = quad(-b, -so)
            check(all(valn[l] == -val[l] for l in LAB), "nu carries Q_b(s_o) onto Q_{-b}(-s_o) label by label")
            check(build(valn, "fwd", "asc") == Ffd, "nu = order reversal on the labels")
            check(conv(Ffa) == build(valn, "bwd", "asc"), "EXACT: converse = negation o time-reversal")
            h = ham(V, Ffa)
            check(h == 2 + (1 if u > 0 else -1), "H(F+) = 2 + sgn(b s_o)")
            check(h == 1 + 2 * c3(V, Ffa), "OCF")
    print("  |s_o| <= %d, both sheets: %d non-degenerate quadruples; degenerate exactly at |s_o| = 1" % (R, nondeg))
    print("  (b=+1: {1,2} the 2-cycle, {-1,0} two fixed points; b=-1: {0,1} fixed points, {-2,-1} the 2-cycle)")
    names = {("so", "se", "te", "to"): None}
    for key in sorted(pats):
        check(len(pats[key]) == 1, "one order pattern per (sheet, side)")
        p = next(iter(pats[key]))
        print("  sheet %+d, %s side: increasing order %s" % (key[0], key[1], " < ".join(p)))
    print("  Theorem A (checked on every quadruple above):")
    print("   (a) map arcs + order (asc or desc): C3+SNK iff b*s_o >= 3 (T_b reverses the pair), TT iff b*s_o <= -3")
    print("   (b) inverse-tree arcs (t_e -D-> s_e, t_o -E-> s_o) + order: SRC+C3 iff b*s_o >= 3, TT otherwise")
    print("   (c) rho(x) = (s_o+s_e) - x swaps s_o<->s_e, t_o<->t_e (AM-fairness), reverses order, fixes the map-arc")
    print("       set: rho: (M,O) ~ (M,O^op) and rho: (F+)^op ~ F-.  So TIME REVERSAL = CONVERSE up to rho.")
    print("   (d) nu(x) = -x carries Q_b(s_o) onto Q_-b(-s_o) with the same labels and the reversed order:")
    print("       F+(nu Q) = F+_desc(Q) ~ F+(Q) (same class): NEGATION FIXES THE CLASS.")
    print("   (e) EXACT labelled identity: (F+(Q))^op = F-(nu Q), i.e. converse = negation o time-reversal.")
    print("   (f) H(F+) = 2 + sgn(b s_o): on positive pairs H = 3 on the 3n+1 sheet, H = 1 on the 3n-1 sheet.")
    print("   (g) (s_o - s_e)(T(s_o) - T(s_e)) = -b s_o: the pair is order-reversed by T_b iff b s_o > 0.")


# ---------------------------------------------------------------- A3 rule census
def section_A3():
    print("=" * 100)
    print("A3  exhaustive census of local rules (2 map directions x 4^4 orientations of the non-map pairs)")
    reps = {"P1 (+,pos)": (1, 5), "P2 (-,pos)": (-1, 7), "P3 (+,neg)": (1, -7), "P4 (-,neg)": (-1, -5)}
    vals = {k: quad(b, so) for k, (b, so) in reps.items()}
    # nu pairs: P1 <-> P4, P2 <-> P3
    check(all(vals["P4 (-,neg)"][l] == -vals["P1 (+,pos)"][l] for l in LAB), "P4 = nu P1")
    check(all(vals["P3 (+,neg)"][l] == -vals["P2 (-,pos)"][l] for l in LAB), "P3 = nu P2")
    # comparison vector: which non-map comparisons depend on sheet / side
    print("  order comparisons of the non-map pairs (value of u < v for pair (u,v)):")
    for name, (u, v) in NONMAP.items():
        row = [vals[k][u] < vals[k][v] for k in reps]
        print("    %s = %s-%s : %s" % (name, u, v, "  ".join("%s:%s" % (k[:2], int(r)) for k, r in zip(reps, row))))
    print("  => pair A (the source pair) is decided by the SHEET alone; B, C, D by the SIDE alone.")
    # verify on a range that every (sheet, side) has one comparison vector
    for b in (1, -1):
        for so in list(range(-301, -1, 2)) + list(range(3, 303, 2)):
            val = quad(b, so)
            if not distinct(val):
                continue
            key = "P1 (+,pos)" if (b == 1 and so > 0) else "P2 (-,pos)" if (b == -1 and so > 0) else \
                  "P3 (+,neg)" if b == 1 else "P4 (-,neg)"
            for name, (u, v) in NONMAP.items():
                check((val[u] < val[v]) == (vals[key][u] < vals[key][v]), "comparison vector constant")
    opts = ("asc", "desc", "L0", "L1")
    rules = []
    for mu in ("fwd", "bwd"):
        for combo in itertools.product(opts, repeat=4):
            rules.append((mu, dict(zip("ABCD", combo))))
    check(len(rules) == 512, "512 rules")
    V = list(LAB)
    table = []
    for mu, rule in rules:
        F = {k: build(vals[k], mu, rule=rule) for k in reps}
        C = {k: cls(V, F[k]) for k in reps}
        table.append((mu, rule, F, C))
    diam = {SRC, SNK}
    def conv_pair(c1, c2):
        return {c1, c2} == diam
    comparisons = [("same side, positive: P1 vs P2", "P1 (+,pos)", "P2 (-,pos)"),
                   ("same side, negative: P3 vs P4", "P3 (+,neg)", "P4 (-,neg)"),
                   ("nu-related: P1 vs P4", "P1 (+,pos)", "P4 (-,neg)"),
                   ("nu-related: P2 vs P3", "P2 (-,pos)", "P3 (+,neg)"),
                   ("one sheet, both sides: P1 vs P3", "P1 (+,pos)", "P3 (+,neg)"),
                   ("one sheet, both sides: P2 vs P4", "P2 (-,pos)", "P4 (-,neg)")]
    print("  rules giving EXACTLY the converse pair {SRC+C3, C3+SNK} on the two quadruples compared:")
    results = {}
    for title, k1, k2 in comparisons:
        hits = [(mu, rule) for mu, rule, F, C in table if conv_pair(C[k1], C[k2])]
        dist = max(len(F[k1] - F[k2]) for mu, rule, F, C in table)
        results[title] = hits
        pure = [h for h in hits if all(r in ("asc", "desc") for r in h[1].values())]
        print("    %-36s %3d rules (%d without fixed-label arcs); max labelled arc distance over all rules = %d"
              % (title, len(hits), len(pure), dist))
        if title.startswith("same side"):
            check(len(hits) == 0 and dist <= 1, "no rule separates the sheets on one side by the converse pair")
    # describe the gauge hits
    short = {TT: "TT", SRC: "SRC", SNK: "SNK", STR: "STR"}
    isL = lambda x: x in ("L0", "L1")
    for title in ("nu-related: P1 vs P4", "one sheet, both sides: P1 vs P3"):
        hits = results[title]
        check(all(any(isL(v) for v in r.values()) for mu, r in hits), "converse pairs off the sheet axis need a gauge arc")
        print("  the %d hits for '%s' (all use fixed-label arcs):" % (len(hits), title))
        for mu, r in hits:
            if mu != "fwd":
                continue
            C = {k: short[cls(V, build(vals[k], mu, rule=r))] for k in reps}
            print("      mu=fwd A=%-4s B=%-4s C=%-4s D=%-4s  classes P1..P4: %s" % (
                r["A"], r["B"], r["C"], r["D"], " ".join(C[k] for k in reps)))
        print("      (and the same 8 with mu=bwd, classes converted)")
    # Type I / Type II characterisation of the nu-related hits
    hits = results["nu-related: P1 vs P4"]
    typeI = [(mu, r) for mu, r in hits if not isL(r["A"]) and not isL(r["B"]) and r["A"] == r["B"]
             and isL(r["C"]) and isL(r["D"]) and r["C"] != r["D"]]
    typeII = [(mu, r) for mu, r in hits if isL(r["A"]) and isL(r["B"]) and r["A"] != r["B"]
              and not isL(r["C"]) and r["C"] == r["D"]]
    check(len(typeI) == 8 and len(typeII) == 8 and len(hits) == 16, "hits = Type I (8) + Type II (8)")
    hits23 = results["nu-related: P2 vs P3"]
    typeIp = [(mu, r) for mu, r in hits23 if not isL(r["A"]) and not isL(r["B"]) and r["A"] != r["B"]
              and isL(r["C"]) and isL(r["D"]) and r["C"] != r["D"]]
    typeII23 = [(mu, r) for mu, r in hits23 if isL(r["A"]) and isL(r["B"]) and r["A"] != r["B"]
                and not isL(r["C"]) and r["C"] == r["D"]]
    check(len(typeIp) == 8 and len(typeII23) == 8 and len(hits23) == 16, "P2/P3 hits = Type I' (8) + Type II (8)")
    key = lambda h: (h[0], tuple(sorted(h[1].items())))
    check({key(h) for h in typeII23} == {key(h) for h in typeII}, "the Type II rules serve both nu-related comparisons")
    allnu = [(mu, r, C) for mu, r, F, C in table if all(C[k2] == CONV[C[k1]] for k1, k2 in
             (("P1 (+,pos)", "P4 (-,neg)"), ("P2 (-,pos)", "P3 (+,neg)")))]
    withdia = [(mu, r) for mu, r, C in allnu if any(c in diam for c in C.values())]
    vacuous = len(allnu) - len(withdia)
    check({key(h) for h in withdia} == {key(h) for h in typeII}, "nu = converse (with diamonds) exactly for Type II")
    print("  the 16 hits for 'nu-related: P2 vs P3' = Type II (the same 8) + Type I' (8): level pairs A, B by OPPOSITE")
    print("  order conventions, cross pairs by the lineage gauge; these put the opposite diamonds on the expanding sides.")
    print("  nu acts as the converse on ALL four patterns, with diamonds present, for exactly the 8 Type II rules")
    print("  (%d further rules satisfy it vacuously, with only the self-converse classes TT and STRONG)." % vacuous)
    print("  In all, 24 distinct gauge rules give the converse pair on some nu-related comparison; none is pure-order.")
    print("  => Type I (8): level pairs A,B by one order convention, cross pairs C={s_o,t_e}, D={s_e,t_o} by a")
    print("     LINEAGE gauge (both from the odd lineage {s_o,t_o} to the even lineage {s_e,t_e}, or both back):")
    print("     the contracting sides P1 (3n+1, x>0) and P4 (3n-1, x<0) get opposite diamonds; the expanding")
    print("     sides get TT / STRONG.  Which sheet gets which diamond flips with the gauge direction.")
    print("     Type II (8): level pairs by a fixed label gauge, cross pairs by one order convention: the class")
    print("     depends on the side only (P1 = P2, P3 = P4); the two sheets are IDENTICAL on each side.")
    # involution dictionary over the 32 pure-order rules
    print("  involution dictionary on the 32 pure-order rules (no fixed-label arc);")
    print("  rho-symmetric = the rule orients C and D by the same convention (rho swaps C <-> D):")
    nuq = {"P1 (+,pos)": "P4 (-,neg)", "P4 (-,neg)": "P1 (+,pos)", "P2 (-,pos)": "P3 (+,neg)", "P3 (+,neg)": "P2 (-,pos)"}
    pure_rules = [(mu, r) for mu, r in rules if all(v in ("asc", "desc") for v in r.values())]
    check(len(pure_rules) == 32, "32 pure-order rules")
    stats = {}
    for mu, r in pure_rules:
        sym = "rho-symmetric" if r["C"] == r["D"] else "mixed C/D"
        tmu = "bwd" if mu == "fwd" else "fwd"
        a1 = a2 = a3 = a4 = True
        for k in reps:
            F = build(vals[k], mu, rule=r)
            Fn = build(vals[nuq[k]], mu, rule=r)
            Ft = build(vals[k], tmu, rule=r)
            Ftn = build(vals[nuq[k]], tmu, rule=r)
            c, cn, ct = cls(V, F), cls(V, Fn), cls(V, Ft)
            a1 &= (cn == c)
            a2 &= (cn == CONV[c])
            a3 &= (ct == CONV[c])
            a4 &= (conv(F) == Ftn)
        for key, val in (("nu fixes the class", a1), ("nu = converse on classes", a2),
                         ("tau = converse on classes", a3), ("EXACT (F_R(Q))^op = F_tauR(nu Q)", a4)):
            stats.setdefault((key, sym), 0)
            stats[(key, sym)] += val
    for key in ("nu fixes the class", "nu = converse on classes", "tau = converse on classes",
                "EXACT (F_R(Q))^op = F_tauR(nu Q)"):
        print("    %-36s rho-symmetric %2d / 16   mixed C/D %2d / 16" % (
            key, stats[(key, "rho-symmetric")], stats[(key, "mixed C/D")]))
    check(stats[("EXACT (F_R(Q))^op = F_tauR(nu Q)", "rho-symmetric")] == 16 and
          stats[("EXACT (F_R(Q))^op = F_tauR(nu Q)", "mixed C/D")] == 16, "converse = nu o tau for every pure-order rule")
    check(stats[("nu fixes the class", "rho-symmetric")] == 16 and stats[("tau = converse on classes", "rho-symmetric")] == 16,
          "rho-symmetric pure-order rules: nu ~ id and tau ~ converse")
    check(stats[("nu = converse on classes", "rho-symmetric")] == 0 and stats[("nu = converse on classes", "mixed C/D")] == 0,
          "no pure-order rule makes negation the converse")
    # positive-side class pairs (sheet +, sheet -) realised by pure-order rules
    pairs = {}
    for mu, r in pure_rules:
        c1 = short[cls(V, build(vals["P1 (+,pos)"], mu, rule=r))]
        c2 = short[cls(V, build(vals["P2 (-,pos)"], mu, rule=r))]
        pairs[(c1, c2)] = pairs.get((c1, c2), 0) + 1
    print("    (3n+1, 3n-1) class pairs on positive pairs over the 32 pure-order rules: %s"
          % ", ".join("%s/%s x%d" % (a, b, n) for (a, b), n in sorted(pairs.items())))
    check(all({a, b} != {"SRC", "SNK"} for (a, b) in pairs), "never the converse pair")
    mixed_dia = [(mu, r) for mu, r in pure_rules if r["C"] != r["D"]
                 and any(cls(V, build(vals[k], mu, rule=r)) in diam for k in reps)]
    check(not mixed_dia, "mixed C/D rules give no diamond")
    print("    mixed C/D rules produce only TT and STRONG (no diamond at all).")
    return results


# ---------------------------------------------------------------- A4
def generic_class(Vvals, arcs_map, mu="fwd", order="asc"):
    """Vvals: dict label->value; arcs_map: list of (u,v) dynamical arcs; others by order."""
    arcs = set()
    mp = set()
    for (u, v) in arcs_map:
        arcs.add((u, v) if mu == "fwd" else (v, u))
        mp.add(frozenset((u, v)))
    L = list(Vvals)
    for u, v in itertools.combinations(L, 2):
        if frozenset((u, v)) in mp:
            continue
        up = Vvals[u] < Vvals[v]
        if order == "desc":
            up = not up
        arcs.add((u, v) if up else (v, u))
    return cls(L, frozenset(arcs)), frozenset(arcs)


def section_A4():
    print("=" * 100)
    print("A4  other natural four-vertex constructions")
    # --- A4.1 the four up/down choices on a pair {2i-1, 2i}, each member moving by i = ceil(n/2)
    print("A4.1 the four up/down choices on a pair {2i-1, 2i} (each member moves by its length i):")
    print("     (odd up, even down) = 3n+1;  (odd down, even up) = THM-4470 all-flipped = 3n-1 shifted by one;")
    print("     (up, up) = Mahler ceil(3n/2) (THM-2228);  (down, down) = floor(n/2)")
    choices = {"3n+1 (odd^, even v)": (1, -1), "flipped (odd v, even ^)": (-1, 1),
               "Mahler (^, ^)": (1, 1), "floor (v, v)": (-1, -1)}
    rows = {}
    for name, (do, de) in choices.items():
        res = {}
        for side, irange in (("pos", range(2, 400)), ("neg", range(-400, -1))):
            seen = set()
            for i in irange:
                so, se = 2 * i - 1, 2 * i
                val = {"so": so, "se": se, "to": so + do * i, "te": se + de * i}
                if len(set(val.values())) < 4:
                    continue
                for mu in ("fwd", "bwd"):
                    c, _ = generic_class(val, [("so", "to"), ("se", "te")], mu, "asc")
                    c2, _ = generic_class(val, [("so", "to"), ("se", "te")], mu, "desc")
                    seen.add((mu, c, c2))
            fw = {x[1] for x in seen if x[0] == "fwd"}
            bw = {x[1] for x in seen if x[0] == "bwd"}
            fwd2 = {x[2] for x in seen if x[0] == "fwd"}
            check(len(fw) == 1 and len(bw) == 1, "class constant per choice/side")
            res[side] = (fw.pop(), bw.pop(), fwd2.pop())
        rows[name] = res
        fair = (do + de == 0)
        print("   %-26s AM-fair=%-5s  positive: fwd %-16s bwd %-16s | negative: fwd %-16s bwd %-16s"
              % (name, fair, res["pos"][0], res["pos"][1], res["neg"][0], res["neg"][1]))
    check(rows["3n+1 (odd^, even v)"]["pos"][:2] == (SNK, SRC), "Collatz pair: forward C3+SNK, backward SRC+C3")
    check(rows["flipped (odd v, even ^)"]["pos"][:2] == (TT, TT), "flipped: transitive")
    check(rows["Mahler (^, ^)"]["pos"][:2] == (TT, STR), "Mahler: TT / STRONG")
    check(rows["floor (v, v)"]["pos"][:2] == (STR, TT), "floor: STRONG / TT")
    print("   => the diamonds occur exactly for the crossing pattern (lower member up, upper member down);")
    print("      forward/backward = C3+SNK / SRC+C3, i.e. time reversal = converse.  For Mahler's map (parallel")
    print("      pattern) order reversal maps the arc set M onto M^op, so negation ~ time reversal and both")
    print("      swap TT <-> STRONG, while the converse fixes both self-converse classes.")

    # --- A4.2 matchings of 2 directed arcs on 4 ordered points: all 12 configurations
    print("A4.2 all 12 ways to place two disjoint directed arcs on 4 ordered points (other pairs ascending):")
    pos = {1: 1, 2: 2, 3: 3, 4: 4}
    conf = {}
    for m in ([(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]):
        for o1 in (0, 1):
            for o2 in (0, 1):
                arcs = [m[0] if o1 == 0 else m[0][::-1], m[1] if o2 == 0 else m[1][::-1]]
                c, _ = generic_class(pos, arcs, "fwd", "asc")
                conf[tuple(arcs)] = c
                print("     arcs %-18s -> %s" % (" , ".join("%d->%d" % a for a in arcs), c))
    check(conf[((3, 1), (2, 4))] == SNK, "Collatz pattern (t_e<s_o<s_e<t_o, arcs 3->1 and 2->4)")
    print("     (Collatz forward, contracting side: t_e<s_o<s_e<t_o with 2->4, 3->1 -> C3+SNK;")
    print("      the only diamond-producing configurations are the crossing matching {13, 24} with antiparallel arcs)")
    dia = [a for a, c in conf.items() if c in (SRC, SNK)]
    check(all(set(frozenset(x) for x in a) == {frozenset((1, 3)), frozenset((2, 4))} for a in dia), "diamonds need {13,24}")
    check(all((a[0][1] - a[0][0]) * (a[1][1] - a[1][0]) < 0 for a in dia), "diamonds need antiparallel arcs")

    # --- A4.3 step decomposition n -> 3n -> 3n+b -> (3n+b)/2 on odd n
    print("A4.3 step quadruple {n, 3n, 3n+b, (3n+b)/2} (odd n; arcs x3, +-1, /2; other pairs by order):")
    for b in (1, -1):
        for side, rng in (("pos", range(3, 801, 2)), ("neg", range(-801, -2, 2))):
            seen = set()
            for n in rng:
                val = {"n": n, "3n": 3 * n, "3n+b": 3 * n + b, "h": (3 * n + b) // 2}
                if len(set(val.values())) < 4:
                    continue
                arcsM = [("n", "3n"), ("3n", "3n+b"), ("3n+b", "h")]
                row = tuple(generic_class(val, arcsM, mu, od)[0] for mu in ("fwd", "bwd") for od in ("asc", "desc"))
                seen.add(row)
            check(len(seen) == 1, "step quadruple class constant")
            r = seen.pop()
            print("   b=%+d %s: fwd/asc %-16s fwd/desc %-16s bwd/asc %-16s bwd/desc %-16s" % (b, side, *r))
    print("   => both sheets give the same class on a given side; no converse pair across sheets.")

    # --- A4.4 inverse-tree fork {T(x), x, 2x, E(x)}
    print("A4.4 inverse-tree fork {T_b(x), x, 2x, E_b(x)} (x in the legal class; arcs 2x->x, E(x)->x, x->T(x)):")
    for b in (1, -1):
        for side, rng in (("pos", range(5, 3000)), ("neg", range(-3000, -4))):
            seen = {}
            for x in rng:
                if (2 * x - b) % 3 != 0:
                    continue
                E = (2 * x - b) // 3
                t = int(Tb(x, b))
                val = {"T": t, "x": x, "2x": 2 * x, "E": E}
                if len(set(val.values())) < 4:
                    continue
                arcsM = [("2x", "x"), ("E", "x"), ("x", "T")]
                c = generic_class(val, arcsM, "fwd", "asc")[0]
                want = SRC if x % 2 else (STR if side == "pos" else TT)
                check(c == want, "fork class: SRC+C3 for odd x; STRONG (x>0) / TT (x<0) for even x")
                key = ("odd x" if x % 2 else "even x") + " -> " + c.split("(")[0]
                seen[key] = seen.get(key, 0) + 1
            print("   b=%+d %s: forward/asc %s" % (b, side, dict(sorted(seen.items()))))
    print("   (3 dynamical arcs; the class is set by the parity of x and the side, identically on both sheets)")
    # nu-check for forks: F(b, x) and F(-b, -x) are order-reversed images
    for x in range(5, 400):
        for b in (1, -1):
            if (2 * x - b) % 3 != 0:
                continue
            val = {"T": int(Tb(x, b)), "x": x, "2x": 2 * x, "E": (2 * x - b) // 3}
            valn = {"T": int(Tb(-x, -b)), "x": -x, "2x": -2 * x, "E": (-2 * x + b) // 3}
            check(all(valn[k] == -val[k] for k in val), "nu carries forks label by label")
            if len(set(val.values())) < 4:
                continue
            arcsM = [("2x", "x"), ("E", "x"), ("x", "T")]
            check(generic_class(valn, arcsM, "fwd", "asc")[1] == generic_class(val, arcsM, "fwd", "desc")[1],
                  "nu = order reversal on forks")
            check(conv(generic_class(val, arcsM, "fwd", "asc")[1]) == generic_class(valn, arcsM, "bwd", "asc")[1],
                  "converse = nu o tau on forks")

    # --- A4.5 length-4 orbit windows {x, Tx, T^2x, T^3x}
    print("A4.5 length-4 orbit windows {x, Tx, T^2x, T^3x} (path arcs + ascending order), x in [10^4, 10^4 + 2^12):")
    for b in (1, -1):
        byword = {}
        for x in range(10 ** 4, 10 ** 4 + 4096):
            xs = [Fr(x)]
            for _ in range(3):
                xs.append(Tb(xs[-1], b))
            w = "".join(str(int(v.numerator % 2)) for v in xs[:3])
            val = {j: xs[j] for j in range(4)}
            c = generic_class(val, [(0, 1), (1, 2), (2, 3)], "fwd", "asc")[0]
            byword.setdefault(w, set()).add(c)
        check(all(len(s) == 1 for s in byword.values()), "window class is a word function")
        print("   b=%+d: " % b + "  ".join("%s:%s" % (w, next(iter(s)).split("(")[0]) for w, s in sorted(byword.items())))
        if b == 1:
            plusw = {w: next(iter(s)) for w, s in byword.items()}
        else:
            minusw = {w: next(iter(s)) for w, s in byword.items()}
    check(plusw == minusw, "length-4 window classes identical on the two sheets word by word")
    print("   => word by word the two sheets give identical classes (generic windows); the sheet is invisible.")

    # --- A4.6 b = 0 on Q (T_0(x) = 3x/2 on odd x) with both pairings, and Mahler ceil/floor
    print("A4.6 b = 0 (T_0(x) = x/2 or 3x/2 on Q): no AM-fair pairing, so the partner offset beta is a free gauge:")
    for beta in (1, -1):
        out = {}
        for side, rng in (("pos", range(5, 401, 2)), ("neg", range(-401, -4, 2))):
            cs = set()
            for so in rng:
                val = quad(0, so, beta)
                if not distinct(val):
                    continue
                cs.add((cls(list(LAB), build(val, "fwd", "asc")), cls(list(LAB), build(val, "bwd", "asc"))))
            check(len(cs) == 1, "b=0 class constant")
            out[side] = cs.pop()
        print("   beta=%+d: positive fwd/bwd %s / %s ; negative fwd/bwd %s / %s" % (beta, *out["pos"], *out["neg"]))
    print("   => for b = 0 the class is beta*sign(s_o): the diamonds appear, but which side carries them is the")
    print("      pairing gauge.  For b = +-1 AM-fairness forces beta = b, and the class becomes b*sign(s_o).")
    print("   Mahler's map g(n) = ceil(3n/2) (THM-2228) and its nu-conjugate floor(3n/2), with either pairing:")
    for gname, g in (("ceil(3n/2)", lambda n: -((-3 * n) // 2)), ("floor(3n/2)", lambda n: (3 * n) // 2)):
        for off in (0, 1):
            out = {}
            for side, rng in (("pos", range(3, 400)), ("neg", range(-400, -2))):
                cs = set()
                for i in rng:
                    a, c = 2 * i - 1 + off, 2 * i + off
                    so, se = (a, c) if a % 2 else (c, a)
                    val = {"so": so, "se": se, "to": g(so), "te": g(se)}
                    if len(set(val.values())) < 4:
                        continue
                    cs.add((cls(list(LAB), build(val, "fwd", "asc")), cls(list(LAB), build(val, "bwd", "asc"))))
                check(len(cs) == 1, "Mahler class constant")
                out[side] = cs.pop()
            print("   %-11s pairs {2i-1+%d, 2i+%d}: positive fwd/bwd %s / %s ; negative fwd/bwd %s / %s"
                  % (gname, off, off, *out["pos"], *out["neg"]))
            if gname == "ceil(3n/2)":
                check(out["pos"] == (TT, STR), "Mahler positive: TT / STRONG")

    # --- A4.7 does the diamond see the multiplier?  maps (q n + r)/2 on odds, n/2 on evens
    print("A4.7 multiplier test: f(n) = n/2 (even), (q n + r)/2 (odd), pairs {2i-1, 2i} (odd member lower) or")
    print("     {2i, 2i+1} (odd member upper), 10 <= i <= 300 (none of these is AM-fair unless (q,r) = (3,+-1)):")
    for q in (1, 3, 5, 7, 9):
        for r in (-1, 1):
            cells = []
            for off, oname in ((0, "odd lower"), (1, "odd upper")):
                cs = set()
                for i in range(10, 301):
                    a, c = 2 * i - 1 + off, 2 * i + off
                    so, se = (a, c) if a % 2 else (c, a)
                    val = {"so": so, "se": se, "to": (q * so + r) // 2, "te": se // 2}
                    if len(set(val.values())) < 4:
                        continue
                    cs.add((cls(list(LAB), build(val, "fwd", "asc")), cls(list(LAB), build(val, "bwd", "asc"))))
                if not cs:
                    check(q == 1, "only q = 1 degenerates")
                    cells.append("%s: degenerate (t_o = t_e)" % oname)
                    continue
                check(len(cs) == 1, "class constant for large i")
                f, bk = cs.pop()
                cells.append("%s: %s/%s" % (oname, f.split("(")[0], bk.split("(")[0]))
                if q >= 3:
                    check((f, bk) == ((SNK, SRC) if off == 0 else (TT, TT)), "q >= 3: diamonds iff odd member lower")
            print("     q=%d r=%+d  %s" % (q, r, " | ".join(cells)))
    print("   => every odd q >= 3 gives the same diamond pair on the pairing whose odd member is the lower one;")
    print("      H = 3 does not see the multiplier 3 (ANALOGY).  q = 3 is special only through AM-fairness, which")
    print("      makes the order-reversing isomorphism the affine reflection rho and forces the pairing offset = b.")


def main():
    section_A1()
    section_A2()
    section_A3()
    section_A4()
    print("=" * 100)
    print("PART A: ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
