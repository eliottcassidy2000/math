"""a10_hyp_c_probe.py (AUDIT) - cheap hostile probe of HYP-9664 (C) at n = 5: can more than 50 of the 155 planes of
PG(4,2) be odd at once (a) for arbitrary tournaments, (b) for tournaments with all 155 lines cyclic?

Randomised greedy: planes in random order, each added if the F_2 system stays consistent (own elimination);
then 1-exchange improvement.  A set of more than 50 simultaneously odd planes with all lines cyclic would REFUTE
the hypothesis (fewer than 630 even 4-sets); without the line condition it would refute its second clause as
worded in the body ('whether more than 50 of the 155 planes can be odd at once').
Any set found is verified on an explicit tournament by direct counts of cyclic triangles."""
import random, sys, os, itertools, time

HERE = os.path.dirname(os.path.abspath(__file__))
N = 32
PTS = list(range(1, N))
PAIRS = [(u, v) for u in PTS for v in PTS if u < v]
PIDX = {p: i for i, p in enumerate(PAIRS)}
NV = len(PAIRS)
out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def arc_term(u, v):
    if u < v:
        return (1 << PIDX[(u, v)], 0)
    return (1 << PIDX[(v, u)], 1)


def pairing_form(s1, s2):
    m = c = 0
    for u in s1:
        for v in s2:
            mm, cc = arc_term(u, v)
            m ^= mm
            c ^= cc
    return m, c


planes = sorted({tuple(sorted({a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c})) for a in PTS for b in PTS for c in PTS
                 if a < b < c and c != a ^ b})
assert len(planes) == 155
lines = sorted({tuple(sorted((a, b, a ^ b))) for a in PTS for b in PTS if a < b})
assert len(lines) == 155
plane_rows = []
for pl in planes:
    rows = []
    for a, b, c in itertools.combinations(pl, 3):
        d = a ^ b ^ c
        if d > c and d in pl and d not in (a, b, c):
            m1, c1 = pairing_form((a, b), (c, d))
            m2, c2 = pairing_form((a, c), (b, d))
            rows += [(m1, c1 ^ 1), (m2, c2 ^ 1)]
    assert len(rows) == 14
    plane_rows.append(rows)
line_rows = []
for (a, b, c) in lines:
    m1, c1 = arc_term(a, b)
    m2, c2 = arc_term(b, c)
    m3, c3 = arc_term(c, a)
    line_rows += [(m1 ^ m2, c1 ^ c2), (m2 ^ m3, c2 ^ c3)]


class Elim:
    def __init__(self, piv=None):
        self.piv = dict(piv) if piv else {}

    def copy(self):
        return Elim(self.piv)

    def add(self, m, c):
        piv = self.piv
        while m:
            hb = m.bit_length() - 1
            pr = piv.get(hb)
            if pr is None:
                piv[hb] = (m, c)
                return True
            m ^= pr[0]
            c ^= pr[1]
        return c == 0

    def add_rows(self, rows):
        ok = True
        for m, c in rows:
            if not self.add(m, c):
                ok = False
        return ok

    def solution(self):
        t = 0
        for hb in sorted(self.piv):
            m, c = self.piv[hb]
            val = c ^ (bin(t & (m ^ (1 << hb))).count("1") & 1)
            if val:
                t |= 1 << hb
        return t


def c3(arc, S):
    n = 0
    for u, v, w in itertools.combinations(S, 3):
        if (arc(u, v) and arc(v, w) and arc(w, u)) or (arc(v, u) and arc(w, v) and arc(u, w)):
            n += 1
    return n


def verify(t):
    def arc(u, v):
        if u < v:
            return bool(t >> PIDX[(u, v)] & 1)
        return not (t >> PIDX[(v, u)] & 1)
    odd = 0
    ev = 0
    for pl in planes:
        bad = 0
        for a, b, c in itertools.combinations(pl, 3):
            d = a ^ b ^ c
            if d > c and d in pl and d not in (a, b, c):
                if c3(arc, (a, b, c, d)) % 2 == 0:
                    bad += 1
        ev += bad
        odd += (bad == 0)
    cyc = sum(1 for l in lines if c3(arc, l) == 1)
    return odd, ev, cyc


rng = random.Random(11)
for label, base_rows in (("all lines cyclic", line_rows), ("arbitrary tournaments", [])):
    base = Elim()
    assert base.add_rows(base_rows)
    best = 0
    best_set = None
    hist = {}
    t0 = time.time()
    trials = 0
    while time.time() - t0 < (150 if label.startswith("all") else 150):
        trials += 1
        order = list(range(155))
        rng.shuffle(order)
        e = base.copy()
        chosen = []
        for p in order:
            e2 = e.copy()
            if e2.add_rows(plane_rows[p]):
                e = e2
                chosen.append(p)
        hist[len(chosen)] = hist.get(len(chosen), 0) + 1
        if len(chosen) > best:
            best = len(chosen)
            best_set = (chosen, e)
    P("%s: %d random greedy maximal sets of simultaneously odd planes; sizes -> count %s; largest %d" % (label, trials, sorted(hist.items()), best))
    chosen, e = best_set
    t = e.solution()
    odd, ev, cyc = verify(t)
    P("   explicit tournament for the largest set (free variables 0), direct recount: odd planes %d, even 4-sets %d, cyclic lines %d/155" % (odd, ev, cyc))
open(os.path.join(HERE, "a10_hyp_c_probe.out"), "w").write("\n".join(out) + "\n")
