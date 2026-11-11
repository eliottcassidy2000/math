"""a10_hyp_c_refute.py (AUDIT) - HYP-9664 (C) is FALSE at n = 5 (and n = 6): an explicit structural counterexample.

HYP-9664 (C): 'for every n >= 4 the Cayley-Dickson sign tournament minimises the number of even zero-sum 4-sets
among the tournaments on PG(n-1, 2) with all lines cyclic, the minimum being 42 [n-1 choose 3]_2, and at most
[n choose 3]_2 - 7 [n-1 choose 3]_2 planes are odd at once'  (630 even 4-sets and 50 odd planes at n = 5).

Construction tested here (the generalisation of type A of Theorem 3.4): fix a subspace U of dimension n - 3 of
Z_2^n (a point for n = 4, a line for n = 5, a plane for n = 6) and ask that ALL the 2^(3(n-3)) planes skew to U
(8, 64, 512) be odd, together with all lines cyclic.  The system is solved by own elimination over F_2 and the
solution is verified on the explicit tournament by DIRECT counts of cyclic triangles (no linear algebra).
The witness for n = 5 is written to a10_witness_n5.txt (465-bit word; pair (u < v) number i in lexicographic order,
bit i = 1 iff u -> v)."""
import sys, os, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def gauss(m, k):
    num = den = 1
    for i in range(k):
        num *= 2 ** (m - i) - 1
        den *= 2 ** (i + 1) - 1
    return num // den


def run(n):
    N = 1 << n
    PTS = list(range(1, N))
    PAIRS = [(u, v) for u in PTS for v in PTS if u < v]
    PIDX = {p: i for i, p in enumerate(PAIRS)}

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
    planes = set()
    for a in PTS:
        for b in PTS:
            if b <= a:
                continue
            for c in PTS:
                if c <= b or c == a ^ b:
                    continue
                planes.add(tuple(sorted({a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c})))
    planes = sorted(planes)
    lines = sorted({tuple(sorted((a, b, a ^ b))) for a in PTS for b in PTS if a < b})
    assert len(planes) == gauss(n, 3) and len(lines) == gauss(n, 2)
    U = set(range(1, 1 << (n - 3)))               # nonzero vectors of the span of the first n - 3 basis vectors
    skew = [pl for pl in planes if not (set(pl) & U)]
    piv = {}

    def add(m, c):
        while m:
            hb = m.bit_length() - 1
            pr = piv.get(hb)
            if pr is None:
                piv[hb] = (m, c)
                return True
            m ^= pr[0]
            c ^= pr[1]
        return c == 0
    ok = True
    for (a, b, c) in lines:
        m1, c1 = arc_term(a, b)
        m2, c2 = arc_term(b, c)
        m3, c3 = arc_term(c, a)
        ok &= add(m1 ^ m2, c1 ^ c2)
        ok &= add(m2 ^ m3, c2 ^ c3)
    rank_lines = len(piv)
    for pl in skew:
        for a, b, c in itertools.combinations(pl, 3):
            d = a ^ b ^ c
            if d > c:
                m1, c1 = pairing_form((a, b), (c, d))
                m2, c2 = pairing_form((a, c), (b, d))
                ok &= add(m1, c1 ^ 1)
                ok &= add(m2, c2 ^ 1)
    P("n = %d: planes %d, planes skew to a fixed subspace of dimension %d: %d = 2^%d; system [all lines cyclic + these planes odd] consistent: %s (rank %d of %d unknowns; lines alone %d)"
      % (n, len(planes), n - 3, len(skew), 3 * (n - 3), ok, len(piv), len(PAIRS), rank_lines))
    if not ok:
        return
    t = 0
    for hb in sorted(piv):
        m, c = piv[hb]
        if c ^ (bin(t & (m ^ (1 << hb))).count("1") & 1):
            t |= 1 << hb

    def arc(u, v):
        if u < v:
            return bool(t >> PIDX[(u, v)] & 1)
        return not (t >> PIDX[(v, u)] & 1)

    def c3(S):
        k = 0
        for u, v, w in itertools.combinations(S, 3):
            if (arc(u, v) and arc(v, w) and arc(w, u)) or (arc(v, u) and arc(w, v) and arc(u, w)):
                k += 1
        return k
    odd = 0
    ev = 0
    dist = {}
    for pl in planes:
        bad = 0
        for a, b, c in itertools.combinations(pl, 3):
            d = a ^ b ^ c
            if d > c and d in pl:
                if c3((a, b, c, d)) % 2 == 0:
                    bad += 1
        dist[bad] = dist.get(bad, 0) + 1
        ev += bad
        odd += bad == 0
    cyc = sum(1 for l in lines if c3(l) == 1)
    scores = sorted({sum(1 for v in PTS if v != u and arc(u, v)) for u in PTS})
    cd_even = 42 * gauss(n - 1, 3)
    cd_odd = gauss(n, 3) - 7 * gauss(n - 1, 3)
    P("   explicit tournament, DIRECT recount: cyclic lines %d/%d; planes by number of even 4-sets %s; odd planes %d (Cayley-Dickson: %d); even 4-sets %d (Cayley-Dickson: %d); scores %s"
      % (cyc, len(lines), sorted(dist.items()), odd, cd_odd, ev, cd_even, scores))
    P("   => HYP-9664 (C) at n = %d: %s" % (n, "REFUTED (fewer even 4-sets than the Cayley-Dickson tournament, with all lines cyclic)" if (cyc == len(lines) and ev < cd_even) else "not refuted by this construction"))
    if n == 5:
        open(os.path.join(HERE, "a10_witness_n5.txt"), "w").write("%d\n" % t)


for n in (4, 5, 6):
    run(n)
open(os.path.join(HERE, "a10_hyp_c_refute.out"), "w").write("\n".join(out) + "\n")
