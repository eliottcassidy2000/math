"""Independent sharp norm-four physical/network and coarea audit.

Universe: 82 sorted distinct positive primitive ternary-unit triples c<=32
in one of the three norm-four sectors, plus the zero control (1,1,1).
Exact closed profiles are compared with all-roof piecewise-linear envelopes
and the independently validated native-head TSV. No producer imports.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def cross(x, y):
    return (x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0])


def egcd(a, b):
    if not b:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x-(a//b)*y


def bezout(w):
    g, x, y = egcd(w[0], w[1])
    d, u, v = egcd(g, w[2])
    need(d == 1, "primitive speed")
    return (u*x, u*y, v)


def envelope(w, v, r):
    lines = [(2*r/max(w), F(0))]
    for i in range(3):
        j, k = [j for j in range(3) if j != i]
        lines.append((r*F(w[j]+w[k], w[j]*w[k]), F(abs(v[i]), w[j]*w[k])))
    def f(x):
        return max(F(0), min(A-B*x for A, B in lines))
    end = max(A/B for A, B in lines if B)
    cuts = {F(0), end}
    cuts.update(A/B for A, B in lines if B)
    for i, (A, B) in enumerate(lines):
        for C, D in lines[i+1:]:
            if B != D and 0 <= (A-C)/(B-D) <= end:
                cuts.add((A-C)/(B-D))
    return f, sorted(cuts)


def clip(poly, A, B, C):
    """Keep A*x+B*y<=C with exact rational edge intersections."""
    out = []
    for P, Q in zip(poly, poly[1:]+poly[:1]):
        p, q = A*P[0]+B*P[1]-C, A*Q[0]+B*Q[1]-C
        if p <= 0:
            out.append(P)
        if (p < 0 < q) or (q < 0 < p):
            t = p/(p-q)
            out.append((P[0]+t*(Q[0]-P[0]), P[1]+t*(Q[1]-P[1])))
    return out


def coarea_constant(v, r):
    coeff = sorted(map(abs, v))
    A, B, C = coeff
    poly = [(-r, -r), (r, -r), (r, r), (-r, r)]
    poly = clip(poly, A, B, C*r)
    poly = clip(poly, -A, -B, C*r)
    area = abs(sum((P[0]*Q[1]-Q[0]*P[1]
                    for P, Q in zip(poly, poly[1:]+poly[:1])), F(0)))/2
    return area/C


def profile(w, v):
    a, b, c = w
    d = next(i for i in range(3) if abs(v[i]) == 2)
    complements = [w[i] for i in range(3) if i != d]
    need(max(complements) == c, "doubled-coordinate complement includes c")
    p = min(complements)
    q, k0, k1 = F(3, 7*c), F(3*(c-p), 28), F(3*(c+p), 28)
    def closed(x):
        return max(F(0), min(q, F(3*(p+c), 14*p*c)-F(2, p*c)*x))
    full, cuts = envelope(w, v, F(3, 14))
    cuts = sorted(set(cuts) | {k0, k1})
    for x in cuts:
        need(full(x) == closed(x), "pointwise fixed-roof identity at every affine breakpoint")
    # Both functions are affine on each listed interval and identically zero
    # after the final endpoint: these checks certify the entire real half-line.
    I = 2*sum(((R-L)*(full(L)+full(R))/2
               for L, R in zip(cuts, cuts[1:])), F(0))
    need(I == F(9, 98), "exact full profile integral")
    live = [k for k in range(1, k1.numerator//k1.denominator+1)
            if k % 3 and closed(k) > 0]
    mass = 2*sum((closed(k) for k in live), F(0))
    need(mass < F(3, 49)+F(4, 7*c), "strict layer-cake bound")
    return mass, d, 2*len(live)


def main():
    print("source_sha256", sha256(Path(__file__).read_bytes()).hexdigest())
    tsv = ROOT/'05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv'
    print("comparison_tsv_sha256", sha256(tsv.read_bytes()).hexdigest())
    with tsv.open(newline="", encoding="utf-8") as handle:
        native = {(int(row['a']), int(row['b']), int(row['c'])): row for row in csv.DictReader(handle, delimiter='\t')}
    rows = []
    for a in range(1, 33):
        for b in range(a+1, 33):
            for c in range(b+1, 33):
                if gcd(gcd(a, b), c) != 1 or any(x % 3 == 0 for x in (a, b, c)):
                    continue
                sectors = []
                if c == 2*a+b:
                    sectors.append(("c=2a+b", (2, 1, -1)))
                if c == a+2*b:
                    sectors.append(("c=a+2b", (1, 2, -1)))
                if 2*b == a+c:
                    sectors.append(("2b=a+c", (1, -2, 1)))
                need(len(sectors) <= 1, "distinct sorted sectors disjoint")
                for name, v in sectors:
                    w = (a, b, c)
                    z = bezout(w)
                    n1 = cross(v, z)
                    need(sum(x*y for x, y in zip(w, z)) == 1 and cross(w, n1) == v,
                         "integral saturated parameterization")
                    mass, d, count = profile(w, v)
                    need(mass <= F(11, 140), "sharp finite base")
                    row = native[w]
                    denominator = int(row['denominator'])
                    projections = [F(int(row['E'+str(i)+'_numerator']), denominator) for i in range(3)]
                    need(mass == F(int(row['mass_numerator']), denominator), "native physical comparison")
                    need(mass == projections[d] == min(projections), "one fixed projection equals physical")
                    need(count == int(row['raw_carriers']), "strict live-support count")
                    rows.append((w, name, mass, d, count))
    need(len(rows) == 82, "strict-speed finite universe exactly82")
    equality = [w for w, _, mass, _, _ in rows if mass == F(11, 140)]
    need(equality == [(2, 11, 20)], "unique sharp primitive equality")
    for row in rows:
        print("head", *row)
    need(profile((1, 1, 1), (1, -2, 1))[0] == 0, "degenerate repeated-speed zero control")
    for R in [F(i, d) for d in (1, 2, 7, 28) for i in range(1, 60)]:
        N = (R.numerator-1)//R.denominator
        count = N-N//3
        need(count < F(2, 3)*R+F(2, 3), "strict positive open-interval residue bound")
    need(F(11, 140)-(F(3, 49)+F(4, 7*33)) == F(1, 32340), "tail margin at33")
    constants = [((1, 1, -1), (2, 5, 7), F(3)),
                 ((1, -2, 1), (2, 11, 20), F(2)),
                 ((1, 2, -2), (4, 5, 7), F(7, 4))]
    for v, w, K in constants:
        z = bezout(w)
        need(cross(w, cross(v, z)) == v, "coarea lattice normalization")
        for r in (F(3, 14), F(1, 2), F(2, 3)):
            area = coarea_constant(v, r)
            full, cuts = envelope(w, v, r)
            I = 2*sum(((R-L)*(full(L)+full(R))/2
                       for L, R in zip(cuts, cuts[1:])), F(0))
            need(area == I == K*r*r, "independent polygon and full-profile coarea")
            print("coarea", v, "r", r, "integral", I)
    print("strict_distinct_head_cases", len(rows), "equality", equality)
    print("degenerate111_mass", 0, "tail_margin_at33", F(1, 32340))
    print("active_gates", GATES)
    print("PASS exact norm4 physical=minE ceiling11/140, unique(2,11,20), coarea constants")


if __name__ == "__main__":
    main()
