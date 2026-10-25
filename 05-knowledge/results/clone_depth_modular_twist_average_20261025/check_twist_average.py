"""Quadratic twists of square-root covers: the twist average is the branch complement.

Theorem C (odd characteristic): for n linear forms spanning the dual of k^(d+1), the 2^(n-1) quadratic twist classes
of the root variety X (x_j^2 = L_j) satisfy  sum_t #X_t(k) = 2^(n-1) #P^d(k)  and  sum_t #X_t^{etale}(k) = 2^(n-1) #U(k),
U = complement of the branch hyperplanes.  Checked here on the hourglass (7 forms), a six-root K3 subset and all nine
magic forms, over F_q for q = 5..43 prime and q = 25, 49, with q = 3, 9 as characteristic-three controls.
Standard library only; gates are explicit exceptions.  Run normally and with -O; outputs must agree byte for byte.
"""
import itertools, json, sys
from itertools import combinations


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


class GF:
    """F_p or F_{p^2} (elements (a,b) = a + b*alpha, alpha^2 = r a non-residue)."""
    def __init__(self, p, k):
        need(k in (1, 2) and p % 2 == 1, 'odd p, degree 1 or 2')
        self.p, self.k, self.q = p, k, p ** k
        if k == 2:
            self.r = next(a for a in range(2, p) if pow(a, (p - 1) // 2, p) == p - 1)
            self.elements = [(a, b) for a in range(p) for b in range(p)]
            self.zero, self.one = (0, 0), (1, 0)
        else:
            self.elements = list(range(p)); self.zero, self.one = 0, 1

    def add(self, x, y):
        return (x + y) % self.p if self.k == 1 else ((x[0] + y[0]) % self.p, (x[1] + y[1]) % self.p)

    def neg(self, x):
        return (-x) % self.p if self.k == 1 else ((-x[0]) % self.p, (-x[1]) % self.p)

    def mul(self, x, y):
        if self.k == 1:
            return (x * y) % self.p
        return ((x[0] * y[0] + x[1] * y[1] * self.r) % self.p, (x[0] * y[1] + x[1] * y[0]) % self.p)

    def chi(self, x):
        """Quadratic character: 0, +1, -1."""
        if x == self.zero:
            return 0
        e, acc, base = (self.q - 1) // 2, self.one, x
        while e:
            if e & 1:
                acc = self.mul(acc, base)
            base = self.mul(base, base); e >>= 1
        need(acc in (self.one, self.neg(self.one)), 'Euler criterion')
        return 1 if acc == self.one else -1

    def scal(self, c, x):  # c in {-1,0,1}
        return self.zero if c == 0 else (x if c == 1 else self.neg(x))

    def form(self, coeffs, pt):
        v = self.zero
        for c, x in zip(coeffs, pt):
            v = self.add(v, self.scal(c, x))
        return v

    def projective_plane(self):
        pts = [(self.one, x, z) for x in self.elements for z in self.elements]
        pts += [(self.zero, self.one, z) for z in self.elements]
        pts += [(self.zero, self.zero, self.one)]
        need(len(pts) == self.q ** 2 + self.q + 1, '#P^2')
        return pts


MAGIC = {'E': (1, 0, 0), 'E+x': (1, 1, 0), 'E-x': (1, -1, 0), 'E+z': (1, 0, 1), 'E-z': (1, 0, -1),
         'E+x+z': (1, 1, 1), 'E-x-z': (1, -1, -1), 'E+x-z': (1, 1, -1), 'E-x+z': (1, -1, 1)}
HOURGLASS = ['E', 'E+x', 'E-x', 'E+z', 'E-z', 'E+x+z', 'E-x-z']
K3SIX = ['E', 'E+x', 'E-x', 'E+z', 'E-z', 'E+x+z']
NINE = list(MAGIC)


def twist_counts(F, names):
    """For every twist class (eps with eps_0 = 0): (full count, etale count). Also #U and the character vectors."""
    forms = [MAGIC[nm] for nm in names]; n = len(forms)
    vectors = [tuple(F.chi(F.form(f, pt)) for f in forms) for pt in F.projective_plane()]
    need(all(any(v) for v in vectors), 'forms span the dual space: no common zero')
    U = [v for v in vectors if all(v)]
    twists = [eps for eps in itertools.product((0, 1), repeat=n) if eps[0] == 0]
    full, etale = [], []
    for eps in twists:
        f = e = 0
        for v in vectors:
            nz = [(-c if s else c) for c, s in zip(v, eps) if c]
            m = n - len(nz)
            if len(set(nz)) == 1:              # all nonzero characters equal after twisting
                fib = 2 ** (n - m - 1)
                f += fib
                if m == 0:
                    e += fib
        full.append(f); etale.append(e)
    return twists, full, etale, U, vectors


def incidence_constant(names):
    """For a 7-line subset: t3 and the polynomial constant 15 - t3 (so #U = q^2 - 6q + 15 - t3 away from 2, 3)."""
    from math import gcd
    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    def normalize(v):
        g = 0
        for c in v:
            g = gcd(g, abs(c))
        v = tuple(c // g for c in v)
        for c in v:
            if c:
                return v if c > 0 else tuple(-x for x in v)
    pts = {}
    for a, b in combinations(names, 2):
        pts.setdefault(normalize(cross(MAGIC[a], MAGIC[b])), set()).update({a, b})
    mult = [len(L) for L in pts.values()]
    need(max(mult) <= 3, 'no four concurrent lines')
    t3 = mult.count(3); t2 = mult.count(2)
    need(t2 + 3 * t3 == len(names) * (len(names) - 1) // 2, 'pair count')
    return t2, t3


def main():
    out = {'fields': []}
    fields = [GF(p, 1) for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)] + [GF(5, 2), GF(7, 2)] + [GF(3, 1), GF(3, 2)]
    for F in fields:
        q = F.q; rec = {'q': q, 'char': F.p}
        P2 = q * q + q + 1
        for label, names in (('hourglass', HOURGLASS), ('k3six', K3SIX), ('nine', NINE)):
            twists, full, etale, U, vectors = twist_counts(F, names)
            n = len(names)
            need(sum(full) == 2 ** (n - 1) * P2, '%s q=%d: twist sum of full counts = 2^(n-1) #P^2' % (label, q))
            need(sum(etale) == 2 ** (n - 1) * len(U), '%s q=%d: twist sum of etale counts = 2^(n-1) #U' % (label, q))
            rec[label] = {'n': n, 'U': len(U), 'twist_classes': len(twists), 'min_full': min(full), 'max_full': max(full),
                          'min_etale': min(etale), 'max_etale': max(etale)}
            if label == 'hourglass':
                if F.p >= 5:
                    need(len(U) == q * q - 6 * q + 12, 'hourglass complement count q^2-6q+12 at q=%d' % q)
                elif F.p == 3:
                    need(len(U) == q * q - 6 * q + 10, 'hourglass complement count q^2-6q+10 in characteristic 3 at q=%d' % q)
                # Walsh decomposition: N_t^{etale} = sum_{J even} chi_J(t) S_J,  S_J = sum_{y in U} prod_{j in J} chi(L_j(y))
                evenJ = [J for r in range(0, n + 1, 2) for J in combinations(range(n), r)]
                S = {}
                for J in evenJ:
                    s = 0
                    for v in U:
                        pr = 1
                        for j in J:
                            pr *= v[j]
                        s += pr
                    S[J] = s
                for eps, e in zip(twists, etale):
                    walsh = sum(((-1) ** sum(eps[j] for j in J)) * S[J] for J in evenJ)
                    need(walsh == e, 'Walsh identity at q=%d' % q)
                need(S[()] == len(U), 'S_empty = #U')
                six = {tuple(names[j] for j in range(n) if j not in J): S[J] for J in evenJ if len(J) == 6}
                outer = {S[J] for J in evenJ if len(J) == 6 and 0 in J}
                central = S[tuple(range(1, 7))]
                rec['hourglass']['S_six_outer_deletions_distinct_values'] = sorted(outer)
                rec['hourglass']['S_six_central_deletion'] = central
                rec['hourglass']['parallel_pair_sums'] = [S[(1, 2)], S[(3, 4)], S[(5, 6)]]
                rec['hourglass']['max_abs_S_by_size'] = {r: max(abs(S[J]) for J in evenJ if len(J) == r) for r in (2, 4, 6)}
                need(len(outer) == 1, 'the six outer deletions have equal character sums (1+6 split)')
        out['fields'].append(rec)
    # the seven-subset hostile and the (q-1)-expansion
    table = {}
    for sub in combinations(NINE, 7):
        t2, t3 = incidence_constant(list(sub))
        chi = 3 - 14 + t2 + 2 * t3
        need(chi == 10 - t3 and 15 - t3 == chi + 5, 'chi(U) = 10 - t3')
        table.setdefault(t3, 0); table[t3] += 1
    need(table == {2: 4, 3: 18, 4: 12, 5: 2}, 'seven-subset triple-point census (Session59 control)')
    t2, t3 = incidence_constant(HOURGLASS)
    need((t2, t3) == (12, 3), 'hourglass incidence 12 doubles, 3 triples')
    out['seven_subsets'] = {'t3_census': table, 'chi_U': {t3: 10 - t3 for t3 in table}, 'p_g': 7,
                            'chi_equals_pg_only_for_t3': [t3 for t3 in table if 10 - t3 == 7]}
    # q^2 - 6q + 12 = (q-1)^2 - 4(q-1) + 7 ; value 7 at q = 1 = chi of the complex complement
    expansion = [12 - 6 + 1, -6 + 2, 1]
    need(expansion == [7, -4, 1], '(q-1)-expansion coefficients 7, -4, 1')
    # Euler-number stratification of the resolved hourglass (Session59 audit): 2^6*7 + 2^5*(2n-2t2-3t3) + 2^4 t2 + 2^3 t3 + nodes
    e_strata = 64 * 7 + 32 * (14 - 24 - 9) + 16 * 12 + 8 * 3 + 24
    need(e_strata == 80 == 12 * (1 + 7) - 16, 'e(Y) = 80 by strata and by Noether')
    out['expansion_in_q_minus_1'] = expansion
    out['euler_stratification'] = {'etale_sheets_times_chiU': 64 * 7, 'total': e_strata}
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == '__main__':
    main()
