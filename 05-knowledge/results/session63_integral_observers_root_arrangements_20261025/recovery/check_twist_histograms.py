#!/usr/bin/env python3
"""Independent finite controls, not an abstract geometric proof.

No producer imports. Prime and quadratic extension arithmetic is implemented
directly. Twist counts are computed independently by square-root tables and
compared with six-bit squareclass bins and their Walsh transform.
"""
import itertools
import json
from collections import Counter
from fractions import Fraction


def require(test, message):
    if not test:
        raise RuntimeError(message)


class Field:
    def __init__(self, p, extension=1):
        self.p, self.extension, self.q = p, extension, p**extension
        self.d = next((d for d in range(2, p) if pow(d, (p-1)//2, p) == p-1), None)

    def const(self, n):
        return n % self.p

    def add(self, x, y):
        p = self.p
        return (x % p + y % p) % p + p * ((x // p + y // p) % p)

    def neg(self, x):
        p = self.p
        return (-x % p) % p + p * ((- (x // p)) % p)

    def mul(self, x, y):
        p = self.p
        if self.extension == 1:
            return x*y % p
        a, b, c, d = x % p, x // p, y % p, y // p
        return (a*c+self.d*b*d) % p + p*((a*d+b*c) % p)

    def power(self, x, e):
        r = 1
        while e:
            if e & 1:
                r = self.mul(r, x)
            x = self.mul(x, x)
            e //= 2
        return r

    def inv(self, x):
        require(x != 0, 'division by zero')
        return self.power(x, self.q-2)

    def dot(self, x, y):
        r = 0
        for a, b in zip(x, y):
            r = self.add(r, self.mul(a, b))
        return r

    def normalize(self, v):
        scale = self.inv(next(x for x in v if x))
        return tuple(self.mul(x, scale) for x in v)


FORMS = [(1,0,0), (1,1,0), (1,-1,0), (1,0,1),
         (1,0,-1), (1,1,1), (1,-1,-1)]


def parity_sign(n):
    return -1 if n.bit_count() % 2 else 1


def check_field(p, extension=1):
    F = Field(p, extension)
    q = F.q
    elements = list(range(q))
    require(len({F.mul(a, a) for a in elements if a}) == (q-1)//2,
            'field square census')
    roots = Counter(F.mul(a, a) for a in elements)
    chi = {a: roots[a]-1 for a in elements}
    ns = next(a for a in elements if chi[a] == -1)
    require(all(F.mul(a, F.inv(a)) == 1 for a in elements if a), 'field inverse')
    forms = [tuple(F.const(a) for a in f) for f in FORMS]
    points = [(1,x,z) for x in elements for z in elements]
    points += [(0,1,z) for z in elements] + [(0,0,1)]
    require(len(points) == q*q+q+1, 'projective census')
    U, triple_points, incidence = [], [], Counter()
    evaluations = []
    for v in points:
        values = [F.dot(v, f) for f in forms]
        m = values.count(0)
        evaluations.append((values, m))
        if m >= 2:
            incidence[m] += 1
        if m == 3:
            triple_points.append(v)
        if m == 0:
            U.append(values)
    require(max(incidence) == 3, 'no fourfold points in tested odd fields')
    expected_incidence = {2:6, 3:5} if p == 3 else {2:12, 3:3}
    require(dict(incidence) == expected_incidence, 'incidence correction')
    expected_u = q*q-6*q+(10 if p == 3 else 12)
    require(len(U) == expected_u, 'branch complement count')

    bins = [0]*64
    for vals in U:
        require(vals[0] == 1, 'E chart normalization')
        mask = sum((chi[v] == -1) << i for i, v in enumerate(vals[1:]))
        bins[mask] += 1
    walsh = [sum(parity_sign(mask & s)*n for mask,n in enumerate(bins))
             for s in range(64)]
    open_counts, projective_counts, node_counts, resolved_counts = [], [], [], []
    for mask in range(64):
        coefficients = [1] + [ns if mask >> i & 1 else 1 for i in range(6)]
        open_n = projective_n = nodes_n = 0
        for vals, m in evaluations:
            j = next(i for i, x in enumerate(vals) if x)
            denominator = F.inv(F.mul(coefficients[j], vals[j]))
            n = 1
            for i in range(7):
                if i != j:
                    n *= roots[F.mul(F.mul(coefficients[i], vals[i]), denominator)]
            projective_n += n
            if m == 0:
                open_n += n
            if m == 3:
                nodes_n += n
        require(open_n == 64*bins[mask], 'root-table versus squareclass bins')
        require(open_n == sum(parity_sign(mask&s)*walsh[s] for s in range(64)),
                'Walsh inversion')
        open_counts.append(open_n)
        projective_counts.append(projective_n)
        node_counts.append(nodes_n)
        resolved_counts.append(projective_n + q*nodes_n)
    require(sum(open_counts) == 64*len(U), 'open mean')
    require(sum(projective_counts) == 64*(q*q+q+1), 'singular projective mean')
    require(sum(node_counts) == 64*len(triple_points), 'node mean')
    require(sum(resolved_counts) == 64*(q*q+(6 if p == 3 else 4)*q+1),
            'resolved mean')
    variance = Fraction(sum(n*n for n in open_counts),64) - len(U)**2
    require(variance == sum(a*a for a in walsh[1:]), 'Walsh variance')
    require(variance == 64*sum(n*n for n in bins)-len(U)**2, 'bin variance')
    covariances = []
    for h in range(64):
        cov = Fraction(sum(open_counts[m]*open_counts[m^h] for m in range(64)),64)-len(U)**2
        require(cov == sum(parity_sign(h&s)*walsh[s]**2 for s in range(1,64)),
                'shift covariance')
        covariances.append(int(cov))
    if q == 5:
        require(len(U) == 7 and bins.count(1) == 7 and bins.count(0) == 57,
                'seven occupied bins at F5')
    return dict(p=p, extension=extension, q=q, quadratic_parameter=F.d if extension==2 else None,
                nonsquare=ns, incidence=dict(incidence), complement=len(U),
                squareclass_bins=bins, twist_open_counts=open_counts,
                walsh_character_sums=walsh, mean_open=len(U), variance_open=int(variance),
                covariance_by_shift=covariances, twist_projective_counts=projective_counts,
                twist_node_counts=node_counts, twist_resolved_counts=resolved_counts,
                mean_projective=q*q+q+1, mean_nodes=len(triple_points),
                mean_resolved=q*q+(6 if p==3 else 4)*q+1,
                count_frequency=dict(sorted(Counter(open_counts).items())))


def main():
    cases=[check_field(p,e) for p,e in [(3,1),(5,1),(7,1),(3,2),(11,1),
                                       (13,1),(17,1),(5,2),(7,2)]]
    print(json.dumps(dict(status='PASS',
                         scope='Exact listed finite fields; resolved counts use ordinary A1 resolution addition q per rational node, not an independent geometric proof.',
                         cases=cases), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
