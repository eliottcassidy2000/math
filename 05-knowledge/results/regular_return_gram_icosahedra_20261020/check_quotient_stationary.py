"""Exact checks of the path-first walk and its ordinary-class stationary mass.
Independent arc-set implementation; no producer checker imports.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import json
from independent_audit import tournament, code_of, degree, complement_path, require


def hpaths(arcs, n, start=None, end=None):
    outs = {u:{v for x,v in arcs if x == u} for u in range(n)}
    def extend(p, unused):
        if not unused:
            if end is None or p[-1] == end:
                yield p
            return
        for v in sorted(outs[p[-1]] & unused):
            if v == end and len(unused) > 1:
                continue
            yield from extend(p+(v,), unused-{v})
    for s in range(n) if start is None else (start,):
        yield from extend((s,), set(range(n))-{s})


def relabel(arcs, pi):
    return frozenset((pi[u],pi[v]) for u,v in arcs)


def reverse(arcs):
    return frozenset((v,u) for u,v in arcs)


def regular_row(arcs, n):
    p_words = list(hpaths(arcs,n))
    row = defaultdict(Fraction)
    columns = {}
    for p in p_words:
        s = complement_path(arcs,p)
        sid = code_of(s,n)
        require(sid not in columns, 'binary source incidence')
        qs = list(hpaths(s,n,p[0],p[-1]))
        require(qs, 'positive column degree')
        targets = set()
        for q in qs:
            u = complement_path(s,q)
            require(degree(u,n) == ((n-1)//2,)*n, 'regular return')
            uid = code_of(u,n)
            require(uid not in targets, 'binary target incidence')
            targets.add(uid)
            row[uid] += Fraction(1,len(p_words)*len(qs))
        columns[sid] = frozenset(targets)
    require(sum(row.values()) == 1, 'row-stochastic')
    return dict(row), columns, len(p_words)


def ordinary_orbit(arcs,n):
    return {code_of(relabel(arcs,pi),n) for pi in permutations(range(n))}


def main():
    small=[]
    for n in (1,3,5):
        bank=[tournament(n,c) for c in range(1 << (n*(n-1)//2))]
        bank=[a for a in bank if degree(a,n) == ((n-1)//2,)*n]
        matrices={code_of(a,n):regular_row(a,n) for a in bank}
        for t,(row,cols,h) in matrices.items():
            for u,(other,_,h_u) in matrices.items():
                require(h*row.get(u,0) == h_u*other.get(t,0),'detailed balance')
        orbit=ordinary_orbit(bank[0],n)
        require(orbit == set(matrices),'one ordinary class')
        small.append(dict(n=n,states=len(bank),H=matrices[min(matrices)][2],
                          detailed_balance_entries=len(bank)**2))

    rows_bank=[(26,28,41,112,100,67,7),(14,28,56,112,97,67,7),
               tuple(sum(1 << ((i+d)%7) for d in (1,2,4)) for i in range(7))]
    representatives=[frozenset((u,v) for u in range(7) for v in range(7) if row[u] & (1 << v)) for row in rows_bank]
    orbits=[ordinary_orbit(a,7) for a in representatives]
    lookup={}
    for k,orb in enumerate(orbits):
        for t in orb:
            require(t not in lookup,'disjoint representative orbits')
            lookup[t]=k
    require(list(map(len,orbits)) == [1680,720,240],'ordinary orbit sizes')
    Q=[]; data=[]
    generators=[tuple([1,0,2,3,4,5,6]),tuple([1,2,3,4,5,6,0])]
    for k,a in enumerate(representatives):
        row,columns,h=regular_row(a,7)
        quotient=[sum(p for t,p in row.items() if lookup[t] == j) for j in range(3)]
        require(sum(quotient) == 1,'quotient row sum')
        Q.append(quotient)
        for pi in generators:
            ar=relabel(a,pi)
            r2,c2,h2=regular_row(ar,7)
            expected={code_of(relabel(tournament(7,t),pi),7):p for t,p in row.items()}
            require(r2 == expected and h2 == h,'relabel covariance')
            require([sum(p for t,p in r2.items() if lookup[t] == j) for j in range(3)] == quotient,'ordinary lumpability')
        converse=reverse(a)
        r2,c2,h2=regular_row(converse,7)
        expected={code_of(reverse(tournament(7,t)),7):p for t,p in row.items()}
        require(r2 == expected and h2 == h,'converse covariance')
        self_converse=code_of(converse,7) in orbits[k]
        require(self_converse,'known regular7 class self-converse')
        data.append(dict(H=h,labelled_orbit_size=len(orbits[k]),aut_order=5040//len(orbits[k]),
                         path_orbit_mass=h*len(orbits[k])//5040,
                         intermediate_degree_distribution=dict(sorted(Counter(map(len,columns.values())).items())),
                         row_support=len(row),self_converse=self_converse))
    h_orbit=[Fraction(d['H']*d['labelled_orbit_size'],5040) for d in data]
    require(h_orbit == [57,25,9],'class stationary weights')
    for i in range(3):
        for j in range(3):
            require(h_orbit[i]*Q[i][j] == h_orbit[j]*Q[j][i],'quotient detailed balance')
    stationary=[x/sum(h_orbit) for x in h_orbit]
    require([sum(stationary[i]*Q[i][j] for i in range(3)) for j in range(3)] == stationary,'quotient stationary vector')
    result=dict(status='PASS',small_orders=small,
                regular7=dict(class_order_H=[d['H'] for d in data],classes=data,
                              quotient_kernel=[[str(x) for x in row] for row in Q],
                              stationary_class_weights=[str(x) for x in h_orbit],
                              stationary_class_distribution=[str(x) for x in stationary],
                              scope='Rows from the three inherited complete regular-order-seven classes; their full relabelling orbits independently rebuilt. No all-order class census.'),
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                dependency_sha256=sha256(Path(__file__).with_name('independent_audit.py').read_bytes()).hexdigest())
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
