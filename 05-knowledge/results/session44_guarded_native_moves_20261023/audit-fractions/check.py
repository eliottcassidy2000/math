#!/usr/bin/env python3
"""Independent fractions audit: arc-set canon, SymPy Berkowitz, full phase scans.

Does not import producer code. SymPy 1.14.0 is the installed independent
symbolic path; exact arithmetic only, and all guards survive python -O.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, gcd
from pathlib import Path
import argparse
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
COUNT = 0
X = sp.Symbol('x')


def need(ok, why):
    global COUNT
    COUNT += 1
    if not ok:
        raise RuntimeError(why)


@lru_cache(None)
def edge_pairs(n):
    return tuple(combinations(range(n), 2))


@lru_cache(None)
def graph_label(n, edges):
    position = {e: i for i, e in enumerate(edge_pairs(n))}
    masks = []
    for sigma in permutations(range(n)):
        renamed = {tuple(sorted((sigma[u], sigma[v]))) for u, v in edges}
        masks.append(sum(1 << position[e] for e in renamed))
    return min(masks)


@lru_cache(None)
def tournament_label(n, arcs):
    position = {e: i for i, e in enumerate(edge_pairs(n))}
    masks = []
    for sigma in permutations(range(n)):
        backwards = {(sigma[v], sigma[u]) for u, v in arcs if sigma[u] > sigma[v]}
        masks.append(sum(1 << position[e] for e in backwards))
    return min(masks)


def own_records(n):
    out = []
    short_edges = edge_pairs(n-1)
    for code in range(2**len(short_edges)):
        compressed = frozenset(e for i, e in enumerate(short_edges) if code & (1 << i))
        backwards = frozenset((u, v+1) for u, v in compressed)
        arcs = frozenset((v, u) if (u, v) in backwards else (u, v) for u, v in edge_pairs(n))
        need(all((i, i+1) in arcs for i in range(n-1)), 'independent native forced arcs')
        native_mask = sum(1 << i for i, e in enumerate(edge_pairs(n)) if e in backwards)
        out.append((code, native_mask, graph_label(n, backwards), graph_label(n-1, compressed), tournament_label(n, arcs)))
    return out


def incidence(records, key):
    rows = sorted({key(row) for row in records})
    columns = sorted({row[-1] for row in records})
    entries = Counter((key(row), row[-1]) for row in records)
    W = sp.Matrix([[entries[g, t] for t in columns] for g in rows])
    return rows, columns, W


def coupling(W):
    rows = list(W*sp.ones(W.cols, 1))
    cols = list(W.T*sp.ones(W.rows, 1))
    K = sp.diag(*(1/x for x in rows))*W*sp.diag(*(1/x for x in cols))*W.T
    return rows, cols, K


def as_fractions(values):
    return [sp.Rational(value) for value in values]


FULL = [34479695250,-101537044560,124233482025,-81343863297,30498542418,-6453737759,695751836,-27925620]
COMP = [3377602800,-8397030240,8429380020,-4360639860,1225690398,-180039440,12012273,-267391]
PAIR = [562933800,-914443020,553227084,-153488811,19184700,-874496]


def spectral_audit(producer):
    evidence = json.loads((producer/'receipt.json').read_text())
    matrices = {}
    for n in (4, 5):
        records = own_records(n)
        need([list(r) for r in records] == evidence['spectra'][str(n)]['word_records'], 'all native record classes independently rebuilt')
        for name, key in [('full', lambda r:r[2]), ('compressed', lambda r:r[3]), ('refinement', lambda r:(r[2], r[3]))]:
            labels, targets, W = incidence(records[1:], key)
            old = evidence['spectra'][str(n)]['cases'][name]
            need(W.tolist() == [as_fractions(row) for row in old['incidence']], 'independently reconstructed incidence '+name)
            rows, cols, K = coupling(W)
            # Berkowitz determinant algorithm, independent of producer Newton traces.
            coefficients = K.berkowitz_charpoly(X).all_coeffs()
            need(coefficients == as_fractions(old['charpoly_descending']), 'independent Berkowitz characteristic '+name)
            need(W.rank() == old['rank'] and sp.trace(K) == sp.Rational(old['trace']), 'independent rank and exact trace')
            matrices[n, name] = (W, K)
        paired = []
        for row in records:
            compressed = frozenset(e for i,e in enumerate(edge_pairs(n-1)) if row[0] & (1 << i))
            opposite = frozenset(set(edge_pairs(n-1))-set(compressed))
            paired.append(row[:3]+(min(row[3], graph_label(n-1, opposite)),)+row[4:])
        _, _, W = incidence(paired, lambda r:r[3])
        old = evidence['spectra'][str(n)]['cases']['paired_compressed']
        need(W.tolist() == [as_fractions(row) for row in old['incidence']], 'independent complement-paired incidence')
        need(coupling(W)[2].berkowitz_charpoly(X).all_coeffs() == as_fractions(old['charpoly_descending']), 'independent paired Berkowitz characteristic')
    roots = {}
    for name, coefficients, p, lo, hi in [('full',FULL,31,635789,635790), ('compressed',COMP,23,642470,642471), ('paired',PAIR,19,573377,573378)]:
        f = sp.Poly.from_list(coefficients, X)
        field = sp.Poly(f, X, modulus=p)
        need(field.is_irreducible, 'independent finite-field factorization '+name)
        need(f.is_irreducible, 'independent rational factorization '+name)
        a, b = sp.Rational(lo, 10**6), sp.Rational(hi, 10**6)
        need(f.count_roots(0,1) == f.degree(), 'all roots real and strictly between zero and one')
        need(f.count_roots(a,b) == 1 and f.count_roots(b,1) == 0, 'largest-root exact isolation')
        roots[name] = dict(degree=f.degree(), modulus=p, interval=[str(a),str(b)])
    need(sp.Rational(635790,10**6) < sp.Rational(642470,10**6), 'native compressed gap strictly smaller')
    return dict(order4_words=8, order5_words=64, independently_rebuilt_incidence_matrices=8,
                independent_spectral_method='SymPy Berkowitz characteristic polynomial, exact factorization and root count',
                largest_residual_roots=roots)


def compression_audit():
    tested = 0
    profile_zero, profile_positive = 0, 0
    for entries in product(range(3), repeat=6):
        W = sp.Matrix(3, 2, entries)
        r = list(W*sp.ones(2,1)); c = list(W.T*sp.ones(3,1))
        if not all(r) or not all(c):
            continue
        J = sp.Matrix([[1,0],[1,0],[0,1]])
        E = sp.Matrix([[r[0]/(r[0]+r[1]),r[1]/(r[0]+r[1]),0],[0,0,1]])
        K = coupling(W)[2]
        Wbar = J.T*W
        rr, _, Kbar = coupling(Wbar)
        need(E*J == sp.eye(2) and Kbar == E*K*J, 'weighted compression identity')
        rhs = sum(r[i]*sum((W[i,t]/r[i]-Wbar[0,t]/rr[0])**2/c[t] for t in range(2)) for i in (0,1))
        need(sp.trace(K)-sp.trace(Kbar) == rhs, 'independent exact profile variance formula')
        equality = all(W[0,t]/r[0] == W[1,t]/r[1] for t in range(2))
        factored = K == J*Kbar*E
        finepoly = K.berkowitz_charpoly(X).as_expr()
        coarsepoly = Kbar.berkowitz_charpoly(X).as_expr()
        only_zeros = sp.expand(finepoly-X*coarsepoly) == 0
        need(equality == (rhs == 0) == factored == only_zeros, 'all four equivalences, not merely equal rank')
        profile_zero += equality
        profile_positive += not equality
        tested += 1
    return dict(nonnegative_3_by_2_tables=tested, zero_defect_tables=profile_zero,
                positive_defect_tables=profile_positive, entry_alphabet=[0,1,2],
                source_partition=[[0,1],[2]])


def converse_audit():
    cases = 0
    for n in range(2,7):
        compressed_edges = edge_pairs(n-1)
        for code in range(2**len(compressed_edges)):
            compressed = {e for i,e in enumerate(compressed_edges) if code & (1 << i)}
            backwards = {(u,v+1) for u,v in compressed}
            arcs = {(v,u) if (u,v) in backwards else (u,v) for u,v in edge_pairs(n)}
            changed_arcs = {(n-1-v,n-1-u) for u,v in arcs}
            reflected_backwards = {(n-1-v,n-1-u) for u,v in backwards}
            reconstructed = {(v,u) if (u,v) in reflected_backwards else (u,v) for u,v in edge_pairs(n)}
            need(changed_arcs == reconstructed, 'literal reflection plus converse preserves fixed-path form')
            compressed_reflection = {(n-2-v,n-2-u) for u,v in compressed}
            need({(u,v+1) for u,v in compressed_reflection} == reflected_backwards, 'compressed reflection transports actual chord positions')
            cases += 1
    return dict(full_word_reflection_cases=cases, native_orders=[2,3,4,5,6],
                weights_filter_boundary='Invariance under the actual reflected-converse bijection is required')


def continuation_audit():
    rows = []
    # Scan every phase tuple, not just producer-generated binary differences.
    for t in (3,5,7):
        arcs = set()
        vertices = list(product(range(t),range(3)))
        for x,y in combinations(vertices,2):
            i,u=x; j,v=y
            if i == j:
                forward=(v-u)%3 == 1
            elif 1 <= (j-i)%t <= (t-1)//2:
                forward=(v-u)%3 in (1,2)
            else:
                forward=(u-v)%3 not in (1,2)
            arcs.add((x,y) if forward else (y,x))
        need(all(sum((x,y) in arcs for y in vertices if x != y) == (3*t-1)//2 for x in vertices), 'independent actual regular host')
        counts, first, accepted = Counter(), Counter(), 0
        for phases in product(range(3), repeat=t):
            P = [(i,(phases[i]+j)%3) for i in range(t) for j in range(3)]
            valid = all((u,v) in arcs for u,v in zip(P,P[1:]))
            differences = tuple((phases[i+1]-phases[i])%3 for i in range(t-1))
            need(valid == all(d in (0,1) for d in differences), 'full phase scan exact native validity iff binary differences')
            if valid:
                a,d = phases[0],(phases[-1]-phases[0])%3
                counts[a,d] += 1; first[a,d,differences[0]] += 1
                accepted += 1
        for a in range(3):
            for d in range(3):
                expected = sum(comb(t-1,k) for k in range(d,t,3))
                need(counts[a,d] == expected, 'independent pinned endpoint counts')
                for bit in (0,1):
                    child = sum(comb(t-2,k) for k in range((d-bit)%3,t-1,3))
                    need(first[a,d,bit] == child, 'independent conditional native first-step counts')
        rows.append(dict(cells=t,phase_tuples=3**t,accepted=accepted))
    Js = [0,1]
    for r in range(2,122): Js.append(Js[-1]+2*Js[-2])
    count = 0
    for r in range(3,121):
        for d in range(3):
            a = sum(comb(r-1,k) for k in range(d,r,3))
            b = sum(comb(r-1,k) for k in range((d-1)%3,r,3))
            pair = tuple(sorted((Fraction(a,a+b),Fraction(b,a+b))))
            J = Js[r]
            expected = (Fraction(1,2),)*2 if a == b else (Fraction(J-1,2*J),Fraction(J+1,2*J))
            need(pair == expected and min(a,b)>0, 'all endpoint residues conditional Jacobsthal formula')
            need(Fraction(7,15) not in pair and Fraction(8,15) not in pair, 'native continuation is not old spectral fraction pair')
            if a != b:
                p,q = pair[0].numerator,pair[0].denominator
                need(q == J and q-2*p == 1, 'lowest-term unequal pair and half-neighbor determinant')
                need(abs(p*pair[1].denominator-pair[1].numerator*q) == J, 'pair determinant is J, not one')
                depth=(J-3)//2
                need(Fraction(depth+1,2*depth+3) == pair[0], 'closed A-spine address, all tested clocks')
            count += 1
        need(Js[r+1]*Js[r-1]-Js[r]**2 == (-1)**r*2**(r-1), 'independent generalized continued fraction determinant')
        need(gcd(Js[r],Js[r+1]) == 1, 'primitive ratio approximants')
    # A acts on k/(2k+1) by k -> k+1; direct algebra is an all-k proof.
    k=sp.Symbol('k')
    z=k/(2*k+1)
    need(sp.cancel((1-z)/(3-4*z)-(k+1)/(2*k+3)) == 0, 'symbolic all-k Berggren A action')
    return dict(actual_full_phase_scans=rows,conditional_cases=count,clock_range=[3,120],
                arbitrary_clock_interpretation='Remaining continuation suffix inside an odd-cell regular host',
                generalized_transfer_determinant=-2)


def main():
    parser=argparse.ArgumentParser()
    default_producer = HERE.parent/'fractions-ternary'
    if not default_producer.exists():
        default_producer = HERE.parent/'session44-fractions-ternary'
    parser.add_argument('--producer',type=Path,default=default_producer)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    pins=json.loads((HERE/'source_pins.json').read_text())
    for record in pins['files']:
        need(sha256((args.producer/record['path']).read_bytes()).hexdigest() == record['sha256'], 'producer source pin '+record['path'])
    out=dict(status='PASS: independent audit of scoped fractions claims',sympy_version=sp.__version__,
             spectra=spectral_audit(args.producer), compression=compression_audit(),
             converse=converse_audit(),native=continuation_audit())
    out['exact_check_count']=COUNT
    text=json.dumps(out,indent=2,sort_keys=True)+'\n'
    target=HERE/'receipt.json'
    if args.write: target.write_text(text)
    else: need(target.read_text() == text, 'independent frozen receipt')
    print(text,end='')


if __name__=='__main__':
    main()
