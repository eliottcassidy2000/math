"""Actual Hamilton endpoint responses, Smith lifts and the 81/27 control."""
import argparse
from collections import Counter
from itertools import permutations, combinations
import json
from math import isqrt, prod, gcd
from pathlib import Path
from precision import smith, replay, product


def need(ok, why):
    if not ok:
        raise ValueError(why)


def hostile(flip):
    # Exact constructor inherited from root_ports_closed_carriers_20261028/
    # check_carriers.py:build_hostile, including the original labels.
    rows = [0]*12
    for u in range(12):
        for v in range(u+1, 12):
            b, q = divmod(u, 4)
            c, r = divmod(v, 4)
            arc = (c-b) % 3 == 1 if b != c else q == 0 or q > 0 and (r-q) % 3 == 1
            if b == c == 2 and flip:
                arc = not arc
            rows[u if arc else v] |= 1 << (v if arc else u)
    return rows


def endpoints(rows):
    n, H = len(rows), []
    full = (1 << n)-1
    for s in range(n):
        dp = [[0]*n for _ in range(full+1)]
        dp[1 << s][s] = 1
        for mask in range(1 << s, full+1):
            for u, value in enumerate(dp[mask]):
                if not value:
                    continue
                todo = rows[u] & ~mask
                while todo:
                    bit = todo & -todo
                    todo -= bit
                    dp[mask | bit][bit.bit_length()-1] += value
        H.append(dp[full])
    return H


def walk_inclusion_exclusion(rows, source, target):
    n, total = len(rows), 0
    internal = [i for i in range(n) if i not in (source, target)]
    for bits in range(1 << len(internal)):
        mask = (1 << source) | (1 << target)
        for j, v in enumerate(internal):
            if bits >> j & 1:
                mask |= 1 << v
        vector = [int(i == source) for i in range(n)]
        for _ in range(n-1):
            new = [0]*n
            for u, value in enumerate(vector):
                todo = rows[u] & mask if value else 0
                while todo:
                    bit = todo & -todo
                    todo -= bit
                    new[bit.bit_length()-1] += value
            vector = new
        total += (-1)**(len(internal)-bits.bit_count()) * vector[target]
    return total


def automorphisms(rows):
    n, degree = len(rows), [r.bit_count() for r in rows]
    out, current = [], []
    def visit():
        u = len(current)
        if u == n:
            out.append(current[:])
            return
        for v in range(n):
            if v in current or degree[u] != degree[v]:
                continue
            if all((rows[u] >> j & 1) == (rows[v] >> current[j] & 1) for j in range(u)):
                current.append(v)
                visit()
                current.pop()
    visit()
    return out


def valuation(n, p):
    need(n != 0, 'finite invariant valuation')
    v = 0
    while n % p == 0:
        v += 1
        n //= p
    return v


def determinant(A):
    M = [r[:] for r in A]
    sign = previous = 1
    for k in range(len(M)-1):
        pivot = next((i for i in range(k, len(M)) if M[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            M[k], M[pivot] = M[pivot], M[k]
            sign = -sign
        p = M[k][k]
        for i in range(k+1, len(M)):
            for j in range(k+1, len(M)):
                value = M[i][j]*p-M[i][k]*M[k][j]
                need(value % previous == 0, 'exact Bareiss division')
                M[i][j] = value // previous
        for i in range(k+1, len(M)):
            M[i][k] = 0
        previous = p
    return sign*M[-1][-1]


def certify(label, rows, paley=False):
    n = len(rows)
    H = endpoints(rows)
    if n <= 7:
        brute = [[0]*n for _ in range(n)]
        for path in permutations(range(n)):
            if all(rows[a] >> b & 1 for a, b in zip(path, path[1:])):
                brute[path[0]][path[-1]] += 1
        need(H == brute, 'complete permutation census')
    probes = [(0, 1), (0, n-1)]
    if n == 12:
        probes.extend([(8, 9), (9, 10)])
    for s, t in probes:
        need(walk_inclusion_exclusion(rows, s, t) == H[s][t], 'walk inclusion exclusion')
    D, U, V, ops = smith(H)
    need(replay(H, ops) == D and product(product(U, H), V) == D, 'Smith certificate')
    diagonal = [D[i][i] for i in range(n)]
    if n > 7:
        # The elementary producer has severe coefficient swell. Obtain a
        # smaller witness with SymPy, then verify it without library algebra.
        from sympy import ZZ
        from sympy.polys.matrices import DomainMatrix
        from sympy.polys.matrices.normalforms import smith_normal_decomp
        dm = DomainMatrix([[ZZ(x) for x in row] for row in H], (n, n), ZZ)
        sd, su, sv = smith_normal_decomp(dm)
        D2, U, V = [[[int(x) for x in row] for row in mat.to_list()] for mat in (sd, su, sv)]
        need(D2 == D, 'independent Smith producer')
        need(product(product(U, H), V) == D, 'compact certificate multiplication')
        need(abs(determinant(U)) == abs(determinant(V)) == 1, 'compact witness unimodular')
        ops = []
    group = automorphisms(rows)
    for g in group:
        need(all(H[g[i]][g[j]] == H[i][j] for i in range(n) for j in range(n)), 'endpoint transport')
    out = dict(label=label, rows=rows, endpoints=H, total_paths=sum(map(sum, H)),
               smith_diagonal=diagonal, smith_operations=ops, row_transform=U,
               column_transform=V, full_automorphism_order=len(group),
               two_adic_exponents=[valuation(d, 2) for d in diagonal],
               three_adic_exponents=[valuation(d, 3) for d in diagonal])
    if paley:
        xvals = {H[i][j] for i in range(n) for j in range(n) if rows[i] >> j & 1}
        yvals = {H[i][j] for i in range(n) for j in range(n) if i != j and not(rows[i] >> j & 1)}
        need(len(xvals) == len(yvals) == 1, 'Paley two-orbit endpoint law')
        x, y = next(iter(xvals)), next(iter(yvals))
        norm = ((x+y)**2+n*(x-y)**2)//4
        predicted_determinant = ((n-1)//2)*(x+y)*norm**((n-1)//2)
        need(prod(diagonal) == predicted_determinant, 'Paley norm determinant')
        r, s = (n-1)//2, x+y
        need(gcd(x, y) == 1 and s % n != 0, 'conditional module theorem hypotheses')
        d = gcd(r*s, norm)
        formula = [1]*r+[d]+[norm]*(r-1)+[r*s*norm//d]
        need(diagonal == formula, 'full Paley Smith formula')
        out.update(forward_endpoints=x, backward_endpoints=y, quadratic_norm=norm)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    objects = []
    for n in (3, 7, 11):
        residues = sorted({a*a % n for a in range(1, n)})
        rows = [sum(1 << ((u+d) % n) for d in residues) for u in range(n)]
        objects.append(certify(f'Paley{n}', rows, True))
    T5 = [sum(1 << ((u+d) % 5) for d in (1, 2)) for u in range(5)]
    objects.append(certify('T5', T5))
    left, right = certify('81-control', hostile(False)), certify('27-control', hostile(True))
    def triples(rows):
        return [list(vs) for vs in combinations(range(12), 3)
                if all(sum(rows[u] >> v & 1 for v in vs if u != v) == 1 for u in vs)]
    need(triples(left['rows']) == triples(right['rows']), 'same labelled cyclic triples')
    need([left['full_automorphism_order'], right['full_automorphism_order']] == [81, 27], 'full transporter control')
    need(left['total_paths'] == right['total_paths'] == 462915, 'same total Hamilton paths')
    need(all(v % 9 == 0 for row in left['endpoints'] for v in row), 'left vanishes modulo nine')
    witness = next((i, j, right['endpoints'][i][j]) for i in range(12) for j in range(12)
                   if right['endpoints'][i][j] % 9)
    need(all(v % 3 == 0 for row in right['endpoints'] for v in row), 'right vanishes modulo three')
    objects.extend([left, right])
    H7 = objects[1]['endpoints']
    code = {0}
    for row in H7:
        r = sum((v % 2) << i for i, v in enumerate(row))
        code |= {w ^ r for w in list(code)}
    need(Counter(w.bit_count() for w in code) == {0: 1, 3: 7, 4: 7, 7: 1}, 'Hamming weight enumerator')
    primes = [7, 2, 19, 5, 59, 29, 89]
    square_selections = {w for w in range(128)
                         if (lambda z: isqrt(z)**2 == z)(prod(p+1 for i, p in enumerate(primes) if w >> i & 1))}
    need(square_selections == code, 'explicit prime-bank isomorphism')
    # Larger exact counts are inherited from Session54, not rerun here.
    q, x, y = 19, 3040229428, 3817640437
    r, s, delta = 9, x+y, x-y
    norm = (s*s+q*delta*delta)//4
    d = gcd(r*s, norm)
    formula = [1]*r+[d]+[norm]*(r-1)+[r*s*norm//d]
    need(gcd(x, y) == 1 and s % q, 'inherited q19 hypotheses')
    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_form
    residues = {a*a % q for a in range(1, q)}
    H19 = [[0 if i == j else x if (j-i) % q in residues else y for j in range(q)] for i in range(q)]
    S19 = smith_normal_form(Matrix(H19), domain=ZZ)
    need([abs(int(S19[i, i])) for i in range(q)] == formula, 'independent q19 Smith computation')
    result = dict(objects=objects, hostile_mod_nine_witness=witness,
                  prime_labels=primes, hamming_code=sorted(code),
                  labelled_cyclic_triples=triples(left['rows']),
                  inherited_q19=dict(a=x, b=y, norm=norm, smith_diagonal=formula,
                    source='../session54_reframed_numerology_20261025/paley-endpoint-observer/RESULT.md'))
    path = Path(__file__).with_suffix('.json')
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    else:
        need(json.loads(path.read_text(encoding='utf-8')) == json.loads(json.dumps(result)), 'stored atom receipt')
    print(json.dumps(dict(objects=[{k: d[k] for k in ['label','total_paths','smith_diagonal',
                     'full_automorphism_order']} for d in objects], witness=witness,
                     paley_norms=[d.get('quadratic_norm') for d in objects[:3]])))


if __name__ == '__main__':
    main()
