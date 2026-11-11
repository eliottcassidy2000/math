"""Exact finite controls for Session259 diagonal growth; no external packages."""
from collections import Counter
from itertools import product, combinations, combinations_with_replacement
from math import prod, isqrt


def require(p, label):
    if not p:
        raise RuntimeError(label)


def conv(a, b, op):
    c = Counter()
    for x, nx in a.items():
        for y, ny in b.items():
            c[op((x, y))] += nx * ny
    return c


def enum(S, k, op, kind):
    it = {'ordered': product, 'weak': combinations_with_replacement,
          'strict': combinations}[kind]
    return Counter(op(t) for t in (it(S, repeat=k) if kind == 'ordered' else it(S, k)))


def run():
    pairs = 0
    for N in range(1, 81):
        J = {(a+b, a*b) for a in range(1, N+1) for b in range(a, N+1)}
        old = {(a+b, a*b) for a in range(1, N) for b in range(a, N)}
        require(J-old == {(N+a, N*a) for a in range(1, N+1)}, 'birth strip')
        for s, p in J-old:
            d = isqrt(s*s-4*p)
            a, b = (s-d)//2, (s+d)//2
            require((a+b, a*b, b) == (s, p, N), 'inverse and birth')
            for c in (2, 3):
                require((c*s)**2-4*c*c*p == c*c*d*d, 'scale discriminant')
                qs, qp = s*s-2*p, p*p
                require(qs*qs-4*qp == s*s*d*d, 'power discriminant')
                require(((c*s)**2-2*c*c*p, (c*c*p)**2) ==
                        (c*c*qs, c**4*qp), 'QD law')
            pairs += 1
    print('FINITE-EXACT: joint birth strips, inverse and power/scale laws, prefixes N=1..80; pairs=', pairs)
    checks = 0
    for op_name, op in [('sum', sum), ('product', prod)]:
        for N in range(1, 9):
            S, oldS = list(range(1, N+1)), list(range(1, N))
            for k in range(1, 5):
                for kind in ('weak', 'strict'):
                    direct = enum(S, k, op, kind)
                    old = enum(oldS, k, op, kind)
                    inc = Counter()
                    for j in (range(1, k+1) if kind == 'weak' else (1,)):
                        rest = enum(oldS, k-j, op, kind)
                        inc.update({op([N]*j+[v]): count for v, count in rest.items()})
                    require(direct == old+inc, (op_name, N, k, kind, 'insertion'))
                    checks += 1
            F = Counter(S)
            C = enum(S, 3, op, 'ordered')
            T = conv(Counter(op([a, a]) for a in S), F, op)
            D = Counter(op([a]*3) for a in S)
            H, E = enum(S, 3, op, 'weak'), enum(S, 3, op, 'strict')
            for z in C:
                require(6*H[z] == C[z]+3*T[z]+2*D[z], 'H3')
                require(6*E[z] == C[z]-3*T[z]+2*D[z], 'E3')
        for S in ([1, 2, 3, 4], [2, 3, 5], [1, 4, 9]):
            for ell in (2, 3, 5, 7):
                C = enum(S, ell, op, 'ordered')
                diag = Counter(op([a]*ell) for a in S)
                require(all(C[z] % ell == diag[z] % ell for z in C), 'prime cyclic law')
    print('FINITE-EXACT: binary/ternary/quaternary insertion identities, both operations, N=1..8; checks=', checks)
    print('FINITE-EXACT: ternary cycle formulas, both operations, N=1..8')
    print('FINITE-EXACT: prime-arity laws ell=2,3,5,7; S=[1,2,3,4], [2,3,5], [1,4,9]; both operations')
    hostile = enum([1, 2], 4, sum, 'ordered')[6]
    require(hostile == 6 and hostile % 4 == 2, 'composite-arity hostile')
    print('HOSTILE: arity 4, S={1,2}, sum target 6: ordered count 6, residue 2 mod 4, no all-equal witness')
    x, y = (1, 6, 6), (2, 2, 9)
    require(sum(x) == sum(y) == 13 and prod(x) == prod(y) == 36, 'triple chart hostile')
    print('HOSTILE: triples (1,6,6), (2,2,9): same sum=13/product=36, distinct births 6 and 9')
    print('All finite controls passed. All-N statements depend on proofs in RESEARCH.md.')

if __name__ == '__main__':
    run()
