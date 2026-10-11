#!/usr/bin/env python3
"""Exact controls: changing generators, affine source offsets, atomic weights.

Finite universes are declared in main(); universal statements have separate
proofs in the companion notes. No universal ROOT or Jacobian claim is tested.
"""
import itertools
import json
from fractions import Fraction as Q
from math import comb

import sympy as s

x, t, v = s.symbols('x t v')
p = t+x*x*t*t
y = x*t*t+x**3*t**3
D = s.expand(p**3-y*y)
u = x*x*t
checks = 0


def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def lift(M):
    a, b, c, d = list(M)
    return (a*x+b)/(c*x+d), t*(c*x+d)**2


def pull(f, M):
    xx, tt = lift(M)
    return s.cancel(f.subs({x: xx, t: tt}, simultaneous=True))


def nu(r, m):
    return Q(1, (m+1)**r)-Q(1, (m+2)**r)


def replay(n, word):
    for a in word:
        z = 3*n+1
        if n == 1 or (z & -z).bit_length()-1 != a:
            raise RuntimeError('non-native or ROOT-padded receipt')
        n = z >> a
    return n


def linear_test(M, basis):
    """Independent coefficient equations in x,t for Phi(H)-K=Delta."""
    delta = s.expand((pull(u, M)-u)/2)
    cols = [s.Poly(pull(h, M), x, t) for h in basis]
    cols += [s.Poly(-h, x, t) for h in basis]
    target = s.Poly(delta, x, t)
    keys = sorted(set(target.monoms()).union(*(set(z.monoms()) for z in cols)))
    mat = s.Matrix([[z.coeff_monomial(k) for z in cols] for k in keys])
    rhs = s.Matrix([target.coeff_monomial(k) for k in keys])
    rank = mat.rank()
    return {'homogeneous_dimension': len(cols)-rank,
            'affine_solvable': rank == mat.row_join(rhs).rank()}


def main():
    mats = {
        'identity': s.eye(2), 'minus_identity': -s.eye(2),
        'upper': s.Matrix([[1, 1], [0, 1]]),
        'lower': s.Matrix([[1, 0], [-1, 1]]),
        'diagonal': s.diag(2, s.Rational(1, 2)),
        'anti': s.Matrix([[0, 1], [-1, 0]]),
        'scaled_anti': s.Matrix([[0, 2], [-s.Rational(1, 2), 0]]),
        'generic': s.Matrix([[1, 1], [1, 2]]),
        'affine': s.Matrix([[2, 1], [0, s.Rational(1, 2)]]),
        'a_zero': s.Matrix([[0, 1], [-1, 2]]),
    }
    for M in mats.values():
        xx, tt = lift(M)
        check(M.det() == 1, 'SL2 matrix')
        check(s.cancel(s.diff(xx, x)*s.diff(tt, t)
                       -s.diff(xx, t)*s.diff(tt, x)) == 1, 'symplectic lift')
        for N in mats.values():
            # Pullbacks reverse function-composition order.
            for f in [x, t]:
                check(s.cancel(pull(pull(f, M), N)-pull(f, M*N)) == 0,
                      'matrix composition retains cotangent multiplier')

    # Basis with exact highest t-degree <=9; different leading exponents
    # give an independent certificate of linear independence.
    basis, leading = [], []
    for a in range(5):
        for b in range(2):
            for c in range(2):
                L = 2*a+3*b+5*c
                if L <= 9:
                    h = s.expand(p**a*y**b*D**c)
                    basis.append(h)
                    pair = (L, L-c)
                    leading.append(pair)
                    check(s.Poly(h, t).LC() == x**(L-c), 'leading coefficient')
                    check(5*(L-c) >= 4*L, 'slope cone lower bound')
    check(len(set(leading)) == len(basis), 'no leading cancellation at fixed L')
    rows = {}
    for name, M in mats.items():
        row = linear_test(M, basis)
        identity = M == s.eye(2) or M == -s.eye(2)
        check(row['homogeneous_dimension'] == (len(basis) if identity else 1),
              'finite homogeneous intersection classification')
        check(row['affine_solvable'] == (name in {
            'identity', 'minus_identity', 'lower', 'anti', 'scaled_anti'}),
              'affine source compatibility classification')
        rows[name] = row

    for b in [s.Rational(1), s.Rational(2), s.Rational(-3), s.Rational(2, 3)]:
        U = s.Matrix([[1, b], [0, 1]])
        L = s.Matrix([[1, 0], [-1/b, 1]])
        W = U*L*U
        source = -(u+b*b*p)/2
        check(W == s.Matrix([[0, b], [-1/b, 0]]), 'three changed-generator flows')
        check(pull(source, W) == source, 'nonconstant affine-source fixed point')
        check(s.expand(pull(p, W)-p) != 0, 'fixed source does not fix p')
        check(pull(pull(x, W), W) == x and pull(pull(t, W), W) == t,
              'birational involution')
        # Independent rational chart with v=xt: cancellation of the pole.
        Hp = -b*b*(t+v*v)/2
        Hpull = -b*b*(v*v+v*v/(b*b*t))/2
        delta = (b*b*t-v*v/t)/2
        check(s.cancel(Hpull-Hp-delta) == 0, 'affine offset cancels t pole')

    # A larger positive family: all fibre-affine symplectic stabilizers of
    # I=u+p. Independent direct substitutions and Weierstrass chord formulas.
    invariant = u+p
    four = [(x, t), (-1/x, x*x*t),
            (1/x, -1-x*x*(1+t)), (-x, -1-t-1/x**2)]
    def sub(f, pair):
        return s.cancel(f.subs({x: pair[0], t: pair[1]}, simultaneous=True))
    for pair in four:
        xx, tt = pair
        check(sub(invariant, pair) == s.expand(invariant), 'four rational source stabilizers')
        check(s.cancel(s.diff(xx, x)*s.diff(tt, t)-s.diff(xx, t)*s.diff(tt, x)) == 1,
              'fibre-affine map is symplectic')
        for f in [x, t]:
            check(sub(sub(f, pair), pair) == f, 'each stabilizer squares to identity')
        for other in four:
            composite = tuple(sub(f, other) for f in pair)
            check(any(all(s.cancel(z-w) == 0 for z, w in zip(composite, candidate))
                      for candidate in four), 'Klein four closure')
            check(all(s.cancel(sub(f, other)-sub(g, pair)) == 0
                      for f, g in zip(pair, other)), 'Klein four commutativity')
    kappa, X, Y = s.symbols('kappa X Y')
    cubic = X*(X-1)*(X+kappa)
    translations = [(s.Integer(0), -kappa), (s.Integer(1), 1+kappa),
                    (-kappa, kappa*(kappa+1))]
    for root, B in translations:
        XP = root+B/(X-root)
        YP = -B*Y/(X-root)**2
        difference = s.together(YP**2-cubic.subs(X, XP))
        numerator = s.fraction(difference)[0]
        check(s.rem(s.Poly(numerator, Y), s.Poly(Y**2-cubic, Y)).as_expr() == 0,
              'two-torsion translation preserves the elliptic cubic')
        check(s.cancel((XP-root)*(X-root)-B) == 0, 'two-torsion chord abscissa')
    for pair, (root, B) in zip([four[2], four[1], four[3]], translations):
        oldX, oldY = -t, x*t*(t+1)
        expectedX = (root+B/(X-root)).subs({X: oldX, kappa: invariant}, simultaneous=True)
        expectedY = (-B*Y/(X-root)**2).subs(
            {X: oldX, Y: oldY, kappa: invariant}, simultaneous=True)
        check(s.cancel(sub(oldX, pair)-expectedX) == 0, 'ambient map induces torsion abscissa')
        check(s.cancel(sub(oldY, pair)-expectedY) == 0, 'ambient map induces torsion ordinate')
    check(s.expand(s.discriminant(cubic, X)-kappa**2*(kappa+1)**2) == 0,
          'exact singular-fibre boundary')

    low, high = -Q(16, 17), Q(482, 295)
    actions = [('H', Q(729, 1024), Q(241, 512)), ('GR', Q(8, 9), Q(0))]
    for r in range(4, 10):
        alpha = Q(2**(r+2), 3**r)
        actions.append((f'R{r}', alpha, -alpha/4))
    for j in range(1, 4):
        K = 2**(18*j)
        actions.append((f'J{j}', Q(8, 3*K), Q(16*(K-1), 81*K)))
    for name, alpha, beta in actions:
        check(0 < alpha <= Q(8, 9), 'contracting bank')
        check(low <= beta/(1-alpha) <= high, 'common fixed-point interval')

    ell, upper = Q(5346, 5605), Q(121, 119)
    compositions = 0
    # Algebraic envelope controls: these arbitrary words are NOT all native.
    for length in range(1, 4):
        for word in itertools.product(actions, repeat=length):
            A, B = Q(1), Q(0)
            for _, alpha, beta in word:
                A, B = alpha*A, alpha*B+beta
            gamma = B/(1-A)
            check(low <= gamma <= high, 'composition anchor is convex average')
            for Y in [Q(56), Q(112), Q(10000)]:
                X = (Y-B)/A
                for z in [Q(0), Q(1, 2), Q(1)]:
                    check(ell <= A*(X+z)/(Y+z) <= upper,
                          'depth-independent integrand distortion')
            compositions += 1

    pairs = []
    for k in range(25):
        n = 155+2048*k
        h = (729*n+669)//1024
        check((h+1)/Q(n+1) != 0, 'positive H endpoint')
        check(Q(h+1, 2) == Q(729, 1024)*Q(n+1, 2)+Q(241, 512), 'H offset')
        pairs.append((n, h, Q(729, 1024)))
        n, h = 251+288*k, 223+256*k
        check(replay(n, (1, 2, 1)) == replay(h, (1, 1, 1, 1, 3)), 'GR receipt')
        pairs.append((n, h, Q(8, 9)))
        for r in range(4, 8):
            odd_part = (7 if r % 2 == 0 else 5)+8*k
            n = (3**r*odd_part-1)//2
            h = 2**(r+1)*odd_part-1
            check(n % 4 == 3, 'native R guard')
            check(replay(n, (1,)) == replay(h, (1,)*r+(3,)), 'R receipt')
            pairs.append((n, h, Q(2**(r+2), 3**r)))
        for j in range(1, 4):
            K = 2**(18*j)
            n, h = (2264*K-23)//27+96*K*k, 223+256*k
            check(replay(n, (1, 2, 18*j, 1)) == replay(h, (1, 1, 1, 1, 3)),
                  'J native common-future receipt')
            pairs.append((n, h, Q(8, 3*K)))
    # The inherited genuine two-run family (no arbitrary guard composition).
    for k in range(25):
        n, mid, h = 13891+26244*k, 10975+20736*k, 8671+16384*k
        for a, b in [(n, mid), (mid, h)]:
            check(a % 4 == 3 and (2*a+1) % 81 == 0, 'two-run native guard')
            check(replay(a, (1,)) == replay(b, (1,)*4+(3,)), 'two-run receipts')
        pairs.append((n, h, Q(64, 81)**2))
    for n, h, A in pairs:
        check(h >= 111 and n > h and n % 2 == h % 2 == 1, 'actual index domain')
        for r in [1, 2, 3, 6]:
            ratio = nu(r, (n-1)//2)/(A**(r+1)*nu(r, (h-1)//2))
            check(upper**(-r-1) <= ratio <= ell**(-r-1), 'exact atomic weight bound')
    check(Q(9, 8) > 1, 'standalone G is outside contracting bank')
    check(replay(251, (1, 2)) == 283, 'neutral-cycle hostile forward edge')
    check((8*283-5)//9 == 251, 'neutral-cycle hostile inverse edge')
    for depth in range(1, 17):
        n, h = 2**(depth+1)-1, 2*3**depth-1
        check(replay(n, (1,)*depth) == h, 'actual unbounded climb family')
        for r in [1, 2, 3, 6]:
            ratio = nu(r, (n-1)//2)/nu(r, (h-1)//2)
            check(ratio >= Q(3**depth, 2**depth+1)**(r+1),
                  'new source prior cannot be a boundedly corrected discounted flow')

    # Historical hostile controls: score fibres and finite moment recovery.
    fibres = {}
    for mask in range(8):
        scores = [0, 0, 0]
        for bit, (i, j) in enumerate(itertools.combinations(range(3), 2)):
            scores[i if (mask >> bit) & 1 else j] += 1
        key = tuple(sorted(scores))
        fibres[key] = fibres.get(key, 0)+1
    check(sorted(fibres.values()) == [2, 6], 'sorted score fibre need not be a code coset')
    regular = 0
    for mask in range(1024):
        scores = [0]*5
        for bit, (i, j) in enumerate(itertools.combinations(range(5), 2)):
            scores[i if (mask >> bit) & 1 else j] += 1
        regular += scores == [2]*5
    check(regular == 24, 'labelled score fibre need not be a code coset either')
    triangle = s.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    swap = s.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    check(swap*triangle*swap.T == triangle.T and triangle != triangle.T,
          'self-converse does not mean self-adjoint')
    for k in range(1, 9):
        even = [Q(comb(k+1, j), 2**k) if j % 2 == 0 else Q(0) for j in range(k+2)]
        odd = [Q(comb(k+1, j), 2**k) if j % 2 else Q(0) for j in range(k+2)]
        check(sum(even) == sum(odd) == 1 and even[0] != odd[0], 'abstract atom hostile')
        for law in [even, odd]:
            moments = [sum(Q(comb(j, a))*law[j] for j in range(a, k+2)) for a in range(k+2)]
            check(sum((-1)**a*m for a, m in enumerate(moments)) == law[0],
                  'finite full moment recovery')
        for a in range(k+1):
            check(sum(comb(j, a)*(even[j]-odd[j]) for j in range(a, k+2)) == 0,
                  'proper truncated observers can forget the atom')

    print(json.dumps({'status': 'FINITE-EXACT; universal proofs in companion notes',
                      'checks': checks, 'polynomial_basis_size': len(basis),
                      'highest_t_degree_cap': 9, 'matrix_cases': rows,
                      'fibre_affine_stabilizers': 4,
                      'elliptic_singular_parameters': ['0', '-1'],
                      'algebraic_words_not_claimed_native': compositions,
                      'authenticated_or_inherited_pairs': len(pairs),
                      'anchor_interval': [str(low), str(high)],
                      'integrand_constants': [str(ell), str(upper)],
                      'r6_conditional_distortion_exact': str((upper/ell)**7),
                      'r6_conditional_distortion_decimal': float((upper/ell)**7),
                      'sorted_score_fibres_n3': {str(k): z for k, z in fibres.items()}},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
