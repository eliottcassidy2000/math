"""Exact repetition fuel and first-hit inverse-word controls.

Reproduce: python -X utf8 -B this_file.py (also with -O).
No convergence search or conjectural stopping criterion is used.
"""
from fractions import Fraction
from itertools import product
from math import comb


def need(condition, label):
    if not condition:
        raise RuntimeError(label)


def valuation(n, p):
    if n == 0:
        return None
    n = abs(n)
    result = 0
    while n % p == 0:
        n //= p
        result += 1
    return result


def inverse_affine(word):
    """Execution order; D(x)=2x and R(x)=(2x-1)/3."""
    need(set(word) <= set('DR'), 'inverse alphabet')
    q, d, b = 1, 1, 0
    for letter in word:
        q *= 2
        b *= 2
        if letter == 'R':
            b += d
            d *= 3
    return q, d, b


def inverse_literal(x, word):
    states = [Fraction(x)]
    for letter in word:
        x = 2 * states[-1]
        if letter == 'R':
            x = (x - 1) / 3
        states.append(x)
    return states


def predicted_integral_repeats(x, word):
    q, d, b = inverse_affine(word)
    r = word.count('R')
    need(r > 0, 'a ternary-consuming word is required')
    delta = (q - d) * x - b
    fuel = valuation(delta, 3)
    return None if fuel is None else fuel // r


def first_hit_inverse(word):
    states = inverse_literal(1, word)
    return all(x.denominator == 1 and x > 1 for x in states[1:])


def forward_literal(x, word):
    """0 is x/2, 1 is (3x+1)/2; keep rational failed controls."""
    states = [Fraction(x)]
    for letter in word:
        x = states[-1]
        states.append(x / 2 if letter == '0' else (3 * x + 1) / 2)
    return states


def forward_affine(word):
    q, d, c = 1, 1, 0
    for letter in word:
        if letter == '1':
            q *= 3
            c = 3 * c + d
        d *= 2
    return q, d, c


def main():
    inverse_checks = forward_checks = 0
    for length in range(1, 8):
        for letters in product('DR', repeat=length):
            word = ''.join(letters)
            q, d, b = inverse_affine(word)
            for x in range(1, 65):
                states = inverse_literal(x, word)
                need(states[-1] == Fraction(q*x-b, d), 'inverse affine replay')
                need((states[-1].denominator == 1) ==
                     all(z.denominator == 1 for z in states), 'integer endpoint restores prefixes')
                if not word.count('R'):
                    continue
                limit = predicted_integral_repeats(x, word)
                for k in range(0, 6):
                    repeated = inverse_literal(x, word * k)
                    actual = all(z.denominator == 1 for z in repeated)
                    need(actual == (limit is None or k <= limit), 'inverse exact fuel')
                    delta = (q-d)*x-b
                    need((q-d)*repeated[-1]-b == Fraction(q,d)**k*delta, 'inverse fixed-point identity')
                    inverse_checks += 1
        for letters in product('01', repeat=length):
            word = ''.join(letters)
            q, d, c = forward_affine(word)
            for x in range(1, 65):
                delta = (d-q)*x-c
                fuel = valuation(delta, 2)
                for k in range(0, 6):
                    states = forward_literal(x, word*k)
                    integral = all(z.denominator == 1 for z in states)
                    need(integral == (fuel is None or k*length <= fuel), 'forward exact fuel')
                    need((d-q)*states[-1]-c == Fraction(q,d)**k*delta, 'forward anchor identity')
                    forward_checks += 1
    print('Exact inverse repetition audits:', inverse_checks)
    print('Exact forward repetition audits:', forward_checks)

    # A complete finite comparison of strict root paths, without assuming Collatz.
    valid = []
    for length in range(0, 17):
        for letters in product('DR', repeat=length):
            word = ''.join(letters)
            if first_hit_inverse(word):
                n = int(inverse_literal(1, word)[-1])
                r, d = word.count('R'), word.count('D')
                need(d <= n.bit_length()-1+r, 'height bounds number of doublings')
                forward = ''.join('0' if c == 'D' else '1' for c in reversed(word))
                route = forward_literal(n, forward)
                need(route[-1] == 1 and all(z > 1 for z in route[:-1]), 'independent first-hit forward replay')
                valid.append((word,n))
    need(len({n for _,n in valid}) == len(valid), 'first-hit word uniqueness')
    print('All binary inverse words through length16:', 2**17-1,
          '; strict first-hit words:', len(valid))

    # Any repeated automaton state crossing an R would allow pumping this block.
    pumped = 0
    for word, _ in valid:
        states = inverse_literal(1, word)
        for i in range(len(word)):
            for j in range(i+1, len(word)+1):
                block = word[i:j]
                if 'R' not in block:
                    continue
                limit = predicted_integral_repeats(int(states[i]), block)
                need(limit is not None, 'no nontrivial cycle inside a strict root path')
                need(not first_hit_inverse(word[:i]+block*(limit+1)+word[j:]),
                     'pumped accepting path becomes invalid')
                pumped += 1
    print('All R-containing subwords of these strict paths:', pumped, '; hostile pump rejected')

    print('Explicit unbounded-R witnesses D^(3^(r-1)) R^r:')
    for r in range(2, 10):
        a = 3**(r-1)
        x = 2**a
        need(valuation(x+1,3) == r, 'LTE fuel witness')
        word = 'D'*a+'R'*r
        need(first_hit_inverse(word), 'explicit rooted witness')
        n = 2**r*(x+1)//3**r-1
        need(inverse_literal(1,word)[-1] == n, 'explicit endpoint')
        need(predicted_integral_repeats(x,'R') == r, 'exact failure after last allowed R')
        print(' r=',r,' initial_doublings=',a,' source_bits=',n.bit_length(),sep='')

    # A sound regular positive control has at most one R, but infinitely many words.
    # D^(2k+1) R D^b, k>=1 and b>=0, is integral and cannot revisit the root.
    controls = 0
    for k in range(1, 32):
        for b in range(33):
            word = 'D'*(2*k+1)+'R'+'D'*b
            need(first_hit_inverse(word), 'sound regular bounded-R control')
            controls += 1
    need(not first_hit_inverse('DR'), 'exclude the root cycle')
    need(all(x.denominator == 1 for x in inverse_literal(1,'DR'*12)), 'root-cycle integrality control')
    need(predicted_integral_repeats(1,'DR') is None, 'zero fuel numerator exceptional cycle')
    print('Sound regular one-R controls:', controls, '; root-cycle exception retained')

    # Two-anchor refills: high closeness to one center cannot be free high closeness
    # to a distinct center. Equality is deliberately left as the refill boundary.
    anchors = [Fraction(-1),Fraction(1),Fraction(-19,11),Fraction(-17),Fraction(3,5)]
    switches = collisions = 0
    def vp_fraction(z,p):
        if z == 0:
            return None
        return valuation(z.numerator,p)-valuation(z.denominator,p)
    for alpha in anchors:
        for beta in anchors:
            if alpha == beta:
                continue
            kappa = vp_fraction(alpha-beta,2)
            for n in range(1,2002,2):
                ka,kb = vp_fraction(n-alpha,2),vp_fraction(n-beta,2)
                if ka is None:
                    need(kb == kappa,'exact-anchor switch')
                elif ka != kappa:
                    need(kb == min(ka,kappa),'ultrametric switch')
                else:
                    need(kb is None or kb >= kappa,'collision boundary')
                    collisions += 1
                switches += 1
    print('Two-anchor precision switches:',switches,'; equality-boundary cases:',collisions)

    # Explicit upper bound for any language with at most K R letters.
    bounds = []
    X = 10**6
    for K in range(4):
        H = X.bit_length()-1+K
        bound = sum(comb(H+r+1,r+1) for r in range(K+1))
        seen = sum(n <= X and word.count('R') <= K for word,n in valid)
        need(seen <= bound,'polylog cardinality bound')
        bounds.append((K,bound))
    print('Polylog word-count bounds at X=1000000 for K=0..3:',bounds)
    print('PASS: arithmetic and finite-language controls; universal coverage remains OPEN')


if __name__ == '__main__':
    main()
