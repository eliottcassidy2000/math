"""Independent exact referee: power-sum Newton identities and sparse algebra.

No producer import. The universal claims are audited in the adjacent report;
these controls reconstruct the finite certificates by different algorithms.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(newline='\n')
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(label)


def certificate(stem):
    here = Path(__file__).resolve().parent
    candidates = [here/(stem+'_certificate.json'),
                  here.parent/'05-knowledge/results'/(stem+'_certificate.json')]
    path = next((p for p in candidates if p.exists()), None)
    need(path is not None, 'adjacent/installed certificate')
    pins = {
        'continuing8_20260906_newton_clusters': '1d35470321fbc24289b62fb49955cafe9e10e7518f26b9ae17f59ee12145e17e',
        'continuing8_20260906_ballot_repair': '151ddbd2240c210ccaf87f5f40f101439af2d58bf4fcebd765a6d6c06cc379d2'}
    need(sha256(path.read_bytes()).hexdigest() == pins[stem], 'frozen complete certificate pin')
    return json.loads(path.read_bytes()), sha256(path.read_bytes()).hexdigest()


def newton_from_power_sums(roots):
    n = len(roots)
    s = [0] + [sum(x**k for x in roots) for k in range(1, n+1)]
    e = [F(1)]
    for k in range(1, n+1):
        e.append(sum((-1)**(j-1)*e[k-j]*s[j] for j in range(1, k+1))/k)
    need(e[-1] == prod(roots), 'independent Newton determinant coefficient')
    return e


def ratios(e):
    d = len(e)-1
    h = [x/comb(d, k) for k, x in enumerate(e)]
    return {k: h[k]**2/(h[k-1]*h[k+1]) for k in range(1, d)}


def circuits(R):
    return {k: R[k]/R[k-1] for k in range(2, len(R)+1)}


def signs(c):
    return [(x > 1)-(x < 1) for x in c.values()]


def changes(w):
    w = [x for x in w if x]
    return sum(x != y for x, y in zip(w, w[1:]))


def leading_semiring(profile):
    """Max-plus elementary-symmetric DP with multiplicities of tied optima."""
    out = [(0, 1)]
    K = len(profile)
    for j, count in enumerate(profile):
        for _ in range(count):
            old = out[:]
            out.append((-10**9, 0))
            for k in range(1, len(out)):
                take = (old[k-1][0]+K-1-j, old[k-1][1])
                skip = old[k] if k < len(old) else (-10**9, 0)
                out[k] = take if take[0] > skip[0] else skip if skip[0] > take[0] else (
                    take[0], take[1]+skip[1])
    return out


def factorial_kappa(profile):
    d, a, out = sum(profile), 0, {}
    for j, m in enumerate(profile):
        for ell in range(1, m):
            k = a+ell
            out[k] = F((ell+1)*(m-ell+1)*k*(d-k),
                       ell*(m-ell)*(k+1)*(d-k+1))
        a += m
        if a < d:
            out[a] = F(a*(d-a), m*profile[j+1]*(a+1)*(d-a+1))
    return out


def cluster_audit(data):
    actual_rows, literal_ties = 0, 0
    for record in data['profiles']:
        profile = record['profile']
        d, K = sum(profile), len(profile)
        bounds = [sum(profile[:j]) for j in range(1, K)]
        kap = factorial_kappa(profile)
        leading = leading_semiring(profile)
        leadkap = ratios([F(c) for power, c in leading])
        need(kap == leadkap, 'factorial formula versus max-plus coefficient DP')
        need([str(kap[k]) for k in range(1, d)] == record['kappa'], 'all pinned limiting ratios')
        for k in range(1, d):
            need(2*leading[k][0]-leading[k-1][0]-leading[k+1][0] == int(k in bounds),
                 'complete dominant boundary-exponent support')
        c0 = circuits(kap)
        margins = [c0[k] for k in range(2, bounds[0])]
        margins += [1/c0[k] for k in range(bounds[-1]+2, d)]
        for a, b in zip(bounds, bounds[1:]):
            margins += [c0[k+1]/c0[k] for k in range(a+2, b-1)]
        need(all(x > 1 for x in margins), 'all strict limiting interval margins')
        eta = min(margins, default=F(2))
        M = 1
        while F(M+1, M)**8 >= eta:
            M *= 2
        tau = F(M+1, M)
        spike = max([1/c0[b] for b in bounds]+[c0[b+1] for b in bounds])
        v = tau**4*spike
        ceilv = -((-v.numerator)//v.denominator)
        T = max(3*M*2**d, ceilv+1)
        eps = F(1, 6*M*d)
        need((M, str(tau), T, str(eps), bounds) ==
             (record['M'], record['tau'], record['separation'],
              record['relative_cluster_width'], record['boundaries']), 'complete threshold reconstruction')
        # Independent geometric-series proof budget, stronger than evaluating the product.
        need((1+eps)**d < 1/(1-d*eps), 'finite binomial/geometric majorant')
        need((1+F(1, 3*M))/(1-F(1, 6*M)) < tau, 'uniform analytic coefficient budget')
        for control in record['controls']:
            roots = list(map(F, control['roots']))
            bases = list(map(F, control['bases']))
            need(all(bases[j]/bases[j+1] >= T for j in range(K-1)), 'all actual separations')
            start = 0
            for base, m in zip(bases, profile):
                need(all(base <= x <= (1+eps)*base for x in roots[start:start+m]), 'all cluster boxes')
                start += m
            e = newton_from_power_sums(roots)
            R = ratios(e)
            c = circuits(R)
            need(signs(c) == control['sign_word'], 'entire actual circuit word')
            need(changes(signs(c)) == 2*K-3, 'exact zero-filtered sign count')
            literal_ties += signs(c).count(0)
            k, prefix = 0, F(1)
            for base, m in zip(bases, profile):
                for ell in range(1, m+1):
                    k += 1
                    E = e[k]/(prefix*base**ell*comb(m, ell))
                    need(1 <= E < tau, 'every actual coefficient in uniform envelope')
                prefix *= base**m
            need(all(c[k] > 1 for k in range(2, bounds[0]+1)), 'all first-band signs')
            need(all(c[k] < 1 for k in range(bounds[-1]+1, d)), 'all final-band signs')
            for b in bounds:
                need(c[b] > 1 and c[b+1] < 1, 'both signs of each spike')
            for a, b in zip(bounds, bounds[1:]):
                need(changes([signs(c)[k-2] for k in range(a+1, b+1)]) == 1,
                     'exactly one inter-band reversal including ties')
                need(all(c[k+1] > c[k] for k in range(a+2, b-1)), 'complete interior ordering')
            reverse = circuits(ratios(newton_from_power_sums([1/x for x in roots])))
            need(all(reverse[k] == 1/c[d+1-k] for k in c), 'literal reciprocal circuit map')
            actual_rows += 1
    for hostile in data['hostiles']:
        c = circuits(ratios(newton_from_power_sums(list(map(F, hostile['roots'])))))
        need(signs(c) == hostile['word'], 'all canonical scope hostiles')
    need(literal_ties > 0, 'exact zero boundary audited')
    return actual_rows, literal_ties


def poly_value(coeff, x):
    v = 0
    for c in coeff[::-1]:
        v = v*x+c
    return v


def convolve(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def factorial_transport_audit(data):
    item = data['factorial_transport']
    brackets = {2: [(-4,-3),(-1,0)],
                3: [(-7,-6),(-3,-2),(-1,0)],
                4: [(-10,-9),(-5,-4),(-2,-1),(-1,0)]}
    target = [F(1)]
    for row in item['factors']:
        m = row['degree']
        # Gregory-Newton coefficients of (x+1)_m, by finite differences.
        values = [prod(i+j for j in range(1,m+1)) for i in range(m+1)]
        coeff = []
        for j in range(m+1):
            coeff.append(F(values[0],factorial(j)))
            values = [b-a for a,b in zip(values, values[1:])]
        need(list(reversed(coeff)) == list(map(F,row['original_descending_coefficients'])),
             'finite-difference factorial row independent reconstruction')
        need(coeff[-1] == 1 and all(c > 0 for c in coeff), 'positive monic PF factor')
        for left, right in brackets[m]:
            need(poly_value(coeff,F(left))*poly_value(coeff,F(right)) < 0,
                 'different exact sign isolation of every negative root')
        need(len(brackets[m]) == m, 'degree exhausts negative simple roots')
        bound = 1+max(coeff[:-1])
        shift, dilation, base = map(F,(row['shift'],row['scale'],row['base']))
        eps = F(item['epsilon'])
        need(bound == F(row['Cauchy_bound']) and shift >= bound/eps,
             'Cauchy-to-relative-width typed bound')
        need(dilation*shift == base, 'actual root-parameter clock')
        # Derivatives at the translating point, instead of symbolic substitution.
        moved = [dilation**(m-k)*sum(coeff[j]*comb(j,k)*shift**(j-k)
                  for j in range(k,m+1)) for k in range(m+1)]
        need(moved[-1] == 1 and all(c > 0 for c in moved), 'literal moved factor')
        target = convolve(target,moved)
    e = list(reversed(target))
    need(e == list(map(F,item['target_descending_coefficients'])), 'all derivative-transport coefficients')
    c = circuits(ratios(e))
    need(signs(c) == item['sign_word'] == [1,-1,1,1,-1,-1,-1], 'literal new factorial consumer word')
    need(changes(signs(c)) == 3, 'factorial consumer has the predicted exact count')
    d = sum(item['profile'])
    eps = F(item['epsilon'])
    M = 1/(6*d*eps)
    need(M.denominator == 1, 'factorial consumer integer threshold parameter')
    tau = (M+1)/M
    kap = factorial_kappa(item['profile'])
    c0 = circuits(kap)
    bounds = [sum(item['profile'][:j]) for j in range(1,len(item['profile']))]
    k, prefix = 0, F(1)
    for row in item['factors']:
        base, m = F(row['base']), row['degree']
        for ell in range(1,m+1):
            k += 1
            need(1 <= e[k]/(prefix*base**ell*comb(m,ell)) < tau,
                 'factorial consumer full coefficient envelope')
        prefix *= base**m
    need(all(c[b] > 1 and c[b+1] < 1 for b in bounds), 'factorial consumer strict boundary pair signs')


# Sparse Z[k,a,b] arithmetic gives an independent symbolic ballot identity.
def const(n):
    return {(0, 0, 0): n} if n else {}


def pa(*items):
    z = {}
    for p in items:
        for e, c in p.items():
            z[e] = z.get(e, 0)+c
    return {e: c for e, c in z.items() if c}


def pm(*items):
    z = const(1)
    for p in items:
        out = {}
        for a, c in z.items():
            for b, v in p.items():
                e = tuple(x+y for x, y in zip(a, b))
                out[e] = out.get(e, 0)+c*v
        z = {e: c for e, c in out.items() if c}
    return z


def ps(p, n):
    return {e: n*c for e, c in p.items() if n*c}


def ballot_audit(data):
    k, a, b = {(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}
    t, u, v = pa(ps(k, 2), a), pa(k, b), pa(k, a, ps(b, -1))
    numerator = pm(t, pa(t, const(-1)), pa(u, const(1)), pa(v, const(1)))
    denominator = pm(u, v, pa(t, const(2)), pa(t, const(1)))
    delta2 = pm(pa(a, ps(b, -2)), pa(a, ps(b, -2)))
    q4 = pa(ps(pm(k, k), 4), ps(pm(pa(a, const(1), ps(delta2, -1)), k), 4),
            pm(a, pa(a, const(2))), ps(pm(pa(ps(a, 2), const(1)), delta2), -1))
    need(pa(ps(denominator, 2), ps(numerator, -2), ps(q4, -1)) == {},
         'universal sparse-polynomial all-shift identity')
    DD = pa(ps(a, 2), const(1))
    original = pa(pm(a, pa(a, const(2))), const(4), ps(pm(DD, delta2), -1))
    transformed = pa(pm(DD, DD), ps(DD, 2), const(13), ps(pm(DD, delta2), -4))
    need(pa(ps(original, 4), ps(transformed, -1)) == {}, 'universal divisor identity')
    p4den = pm(k,pa(ps(k,2),const(1)),pa(ps(k,3),const(-1)),
               pa(ps(k,3),const(1)),pa(ps(k,4),const(1)),pa(ps(k,4),const(3)))
    p4num = pm(pa(ps(k,4),const(-1)),pa(ps(k,2),const(-1)),
               pa(ps(k,4),const(-3)),pa(k,const(1)),
               pa(ps(k,3),const(2)),pa(ps(k,3),const(4)))
    quartic = {(4,0,0):432,(3,0,0):404,(2,0,0):-129,(1,0,0):-101,(0,0,0):24}
    need(pa(p4den,ps(p4num,-1),ps(quartic,-1)) == {}, 'universal p4 quartic numerator identity')
    for root in (F(0),F(-1,2),F(1,3),F(-1,3),F(-1,4),F(-3,4)):
        need(sum(c*root**e[0] for e,c in quartic.items()) != 0, 'quartic numerator has no denominator cancellation')
    alternatives = []
    for D in (1, -1, 13, -13):
        square = F(D+2, 4)+F(13, 4*D)
        need(square in (4, -3), 'complete signed divisors of prime thirteen')
        if square < 0:
            continue
        aa = (D-1)//2
        for delta in (-2, 2):
            bb = (aa-delta)//2
            linear = aa+1-delta**2
            need(F(aa*(aa+2)-(2*aa+1)*delta**2, 4) == -1, 'constant minus one')
            need(aa-2*bb == delta, 'integer shift reconstruction')
            alternatives.append((aa, bb, linear))
    need(sorted((a,b) for a,b,c in alternatives if c <= -1) == [(0,-1),(0,1)],
         'global positive-metallic classification')
    # An independent, different finite rectangle, not needed for the global proof.
    literal = 0
    for aa in range(-4, 9):
        for bb in range(-4, 7):
            for kk in range(max(1, 1-bb, 1-aa+bb), 14):
                h = [comb(2*j+aa, j+bb) for j in (kk-1, kk, kk+1)]
                den = (kk+bb)*(kk+aa-bb)*(2*kk+aa+2)*(2*kk+aa+1)
                delt = aa-2*bb
                Q = kk**2+(aa+1-delt**2)*kk+F(aa*(aa+2)-(2*aa+1)*delt**2,4)
                need(den > 0 and 1-F(h[1]**2,h[0]*h[2]) == 2*Q/den,
                     'literal binomials versus universal quadratic')
                literal += 1
    fuss_cases = 0
    for order in range(2, 11):
        for kk in range(2, 13):
            def ff(j):
                return F(comb(order*j,j),(order-1)*j+1)
            actual = ff(kk)/ff(kk-1)
            formula = F(prod(order*(kk-1)+j for j in range(1,order+1)),
                        kk*prod((order-1)*(kk-1)+j for j in range(1,order)))
            formula *= F((order-1)*(kk-1)+1,(order-1)*kk+1)
            need(actual == formula, 'different Fuss factorial bank')
            if order == 4:
                num = 432*kk**4+404*kk**3-129*kk**2-101*kk+24
                den = kk*(2*kk+1)*(3*kk-1)*(3*kk+1)*(4*kk+1)*(4*kk+3)
                need(1-ff(kk)**2/(ff(kk-1)*ff(kk+1)) == F(num,den), 'explicit quartic Fuss response')
            fuss_cases += 1
    hostile = data['reciprocal_hostile']
    R = ratios(newton_from_power_sums(list(map(F,hostile['roots']))))
    need(list(R.values()) == list(map(F, hostile['normalized_ratios'])), 'literal reciprocal normalized row')
    need(signs(circuits(R)) == [1,-1,1,-1], 'simultaneous antipalindrome and maximal alternation')
    return literal, fuss_cases


def main():
    clusters, csha = certificate('continuing8_20260906_newton_clusters')
    ballot, bsha = certificate('continuing8_20260906_ballot_repair')
    rows, ties = cluster_audit(clusters)
    factorial_transport_audit(clusters)
    literal, fuss = ballot_audit(ballot)
    print('Independent Newton/ballot audit: exact controls PASS')
    print('Cluster certificate SHA256', csha)
    print('Ballot certificate SHA256', bsha)
    print('Profiles',len(clusters['profiles']),'power-sum rows',rows,'exact ties',ties)
    print('Independent ballot rows',literal,'Fuss ratios',fuss)
    print('Universal separation/count and bronze proofs reviewed separately; finite controls are not quantifiers.')
    print('PASS',GATES,'always-active exact gates; raw LF')


if __name__ == '__main__':
    main()
