#!/usr/bin/env python3
"""Exact identities and hostile controls for the minimizing-sequence proof.

No finite bank is used to infer compactness or a universal inequality.
Normal and optimized execution retain every gate. The analytic proof is
05-knowledge/results/long_frontier_sep06_minimizers.md.
"""
import hashlib
import sympy as S

u, v = S.sqrt(2), S.sqrt(3)
h, z = 1/u, 1/v
A, B = 2-u, u-1
cstar = (13-8*u)/3
K = 4*v/(3*(1+u)*(1+v))
gamma = (2-u)*(3-v)/4
C0 = (1+u+v-u*v)/2
GATES = 0
TRACE = hashlib.sha256()


def gate(ok, label):
    global GATES
    if ok != True:
        raise RuntimeError(label + ': ' + str(ok))
    GATES += 1
    TRACE.update((label+'\n').encode())


def zero(expr, label):
    gate(S.simplify(expr) == 0, label)


def identities_and_boundaries():
    p3, p4, t = S.symbols('p3 p4 t')
    E = (1-p4)/2
    D = (5-8*p3+3*p4)/6
    d2 = 2-u*t
    g = B-p3+A*p4
    C = C0-gamma*t
    zero(g-3*E*(D/E-cstar)/4, 'exact g-to-quotient identity')
    zero(1-C-p3+C*p4-3*E*(D/E-cstar-K*d2)/4,
         'exact F-to-excess identity')
    zero((1-cstar)/(2-u)-(u-S.Rational(2,3)), 'one-atom boundary constant')
    epsilon, lam = S.symbols('epsilon lam', real=True)
    delta, w2 = lam*epsilon, (1-lam)*epsilon
    tt = u*S.sqrt(1-delta-w2/2)
    pp3 = (tt**3+3*tt*w2)/4
    pp4 = (tt**4+6*tt**2*w2+w2**2)/8
    zero(pp3.subs(epsilon,0)-h, 'two-atom cubic constant')
    zero(S.diff(pp3,epsilon).subs(epsilon,0)
         -(-3*h*lam/2+3*h*(1-lam)/4), 'uniform two-atom cubic linear part')
    zero(pp4.subs(epsilon,0)-S.Rational(1,2), 'two-atom quartic constant')
    zero(S.diff(pp4,epsilon).subs(epsilon,0)-(-lam+1-lam),
         'uniform two-atom quartic linear part')
    zero(S.diff(2-u*tt,epsilon).subs(epsilon,0)-(lam+(1-lam)/2),
         'uniform two-atom distance linear part')
    ldust, lsplit = (28*u-32)/3, (64-44*u)/3
    zero(S.Rational(16,3)*(3*h/2-A)-ldust, 'dust boundary response')
    zero(S.Rational(32,3)*(A-3*h/4)-lsplit, 'split boundary response')
    zero(ldust-lsplit-(24*u-32), 'ordered local boundary constants')
    gate(2 > S.Rational(4,3)**2, 'dust constant exceeds split constant')
    gate(125**2 > 88**2*2, 'split constant exceeds one half')
    gate(2 > S.Rational(7,6)**2, 'one-atom constant exceeds one half')
    gate(3 < 2**2 and 2 > 1, 'K3 below four ninths by factor bounds')
    gate(S.Rational(4,9) < S.Rational(1,2), 'boundary gap is strict')
    C3 = C0-2*gamma*z
    zero(C3-(3-v)/2, 'three-atom limiting coefficient')
    zero(1-2*C3*z-(2-v), 'third-atom secant derivative margin')
    gate(2**2 > 3, 'positive third-atom derivative margin')
    zero(1-C3-3*z**3+C3*3*z**4, 'three-atom objective zero')
    R3 = (((5-8*3*z**3+3*3*z**4)/6)/S.Rational(1,3)-cstar)/(2-2*u*z)
    zero(R3-K, 'three-atom ratio is K3')
    print('BOUNDARIES one=sqrt(2)-2/3; two_liminf=(64-44sqrt(2))/3; both>1/2>K3')
    print('LOCAL_TWO_ATOM R=[Ldust*delta+Lsplit*w2/2]/[delta+w2/2]+o(1)')


def actual_mixed_dust_controls():
    # Exactly the declared seven actual lists; tails stored by multiplicity.
    universe = (4,5,6,8,12,20,50)
    for n in universe:
        N = n**4
        qa, qb, qc = 3*N+9, 6*(n-1), n*n+(n-1)**2-N
        a = (-qb+S.sqrt(qb*qb-4*qa*qc))/(2*qa)
        pos, neg = S.Rational(n,N), -(n+3*a-1)/N
        zero(qa*a*a+qb*a+qc, 'mixed dust defining quadratic n='+str(n))
        zero(3*a+N*pos+N*neg-1, 'mixed dust p1 n='+str(n))
        zero(3*a*a+N*pos*pos+N*neg*neg-1, 'mixed dust p2 n='+str(n))
        # Endpoint signs prove 1/2<a<z on the positive monotone branch.
        gate(S.Rational(3,4)+(n*n+(n+S.Rational(1,2))**2)/N<1,
             'mixed dust lower endpoint n='+str(n))
        gate(n*n>0, 'mixed dust upper endpoint n='+str(n))
        gate(S.Rational(n+1,N)<S.Rational(1,2), 'dust smaller than macros n='+str(n))
        zero(N*pos-n, 'positive dust first moment n='+str(n))
        zero(-N*neg-(n+3*a-1), 'negative dust first moment n='+str(n))
        # Actual distance formula and count consequence, using all signs.
        Delta = 2-2*v*a
        tail2 = N*(pos**2+neg**2)
        zero(Delta-3*(a-z)**2-tail2, 'actual square distance decomposition n='+str(n))
        zero(v-1-v*Delta/2-(3*a-1), 'negative count target n='+str(n))
    print('ACTUAL_HOSTILE n='+','.join(map(str,universe))+'; 3 macros, n^4 positive and n^4 negative dust entries')
    print('SEPARATE_DUST positive_sum=n; negative_magnitude=n+3a-1; signed_tail=1-3a')
    print('UNIVERSE seven exact algebraic lists; no exclusions; repeated zero/distance boundary controls are formal only')


def main():
    identities_and_boundaries()
    actual_mixed_dust_controls()
    print('CLAIM analytical iff: R->K3 iff permutation_l2_distance_to_three_equal_positive_atoms->0')
    print('COUNT N_minus*Delta3>=max(sqrt(3)-1-sqrt(3)*Delta3/2,0)^2')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())


if __name__ == '__main__':
    main()
