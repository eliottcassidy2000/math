#!/usr/bin/env python3
"""MRB contour and branch-jump audit. Requires mpmath (no network).

The C enclosure uses interval arithmetic plus a proved Euler/Cauchy remainder.
Quadrature, theta, and root results are cross-checks, not interval certificates.
Run: python session31-mrb/certify_mrb.py > session31-mrb/certify_mrb.out
"""
from mpmath import mp, iv
import json

mp.dps = 100
iv.dps = 120
M, N = 256, 100

def fi(n):
    n = iv.mpf(n)
    return iv.exp(iv.log(n) / n) - 1

def ai(n):
    n = iv.mpf(n)
    return 2/iv.pi * iv.exp(-iv.log(n)/n) * iv.sin(iv.pi/n)

def f(z):
    return mp.expm1(mp.log(z)/z)

def a(x):
    x=mp.mpf(x)
    return 2/mp.pi * mp.exp(-mp.log(x)/x) * mp.sin(mp.pi/x)

# Identity S_M = sum_{j<N} (-1)^j Delta^j f(M)/2^(j+1) + remainder.
row = [fi(M+j) for j in range(N)]
S = iv.mpf(0)
for j in range(N):
    S += (-1)**j * row[0] / 2**(j+1)
    row = [row[k+1]-row[k] for k in range(len(row)-1)]
P = sum(((-1)**n * fi(n) for n in range(1,M)),iv.mpf(0))
L = iv.log(iv.mpf(3)*M/2) + iv.pi/6
b = 2*L/M
B = 2*iv.factorial(N)*iv.exp(b)*(L/iv.mpf(M)**(N+1) + (L/N+iv.mpf(1)/N**2)/iv.mpf(M)**N)
C = P + (-1)**M*S + iv.mpf([-1,1])*B
A8=ai(8)
D=A8-C
q=iv.exp(-iv.pi*iv.sqrt(2))
R=D/q**4
n=105719004
neighbors = {str(j): ai(j) for j in [n-1,n,n+1]}
# Compare D with both consecutive midpoints. Strict positivity certifies n.
mid_low = D-(ai(n)+ai(n+1))/2
mid_high = (ai(n-1)+ai(n))/2-D
remaining = D-ai(n)

print('INTERVAL CERTIFICATES (interval arithmetic + analytic tail bound)')
for name,value in [('Euler remainder bound',B),('C',C),('A8',A8),('D=A8-C',D),('q^4',q**4),('D/q^4',R),('D minus lower midpoint',mid_low),('upper midpoint minus D',mid_high),('D-A(105719004)',remaining),('A(105719003)-D',ai(n-1)-D)]:
    print(name,'=',value)
assert bool(B.b < iv.mpf('2.654e-84'))
assert bool(D.a > iv.mpf('1.891807127880073527054533938441983110587504636897727612419589120e-8'))
assert bool(D.b < iv.mpf('1.891807127880073527054533938441983110587504636897727612419589122e-8'))
assert bool(mid_low.a > 0) and bool(mid_high.a > 0)
assert bool(R.a > 0) and bool(R.b < 1)
assert bool((ai(n-1)-D).a > 0) and bool((D-ai(n)).a > 0)
assert bool((A8-C).a > 0) and bool((C-ai(9)).a > 0)
assert bool((C-(ai(8)+ai(9))/2).a > 0)
assert bool(((ai(7)+ai(8))/2-C).a > 0)

print('\nINDEPENDENT NUMERICAL CROSS-CHECKS (not interval certificates)')
C_nsum=mp.nsum(lambda n:(-1)**n*f(n),[1,mp.inf])
g=lambda t:mp.im(f(1+1j*t))
C_ap=mp.quad(lambda t:g(t)/mp.sinh(mp.pi*t),[0,1,2,4,8,16,mp.inf])
k=mp.sqrt(2)-1
qm=mp.exp(-mp.pi*mp.sqrt(2))
mu=mp.sqrt(2+mp.sqrt(2))
T=4*mp.sqrt(2)
cut_tail=mp.quad(lambda t:g(t)/mp.sinh(mp.pi*t),[T,T+2,T+8,mp.inf])
results={
 'C_nsum':C_nsum,'C_AbelPlana':C_ap,'difference_AP_nsum':C_ap-C_nsum,
 'k_r2':k,'Kprime_over_K':mp.ellipk(1-k*k)/mp.ellipk(k*k),
 'modular_lambda':(mp.jtheta(2,0,qm)/mp.jtheta(3,0,qm))**4,
 'lambda_minus_k_squared':(mp.jtheta(2,0,qm)/mp.jtheta(3,0,qm))**4-k*k,
 'honeycomb_mu':mu,'A8_minus_honeycomb_formula':a(8)-mp.root(2,8)/(mp.pi*mu),
 'A8_minus_CM_formula':a(8)-mp.power(2,-mp.mpf(1)/8)*mp.sqrt(k)/mp.pi,
 'cutoff_T':T,'vertical_tail_at_T':cut_tail,'vertical_tail_over_q4':cut_tail/qm**4,
 'm_star_with_A_m_equal_C':mp.findroot(lambda x:a(x)-C_nsum,8),
 'A_prime_8':mp.diff(a,8),
}
for name,value in results.items(): print(name,'=',mp.nstr(value,90))
print('\nODD LAPLACE MODES F_l = 2 int exp(-pi*l*t) Im f(1+it) dt')
for ell in [1,3,5,7,9,21,201]:
    omega=mp.pi*ell
    mode=2/omega*mp.quad(lambda u:mp.exp(-u)*g(u/omega),[0,1,4,10,50,200,mp.inf])
    print('ell',ell,'mode',mp.nstr(mode,60),'mode*pi^2*ell^2/2',mp.nstr(mode*omega**2/2,50))
print('\nBRANCH-JUMP-PRESERVING HOSTILE CONTROLS')
for K in [1,2,3,5,10]:
    h=lambda z:mp.exp(-z)*mp.fprod(z-j for j in range(1,K+1))
    direct=mp.fsum((-1)**n*h(n) for n in range(K+1,1001))
    exact=(-1)**(K+1)*mp.factorial(K)/(mp.e+1)**(K+1)
    ap=mp.quad(lambda t:mp.im(h(1+1j*t))/mp.sinh(mp.pi*t),[0,1,2,4,8,16,mp.inf])
    assert abs(direct-exact)<mp.mpf('1e-95') and abs(ap-exact)<mp.mpf('1e-95')
    print('K',K,'shift',mp.nstr(exact,60),'direct error',mp.nstr(direct-exact,5),'AP error',mp.nstr(ap-exact,5))
print('\nALL ASSERTIONS PASSED')
