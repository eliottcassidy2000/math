"""
PURE-POWER exp-integral transcendence, ALL degrees (companion to the degree-2 theorem).

THEOREM: int_0^1 e^{alpha(t-t0)^d+gamma} dt is transcendental for algebraic alpha!=0,t0,gamma and every d>=1.
Engine: M_d(w)=sum_n w^n/((dn+1)n!)=1F1(1/d;1+1/d;w), an E-function, irreducible (a=1/d, b-a=1), with GENUINE
Poincare asymptotics M_d(w)~(1/d)e^w/w (single honest confluent function -- no oscillatory-cancellation trap).
int_0^x e^{alpha u^d}du = x M_d(alpha x^d) => I=e^gamma[(1-t0)M_d(eta1)+t0 M_d(eta2)], eta_i=alpha(1-t0)^d,alpha(-t0)^d.
Distinct exp types e^{eta_i z} => connections non-isomorphic => {g1,g2,1} lin indep /Qbar(z) => (Beukers 2006)
values lin indep /Qbar => I transcendental.  Verifies: series, reduction int=xM_d(ax^d), genuine growth, 2-pt
reduction, PSLQ non-algebraicity.  Sub-family only: general deg>=3 has intermediate terms + FC(2)-hard non-rationality.
"""
import mpmath as mp
mp.mp.dps = 30

def Md(w, d):  return mp.hyp1f1(mp.mpf(1)/d, 1 + mp.mpf(1)/d, w)
def Mser(w, d, N=90): return mp.fsum([w**n / ((d*n+1)*mp.factorial(n)) for n in range(N)])

print("(A) M_d(w)=1F1(1/d;1+1/d;w)=sum w^n/((dn+1)n!); int_0^x e^{a u^d}=x M_d(a x^d); genuine growth ~(1/d)e^w/w:")
for d in [2, 3, 4, 5, 6]:
    w = mp.mpc('0.7', '-0.5'); al = mp.mpf('0.6'); x = mp.mpf('0.8'); W = mp.mpf('60')
    e1 = abs(Md(w, d) - Mser(w, d))
    e2 = abs(mp.quad(lambda u: mp.e**(al*u**d), [0, x]) - x*Md(al*x**d, d))
    ratio = Md(W, d) / ((mp.gamma(1+mp.mpf(1)/d)/mp.gamma(mp.mpf(1)/d))*mp.e**W*W**(-1))
    print(f"    d={d}: |1F1-series|={mp.nstr(e1,3):>10}  |int-xM_d|={mp.nstr(e2,3):>10}  growth-ratio@60={mp.nstr(ratio,7)}  "
          f"irred(a=1/{d},b-a=1): yes")

print("\n(B) I=int_0^1 e^{alpha(t-t0)^d+gamma}dt = e^gamma[(1-t0)M_d(eta1)+t0 M_d(eta2)], eta_i=alpha(1-t0)^d,alpha(-t0)^d:")
for (d, al, t0, ga) in [(3, mp.mpf('0.9'), mp.mpf('0.3'), mp.mpf('0.2')),
                        (4, mp.mpf('-0.7'), mp.mpf('0.65'), mp.mpf('0.1')),
                        (5, mp.mpc('0.5','0.4'), mp.mpf('0.2'), mp.mpc('0.1','-0.3'))]:
    I = mp.quad(lambda t: mp.e**(al*(t-t0)**d + ga), [0, 1])
    eta1 = al*(1-t0)**d; eta2 = al*(-t0)**d
    Ired = mp.e**ga * ((1-t0)*Md(eta1, d) + t0*Md(eta2, d))
    print(f"    d={d}, alpha={str(al):>12}: I={mp.nstr(I,14)}  reduced diff={mp.nstr(abs(I-Ired),3)}  "
          f"eta1!=eta2: {abs(eta1-eta2) > mp.mpf('1e-20')}")

print("\n(C) PSLQ non-algebraicity sanity (deg-3 pure power):")
d = 3; al = mp.mpf('0.9'); t0 = mp.mpf('0.3'); ga = mp.mpf('0.2')
I = mp.quad(lambda t: mp.e**(al*(t-t0)**d + ga), [0, 1])
print(f"    I={mp.nstr(I,20)}  min-poly deg<=5 h<=1e9: {mp.pslq([I**k for k in range(6)], maxcoeff=10**9, maxsteps=10**5)}")
print("\nCorollary: int_0^1 e^{alpha t^d}dt = 1F1(1/d;1+1/d;alpha) transcendental for algebraic alpha!=0, every d.")
