"""
DEGREE-2 exp-integral transcendence, rigorous proof -- numerical verification of every non-cited step.

THEOREM. For alpha,beta,gamma in Qbar (alpha!=0), I = int_0^1 e^{alpha t^2+beta t+gamma} dt is transcendental.

Proof skeleton (see 07-reflections/degree-2-exp-integral-is-transcendental-...-opus-S4.md):
  (*)  I = e^delta[(1+c)M(eta1) - c M(eta2)],  c=beta/2alpha, delta=gamma-beta^2/4alpha,
       M(z)=1F1(1/2;3/2;z)=sum z^n/((2n+1)n!) an E-FUNCTION, eta1=alpha(1+c)^2, eta2=alpha c^2.
  g1=e^{del z}M(eta1 z), g2=e^{del z}M(eta2 z), g0=1 are E-functions, lin. indep. over Qbar(z)
       because g_i ~ e^{(del+eta_i)z}/(2 eta_i z) have DISTINCT rates del+eta1!=del+eta2 (eta1!=eta2).
  Refined Siegel-Shidlovskii (Beukers 2006) => values at z=1 are Qbar-linearly independent
       => I=(1+c)g1(1)-c g2(1) is transcendental.  Degenerate c in {0,-1,-1/2} => single value e^lambda M(eta).

This script verifies: (A) the E-function series, (B) reduction (*), (C) distinct exponential rates (rigidity),
(D) degenerate collapses, (E) PSLQ non-algebraicity sanity.  The ONE non-elementary input (Beukers 2006) is cited.
"""
import mpmath as mp
mp.mp.dps = 40

def M(z):                      # 1F1(1/2;3/2;z)
    return mp.hyp1f1(mp.mpf(1)/2, mp.mpf(3)/2, z)

def Mseries(z, N=100):         # E-function form: sum z^n/((2n+1) n!), coefficient b_n = 1/(2n+1)
    return mp.fsum([z**n / ((2*n+1) * mp.factorial(n)) for n in range(N)])

def I_direct(al, be, ga):
    return mp.quad(lambda t: mp.e**(al*t*t + be*t + ga), [0, 1])

def I_reduced(al, be, ga):     # (*)
    c = be/(2*al); delta = ga - be*be/(4*al)
    eta1 = al*(1+c)**2; eta2 = al*c*c
    return mp.e**delta * ((1+c)*M(eta1) - c*M(eta2))

print("="*78)
print("(A) E-FUNCTION FORM  1F1(1/2;3/2;z) = sum z^n/((2n+1) n!)   [b_n = 1/(2n+1)]")
for z in [mp.mpc('0.6','-0.4'), mp.mpf('2.3'), mp.mpc('-1.1','0.7')]:
    print(f"    z={str(z):>18}:  |1F1 - series| = {mp.nstr(abs(M(z)-Mseries(z)),3)}")

print("="*78)
print("(B) REDUCTION (*)  I = e^delta[(1+c)M(eta1) - c M(eta2)]   (generic + complex coeffs)")
for (al,be,ga) in [(mp.mpf('0.7'),mp.mpf('0.9'),mp.mpf('-0.2')),
                   (mp.mpc('0.7','0.3'),mp.mpc('-0.4','0.8'),mp.mpc('0.1','0.2')),
                   (mp.mpf('-1.3'),mp.mpf('2.1'),mp.mpf('0.5'))]:
    d = abs(I_direct(al,be,ga) - I_reduced(al,be,ga))
    print(f"    alpha={str(al):>16}: |I_direct - I_reduced| = {mp.nstr(d,3)}")

print("="*78)
print("(C) RIGIDITY = DISTINCT EXPONENTIAL RATES.  g_i(z)=e^{del z}M(eta_i z) ~ e^{(del+eta_i)z}/(2 eta_i z).")
al=mp.mpf('0.7'); c=mp.mpf('0.35'); delta=mp.mpf('0.1')      # generic c
eta1=al*(1+c)**2; eta2=al*c*c
def g(eta,z): return mp.e**(delta*z)*M(eta*z)
for z in [mp.mpf('40'), mp.mpf('80')]:
    r1=mp.re(mp.log(g(eta1,z)))/z; r2=mp.re(mp.log(g(eta2,z)))/z
    print(f"    z={float(z):>4g}: rate1={mp.nstr(r1,10)} (->{mp.nstr(delta+eta1,10)})   rate2={mp.nstr(r2,10)} (->{mp.nstr(delta+eta2,10)})")
print(f"    del+eta1={mp.nstr(delta+eta1,12)} != del+eta2={mp.nstr(delta+eta2,12)}  =>  g1,g2,1 lin indep over Qbar(z).")

print("="*78)
print("(D) DEGENERATE PARAMETERS collapse to a single confluent value e^lambda M(eta):")
al=mp.mpf('1.3'); ga=mp.mpf('0.4')
# c=0 (beta=0): I=e^gamma M(alpha)
be=mp.mpf('0'); print(f"    c=0    : |I - e^gamma M(alpha)|   = {mp.nstr(abs(I_direct(al,be,ga)-mp.e**ga*M(al)),3)}")
# c=-1 (beta=-2 alpha): I=e^delta M(alpha),  delta=ga-be^2/4al
be=-2*al; delta=ga-be*be/(4*al); print(f"    c=-1   : |I - e^delta M(alpha)|   = {mp.nstr(abs(I_direct(al,be,ga)-mp.e**delta*M(al)),3)}")
# c=-1/2 (beta=-alpha): I=e^delta M(alpha/4)
be=-al; delta=ga-be*be/(4*al); print(f"    c=-1/2 : |I - e^delta M(alpha/4)| = {mp.nstr(abs(I_direct(al,be,ga)-mp.e**delta*M(al/4)),3)}")

print("="*78)
print("(E) PSLQ non-algebraicity sanity (cannot prove transcendence numerically; rules out low-height algebraic):")
al,be,ga=mp.mpf('0.7'),mp.mpf('0.9'),mp.mpf('-0.2')
I=I_direct(al,be,ga)
rel=mp.pslq([I**k for k in range(6)], maxcoeff=10**9, maxsteps=10**5)
print(f"    I={mp.nstr(I,20)};  min-poly deg<=5 height<=1e9:  {rel}")
print("    (None/huge-height => no low-degree low-height algebraic relation; consistent with the theorem.)")
print("="*78)
print("Cited (not numerically checkable): Beukers, Ann. of Math. 163 (2006) 369-379 -- refined linear S-S,")
print("no exceptional point => Qbar(z)-lin-indep E-functions have Qbar-lin-indep values at every nonzero alg xi.")
