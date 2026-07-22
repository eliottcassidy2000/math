#!/usr/bin/env python3
"""Historical S220 modular/LRC synthesis; refuted by MISTAKE-233.

The script reproduces several independent controls: a numerical dilation
sample, exact cusp counts, standard genus data for five modular curves, and
two Fourier coefficients of the level-14 newform.  It does not construct an
LRC-to-modular-form map.  Modular cusps are not additive clock subgroups,
genus is not S218 entropy, and the non-CM level-14 newform has coefficient
field Q rather than a CM period field Q(sqrt(-7)).
"""
from math import gcd

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)
def frac_norm(x): x%=1.0; return min(x,1-x)
def fS(S,t): return min(frac_norm(v*t) for v in S)
def M_of(S,n=600000):
    best=0.0
    for i in range(n):
        t=(i+0.5)/n; f=fS(S,t)
        if f>best: best=f
    return best
def phi(n):
    r=n;m=n;p=2
    while p*p<=m:
        if m%p==0:
            while m%p==0:m//=p
            r-=r//p
        p+=1
    if m>1:r-=r//m
    return r
def divisors(N): return [d for d in range(1,N+1) if N%d==0]
def num_cusps(N): return sum(phi(gcd(d,N//d)) for d in divisors(N))
def prime_factors(n):
    n=abs(n);s=set();d=2
    while d*d<=n:
        while n%d==0:s.add(d);n//=d
        d+=1
    if n>1:s.add(n)
    return sorted(s)

# ==========================================================================
print("STATUS: HISTORICAL / LRC MODULAR-BRIDGE CLAIMS REFUTED (MISTAKE-233).")
print("SURVIVES: independent dilation, cusp-count, genus, and level-14 newform controls.")
print("NO LRC obstruction or predicate-preserving modular map is verified here.")
sep("P1  LRC DILATION CONTROL (independent of modular level)")
deep=list(range(1,13))+[182]
Md=M_of(deep); print(f"  deep well {{1..12,182}}: M={Md:.5f} (~14/183={14/183:.5f})")
for c in (2,3,5):
    print(f"  c={c}: M(c*deep)={M_of([c*v for v in deep]):.5f} == M? {abs(M_of([c*v for v in deep])-Md)<1e-4}")
print("  => Exact theorem: M(cS)=M(S). This gives no Gamma_0 level identification.")

# ==========================================================================
sep("P2  MODULAR CUSP COUNTS (divisor-indexed here; NOT additive clock subgroups)")
for N in (12,14):
    print(f"  X0({N}): #cusps={num_cusps(N)}; divisor labels {divisors(N)}; prime factors {prime_factors(N)}")
print("  Factorizations: 12=2^2*3; 14=2*7.")
print(f"  gcd(12,14)={gcd(12,14)} ; lcm={12*14//gcd(12,14)} (the THM-2057 double-kill modulus).")

# ==========================================================================
sep("P3  CLASSICAL WEIGHT-2 DIMENSIONS ON THE SELECTED FAMILY X0(2p)")
genus={6:0,10:0,14:1,22:2,26:2}   # X0(2p), p=3,5,7,11,13 (standard)
print("  p : n=2p : X0(n) #cusps : dim Eisenstein(=#cusps-1) : dim S_2(=genus)")
for p in (3,5,7,11,13):
    N=2*p; c=num_cusps(N); g=genus[N]
    tag=" <- first positive genus in this selected X0(2p) list" if p==7 else ""
    print(f"  {p:2d} :  {N:2d}  :   {c}       :        {c-1}                          :   {g}{tag}")
print("  X0(12) has genus 0; X0(14) has genus 1 and dim S_2=1.")
print("  These dimensions do not identify an LRC floor or obstruction; X0(11) already has genus 1.")

# ==========================================================================
sep("P4  LEVEL-14 NEWFORM DATA (no cusp-label or discriminant -7 bridge)")
# 14a q-expansion (standard): q - q^2 - 2q^3 + q^4 + 2q^6 + q^7 - q^8 + ... ; a_2=-1, a_7=+1
a={1:1,2:-1,3:-2,4:1,5:0,6:2,7:1,8:-1,9:1}
print(f"  f14 = 14a: Fourier coefficients a_2={a[2]}, a_7={a[7]}.")
print(f"  Atkin-Lehner w_2=+1, w_7=-1, root number +1, RANK 0 => L(14a,1)>0 (favorable sign).")
print("  coefficient field Q; elliptic isogeny class 14a is non-CM. There is no Q(sqrt(-7)) CM period bridge.")
print("  Classical identity: L(f14 x f14,s) = zeta(s) * L(sym^2 f14,s); no LRC consequence is supplied.")

# ==========================================================================
sep("P5  CORRECTED SCOPE")
print("""  * HYP-3587 proposes a weight-2 X0(2p) second-moment attachment, but no LRC map is proved here.
    THM-515's L(S) is a singular integral with no demonstrated modular split or Euler product.
  * The disc-(-7) binary theta in S219 is independent of f14; the level-14 newform is non-CM.
  * Neither VALUE nor EXISTENCE has been mapped to a modular invariant here. Genus, rank, L-values, and
    periods do not supply the LRC floor constant or a counterexample obstruction.
  * MISTAKE-087: the Phi_6 covering-min construction is NON-extremal ('covering-min lives in class-number-1
    land' retracted to 'a particular non-minimal covering has class-number-1 structure').""")

sep("SUMMARY -- HISTORICAL BRIDGE REFUTED")
print("""  Exact facts coexist: LRC dilation invariance and THM-2057's 12a/14a clocks; modular cusp counts and
  genera; and the non-CM level-14 newform. No construction identifies these objects or preserves the LRC
  witness predicate. Therefore f14 and sym^2 f14 are not established LRC(14) obstructions, modular cusps
  are not clock subgroups, and genus is not an arithmetic entropy. See MISTAKE-233.""")
