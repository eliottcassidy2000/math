#!/usr/bin/env python3
"""scaled_cores_and_clocks_meet_the_modular_split_boxeph_S220.py -- boxeph-2026-07-21-S220

Merge codex's SCALED ZETA-CORES (THM-2057) and the CLOCK structure (Z/q) into the modular-form arithmetic
entropy -- with the CORRECT modular attachment (per the S219 mining correction).

CORRECTION carried from S219: the Eisenstein (+) cusp split is NOT L(S) (THM-515 is a singular integral with
no Euler product, unsplit). It is the WEIGHT-2 second moment on X0(2p) (HYP-3587): floor = (#cusps-1)
Eisenstein (bulk 3/pi^2) (+) genus cusp forms (the obstruction). And the cusps of X0(N) ARE the divisors of
N = the sub-CLOCKS.

Joins:
  P1 SCALING: M is dilation-invariant; a scaled zeta-core carries the same LRC data (the Gamma0 level).
  P2 CLOCKS = CUSPS: divisors of N = cusps of X0(N) = sub-clocks. 12=2^2*3 (X0(12), 6 cusps, genus 0);
     14=2*7 (X0(14), cusps {1,2,7,14}, genus 1). Primes {2,3,7} = S215 Eisenstein(2,3)+apex(7).
  P3 THE SPLIT: floor 2nd moment on X0(2p) = (#cusps-1) EISENSTEIN (+) genus CUSP. Family genus 0,0,1,2,2
     (p=3,5,7,11,13): the FIRST cusp form is at p=7 (X0(14), f14=14a) = the first hard case. 12-clock genus 0
     (cuspless argmax floor) vs 14-clock genus 1 (the f14 obstruction). genus = the S218 deep/hidden entropy.
  P4 f14=14a: a_2=-1 (2-cusp), a_7=+1 (7-cusp) = the 2*7 clocks; period field Q(sqrt-7) = the S215 disc -7
     (apex, as f14's PERIOD field -- NOT a weight-1 theta). sym^2 f14 = the GL(3) 2nd-moment obstruction.
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
sep("P1  SCALING: M dilation-invariant; a scaled zeta-core carries the same LRC data")
deep=list(range(1,13))+[182]
Md=M_of(deep); print(f"  deep well {{1..12,182}}: M={Md:.5f} (~14/183={14/183:.5f})")
for c in (2,3,5):
    print(f"  c={c}: M(c*deep)={M_of([c*v for v in deep]):.5f} == M? {abs(M_of([c*v for v in deep])-Md)<1e-4}")
print("  => L(cS)=L(S): the scaled zeta-core (codex THM-2057) is the SAME LRC object on a refined clock.")

# ==========================================================================
sep("P2  CLOCKS = CUSPS of X0(N): the divisors of N are the cusps = the sub-clocks")
for N in (12,14):
    print(f"  {N}-clock = X0({N}): #cusps={num_cusps(N)} = divisors {divisors(N)} = sub-clocks ; primes {prime_factors(N)}")
print(f"  12=2^2*3 (Eisenstein/argmax primes 2,3; Phi_6 argmax t*=14/183) ; 14=2*7 (apex Paley-7, Phi_7)")
print(f"  gcd(12,14)={gcd(12,14)} (chirality prime) ; lcm={12*14//gcd(12,14)} (double-kill). X0(14) cusps {{1,2,7,14}}: apex-7 cusp hardest.")

# ==========================================================================
sep("P3  THE SPLIT (weight-2 X0(2p) 2nd moment): (#cusps-1) EISENSTEIN (+) genus CUSP; first cusp at p=7")
genus={6:0,10:0,14:1,22:2,26:2}   # X0(2p), p=3,5,7,11,13 (standard)
print("  p : n=2p : X0(n) #cusps : dim Eisenstein(=#cusps-1, the FLOOR/bulk) : genus(=#cusp forms, OBSTRUCTION)")
for p in (3,5,7,11,13):
    N=2*p; c=num_cusps(N); g=genus[N]
    tag=" <- FIRST cusp form: X0(14) genus 1 = f14=14a = the first HARD case (apex 7)" if p==7 else ""
    print(f"  {p:2d} :  {N:2d}  :   {c}       :        {c-1}                          :   {g}{tag}")
print("  12-clock X0(12): genus 0 => CUSPLESS (pure Eisenstein, the argmax/floor, no obstruction).")
print("  14-clock X0(14): genus 1 => ONE cusp form f14 (the LRC(14) obstruction). genus = the S218 DEEP/hidden")
print("  entropy: genus 0 (argmax) = zero deep entropy = rigid ; genus 1 (LRC) = the hidden cuspidal obstruction.")

# ==========================================================================
sep("P4  f14 = 14a: a_2=-1 (2-cusp), a_7=+1 (7-cusp) = the 2*7 clocks; period field Q(sqrt-7) = S215 disc -7")
# 14a q-expansion (standard): q - q^2 - 2q^3 + q^4 + 2q^6 + q^7 - q^8 + ... ; a_2=-1, a_7=+1
a={1:1,2:-1,3:-2,4:1,5:0,6:2,7:1,8:-1,9:1}
print(f"  f14 = 14a: a_2={a[2]} (mirror/2-cusp), a_7={a[7]} (apex/7-cusp) -- encodes the 2*7 = the two clocks;")
print(f"  Atkin-Lehner w_2=+1, w_7=-1, root number +1, RANK 0 => L(14a,1)>0 (favorable sign).")
print(f"  period field of f14 = Q(sqrt-7) = the S215 apex disc -7 (Paley-7) -- as f14's PERIOD field, NOT a weight-1 theta.")
print(f"  the GL(3) 2nd-moment obstruction = sym^2 f14: L(f14 x f14,s) = zeta(s) * L(sym^2 f14, s) = EISENSTEIN bulk (zeta, 3/pi^2) * CUSP obstruction (sym^2).")

# ==========================================================================
sep("P5  HONEST caveats (S219 correction + the mined negatives)")
print("""  * The Eisenstein(+)cusp split is the WEIGHT-2 X0(2p) second moment (HYP-3587), NOT L(S): THM-515's L(S)
    is a singular INTEGRAL with NO Euler product and is not split. (Corrects S219's attachment.)
  * The disc-(-7) binary form (S215/S219) is REAL but enters here as f14's PERIOD FIELD Q(sqrt-7); there is
    NO weight-1 dihedral/CM theta for LRC in the repo -- 'weight-1' elsewhere is combinatorial (false friend).
  * VALUE is modular, EXISTENCE is combinatorial: modular forms give the STRUCTURE+SIGN (genus=hardness,
    rank 0=favorable sign, Q(sqrt-7)=period, Z/6 torsion=units), but the floor CONSTANT is NOT L(14a,1), NOT
    L(sym^2), NOT any period -- it lives in the descent (covering-min product). (opus 3 negatives.)
  * MISTAKE-087: the Phi_6 covering-min construction is NON-extremal ('covering-min lives in class-number-1
    land' retracted to 'a particular non-minimal covering has class-number-1 structure').""")

sep("SUMMARY -- scaled cores x clocks x the CORRECT modular split")
print("""  codex's SCALED ZETA-CORE (dilation-invariant, THM-2057) reduced on the 12a/14a CLOCKS = the WEIGHT-2
  modular floor on X0(2p): the CUSPS of X0(N) ARE the divisor sub-clocks; the floor = (#cusps-1) EISENSTEIN
  (bulk 3/pi^2) (+) genus CUSP forms. The 12-clock (X0(12), 2^2*3, argmax Phi_6) is genus 0 = CUSPLESS =
  zero deep entropy (rigid); the 14-clock (X0(14), 2*7, apex Paley-7) is genus 1 = the ONE cusp form f14=14a
  = the LRC(14) OBSTRUCTION at the apex-7 cusp. f14 spells 2*7 (a_2=-1,a_7=+1), period field Q(sqrt-7) (the
  S215 disc -7). So genus (cusp forms) = the S218 deep/hidden arithmetic entropy; the difficulty of LRC(14)
  is the first cuspidal genus, at the apex prime 7 -- VALUE modular, EXISTENCE/constant combinatorial.""")
