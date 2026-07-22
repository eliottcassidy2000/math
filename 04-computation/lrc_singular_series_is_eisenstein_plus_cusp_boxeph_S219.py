#!/usr/bin/env python3
"""lrc_singular_series_is_eisenstein_plus_cusp_boxeph_S219.py -- boxeph-2026-07-21-S219

STATUS: HISTORICAL / LRC INHERITANCE REFUTED BY MISTAKE-233.

The finite binary-theta coefficient identities in parts A and B survive.
THM-515 is not shown to be a modular theta series, so parts C and D and the
summary do not establish an Eisenstein/cusp split, entropy, or obstruction for
LRC.  They remain below only to preserve the historical failure mechanism.

Tie previous MODULAR-FORM work to cutting-edge LRC. The LRC good-set measure / singular series is a lattice
THETA (THM-515); theta functions of quadratic forms are MODULAR FORMS, which split canonically as

   theta  =  EISENSTEIN  (+)  CUSP.

And this split IS my S217/S218 arithmetic-entropy split:
   EISENSTEIN part  <->  the GENUS / singular-series MAIN TERM / the FLOOR   =  LOCAL (congruence-visible),
                         computable from Legendre data, ZERO deep entropy.
   CUSP part        <->  the DEEP class-group fluctuation                     =  HIDDEN arithmetic entropy.

Apex-prime Heegner h=1 (disc -7 = Paley-7 spectrum, LRC(14)=2*7) => the weight-1 theta is PURE EISENSTEIN
(no cusp forms) => zero deep entropy => RIGID. h>1 (disc -23) => cusp forms appear = the hidden fluctuation
where a counterexample must hide. Verified.
"""
from math import gcd

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def rep_count(form,n,B=80):
    a,b,c=form; cnt=0
    for x in range(-B,B+1):
        for y in range(-B,B+1):
            if a*x*x+b*x*y+c*y*y==n: cnt+=1
    return cnt
def divisor_char_sum(n,p):   # Eisenstein coeff: sum_{d|n} chi_{-p}(d), chi = Legendre(.,p) (0 if p|d)
    return sum(legendre(d,p) for d in range(1,n+1) if n%d==0)

# ==========================================================================
print("STATUS: HISTORICAL / LRC INHERITANCE REFUTED (MISTAKE-233).")
print("SURVIVES: finite binary-theta identities in A-B; C-D/SUMMARY have no proved LRC map.")
sep("A  disc -7 binary theta: PURE EISENSTEIN  r_Q(n)=2*sum_{d|n}(d/7)")
Q7=(1,1,2)  # principal form of disc -7
ok=True
for n in range(1,41):
    r=rep_count(Q7,n); eis=2*divisor_char_sum(n,7)
    if r!=eis: ok=False; print(f"   MISMATCH n={n}: r={r} eis={eis}")
print(f"  r_Q(n) == 2*sum_{{d|n}}(d/7) for n=1..40 (pure weight-1 EISENSTEIN, NO cusp)? {ok}")
print("  sample r_Q(n):", [rep_count(Q7,n) for n in range(1,16)])
print("  => This binary theta is Eisenstein-only. No LRC or entropy consequence follows.")

# ==========================================================================
sep("B  disc -23: class average is EISENSTEIN; class difference is a nonzero CUSP form")
Qprin=(1,1,6); Qnon=(2,1,3)   # principal + a non-principal class (its inverse (2,-1,3) reps the same n)
eis_ok=True; cusp_nonzero=0
print("   n :  r_prin  r_non  | (r_prin+2 r_non)/2  vs  sum_{d|n}(d/23)  |  cusp = r_prin - r_non")
for n in range(1,31):
    rp=rep_count(Qprin,n); rn=rep_count(Qnon,n)
    eis=divisor_char_sum(n,23); genus=(rp+2*rn)//2 if (rp+2*rn)%2==0 else None
    cusp=rp-rn
    if genus!=eis: eis_ok=False
    if cusp!=0: cusp_nonzero+=1
    if n<=16 or cusp!=0 and n<=24:
        print(f"  {n:2d} :   {rp:3d}    {rn:3d}  |      {genus}            {eis}          |   {cusp:+d}")
print(f"  EISENSTEIN identity (r_prin+2 r_non)/2 == sum_{{d|n}}(d/23) for all n? {eis_ok}  (= the genus/ideal count/main term)")
print(f"  CUSP difference r_prin - r_non is NONZERO ({cusp_nonzero} of 30 n).")

# ==========================================================================
sep("C  HISTORICAL / REFUTED DICTIONARY (retained to show MISTAKE-233)")
print("""  theta of a disc-D binary form  =  EISENSTEIN  (+)  CUSP :
    EISENSTEIN coeff = sum_{d|n} chi_D(d)  = the GENUS ideal count = the singular-series MAIN TERM = the FLOOR.
      -> computable from LOCAL Legendre data (S215); it is the genus part (S218) => ZERO deep entropy.
    CUSP coeff = the class-group deviations (r_i - r_j) = a weight-1 CM/dihedral cusp form.
      -> invisible to any congruence; the DEEP class group (S218) => the HIDDEN arithmetic entropy.
  h=1 (Heegner -3,-7,-11 = Paley 3,7,11): NO cusp forms => pure Eisenstein => zero deep entropy => RIGID.
  h>1 (-23,-47): cusp forms present => the deep fluctuation = where hidden information lives.""")

# ==========================================================================
sep("D  REFUTED LRC inheritance: THM-515 supplies no modular transformation law")
print("""  THM-515: L(S) = sum_{relation lattice} prod h  = a positive-definite lattice THETA (a modular form).
  By the same Eisenstein (+) cusp split:
    * the EISENSTEIN part = the singular-series MAIN TERM = the archimedean/genus expected loneliness = the
      FLOOR (the AP extremal value; the Fejer-Bochner SOS minorant; the Phi_6-Eisenstein argmax t*=14/183).
    * the CUSP part = the DEEP fluctuation of the singular series = the hidden arithmetic entropy = exactly
      WHERE a counterexample would have to hide (a cuspidal dip below the Eisenstein floor).
  LRC(14)=2*7 -> the apex weight-1 theta is at disc -7 (Paley-7), h(-7)=1 => PURE EISENSTEIN, no cusp =>
  the apex has ZERO deep entropy: the floor has no cuspidal dip to exploit at the 7-part. The remaining
  fluctuation is the non-apex (composite) structure -- consistent with kps-S17's 'the apex prime's arithmetic
  is the floor' and my S216/S217 (corrected): the RIGIDITY of 7 = no cusp forms at the apex = Eisenstein-only.""")

sep("HISTORICAL SUMMARY -- NOT A VALID LRC CONCLUSION")
print("""  Previous modular-form work ties to cutting-edge LRC through ONE split: theta = EISENSTEIN (+) CUSP.
   * VERIFIED: disc -7 (Paley-7, LRC(14)) theta = pure Eisenstein (r_Q=2 sum(d/7), no cusp) = h=1 rigid;
     disc -23 (h=3) theta = Eisenstein genus (sum(d/23)) + a NONZERO cusp (the deep class-group part).
   * The Eisenstein/cusp split IS the S218 genus/deep arithmetic-entropy split: Eisenstein = local floor
     (zero deep entropy), cusp = the hidden deep entropy.
   * The LRC singular series (THM-515) is a theta, so it inherits the split: Eisenstein = the FLOOR (the AP
     value / Fejer-Bochner / Eisenstein argmax), cusp = the deep fluctuation where a counterexample hides.
   * Apex Heegner h=1 (disc -7) => no cusp at the 7-part => the apex floor is cuspless => the RIGIDITY of 7.
  So: modular forms give the LRC floor (Eisenstein) and locate the difficulty (cusp = deep entropy), and the
  Paley/Heegner apex-prime rigidity is the statement that the apex weight-1 theta is Eisenstein-only.""")
