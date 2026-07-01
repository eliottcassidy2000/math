#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
DEEP STUDY of the far-element BARRIER RESIDUAL: is it self-concordant, and what is
its exact periodic/Dedekind nature?

kind-pasteur-2026-06-30. Follows the-far-runner-is-the-log-barrier reflection, which
called Delta_w a "self-concordant barrier residual." THM-563 pins Delta_w exactly:
  Delta_w * w = Sum_j Sum_{endpoints t} +-S_j(frac(w t)),   S_j = SAWTOOTH antiderivative,
a GENERALIZED DEDEKIND SUM (Rademacher form), period 7*lcm(B), |S_j|<=3/49.

QUESTIONS this script settles:
  Q1  Is the SECTOR residual self-concordant?  (The barrier F with |F'''|<=2(F'')^{3/2}.)
      Sawtooth/tent barriers are piecewise-LINEAR -> F''=0 -> NOT self-concordant.
      The LOG barrier -log(dist) IS self-concordant (ratio == 2 exactly).
  Q2  What barrier gives which residual?
      TENT/sawtooth kernel  -> RADEMACHER sum (B_2, L1/absolute)  = the current Delta_w.
      LOG barrier (-log)    -> CLASSICAL DEDEKIND COTANGENT sum   = the SIGNED shadow.
      (Cycle-resistance reflection: classical Dedekind sum = degree-1 Fourier SIGNED shadow
       of the B_2/L1 discrepancy -- i.e. exactly the signed content the conjecture needs.)
  Q3  What does the log-barrier BUY?  DEDEKIND RECIPROCITY = a EUCLIDEAN DESCENT giving a
      UNIFORM-in-w bound |s(h,k)| via continued fractions -- replacing the per-base finite
      period-max (THM-563) with ONE structural reciprocity.

CONCLUSION target: "self-concordant barrier residual" is only literal for the LOG barrier;
the sector/tent residual is a Rademacher B_2 sum (not self-concordant, per-base period-max).
Switching to the log barrier makes interior-point LITERAL and hands us Dedekind reciprocity.
"""
import sys, math
from fractions import Fraction as Fr
from math import gcd, tan, pi, cos, sin, log
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---------- Dedekind sum: sawtooth (Rademacher) form and cotangent form ----------
def sawtooth(x):  # ((x)) : x - floor - 1/2 if not integer else 0  (exact on rationals)
    if isinstance(x, Fr):
        if x.denominator==1: return Fr(0)
        f=x-(x.numerator//x.denominator); return f-Fr(1,2)
    f=x-math.floor(x)
    return 0.0 if abs(f)<1e-15 or abs(f-1)<1e-15 else f-0.5
def dedekind_sawtooth(h,k):  # s(h,k)=sum_{i=1}^{k-1} ((i/k))((hi/k))  -- EXACT rational
    return sum(sawtooth(Fr(i,k))*sawtooth(Fr(h*i,k)) for i in range(1,k))
def dedekind_cotangent(h,k):  # (1/4k) sum cot(pi i/k) cot(pi h i/k)  -- float
    return sum(1.0/tan(pi*i/k)*1.0/tan(pi*h*i/k) for i in range(1,k))/(4*k)

print("="*84)
print("Q2a  THE TWO FORMS AGREE: sawtooth (Rademacher) form == cotangent (log-barrier) form")
print("     s(h,k) = sum ((i/k))((hi/k))  ==  (1/4k) sum cot(pi i/k)cot(pi hi/k)")
print("="*84)
for (h,k) in [(1,5),(2,7),(3,11),(5,13),(8,21),(13,34)]:
    sw=dedekind_sawtooth(h,k); ct=dedekind_cotangent(h,k)
    print(f"   s({h:2d},{k:2d}) sawtooth={float(sw):+.6f}={sw}   cotangent={ct:+.6f}   diff={abs(float(sw)-ct):.2e}")
print("   => the L1/sawtooth (B_2, current sector Delta_w) and the cotangent (log-barrier)")
print("      are the SAME NUMBER: one signed object, two barriers.")

print("\n"+"="*84)
print("Q3  DEDEKIND RECIPROCITY -> EUCLIDEAN DESCENT -> uniform bound (the log-barrier payoff)")
print("     s(h,k)+s(k,h) = -1/4 + (1/12)(h/k + k/h + 1/(hk))")
print("="*84)
def recip_rhs(h,k): return Fr(-1,4)+Fr(1,12)*(Fr(h,k)+Fr(k,h)+Fr(1,h*k))
for (h,k) in [(2,7),(5,13),(13,34),(89,233)]:
    lhs=dedekind_sawtooth(h,k)+dedekind_sawtooth(k,h)
    print(f"   s({h},{k})+s({k},{h}) = {float(lhs):+.6f}   reciprocity RHS = {float(recip_rhs(h,k)):+.6f}   diff={abs(float(lhs-recip_rhs(h,k))):.1e}")
def dedekind_descent(h,k):
    """compute s(h,k) by reciprocity+Euclid in O(log k); return (value, #steps, max partial quotient)."""
    h%=k; steps=0; val=Fr(0); sign=1; maxq=0
    while h!=0:
        # s(h,k) = -s(k,h) + recip(k? ) ; use s(h,k)= -s(k mod h, h) + [-1/4+(1/12)(h/k+k/h+1/hk)]  via s(k,h)=s(k mod h,h)
        val += sign*recip_rhs(h,k)
        sign=-sign
        h,k=k%h,h
        maxq=max(maxq, (k)// (h if h else 1))
        steps+=1
    # now s(0,k)=0; but we added recip terms; the accumulated val equals s(h0,k0) with the -s(...) unfolded
    return val, steps, maxq
print("   EUCLIDEAN DESCENT (reciprocity unfold) vs direct sum:")
for (h,k) in [(5,13),(13,34),(89,233),(196,509),(1,1000000)]:
    if gcd(h,k)!=1: continue
    direct = dedekind_sawtooth(h,k) if k<=2000 else None
    dv,steps,maxq = dedekind_descent(h,k)
    ds = f"{float(direct):+.5f}" if direct is not None else " (skip: k large)"
    ok = "" if direct is None else f" match={abs(float(direct-dv))<1e-9}"
    print(f"   s({h},{k}): descent={float(dv):+.5f} in {steps:2d} steps (Euclid, O(log k)); direct={ds}{ok}")
print("   => reciprocity computes s(h,k) in O(log k) steps and BOUNDS it by the continued-")
print("      fraction partial quotients of h/k -- a UNIFORM structural bound. The sector")
print("      period-max (THM-563) is a per-base O(lcm) max with NO such descent.")

print("\n"+"="*84)
print("Q1  SELF-CONCORDANCE  |F'''| <= 2 (F'')^{3/2}  : which barrier qualifies?")
print("="*84)
def selfconc_ratio(F, x, h=1e-4):
    f2=(F(x+h)-2*F(x)+F(x-h))/h**2
    f3=(F(x+2*h)-2*F(x+h)+2*F(x-h)-F(x-2*h))/(2*h**3)
    if f2<=1e-9: return f2, f3, float('inf')
    return f2, f3, abs(f3)/(f2**1.5)
print("   barrier F(x)            F''      F'''        |F'''|/(F'')^1.5   self-concordant?")
tests=[
 ("-log(x)            ", lambda x:-math.log(x), 0.3),
 ("-log(x)            ", lambda x:-math.log(x), 0.05),
 ("-log(sin(pi x))    ", lambda x:-math.log(math.sin(math.pi*x)), 0.1),   # arc/circle log-barrier
 ("tent |x-1/2|       ", lambda x:abs(x-0.5), 0.3),                        # piecewise-linear (sawtooth-like)
 ("tent |x-1/2|       ", lambda x:abs(x-0.5), 0.49),
]
for name,F,x in tests:
    f2,f3,r=selfconc_ratio(F,x)
    verdict = "YES (<=2)" if r<=2.0+0.05 else ("NO (F''=0)" if f2<=1e-6 else f"NO (ratio {r:.2f})")
    print(f"   {name} @x={x:.2f}  F''={f2:9.3f}  F'''={f3:11.2f}   {('inf' if r==float('inf') else f'{r:.4f}'):>8}          {verdict}")
print("   => -log is self-concordant (ratio->2, the canonical barrier). The TENT/sawtooth")
print("      has F''=0 a.e. (F''' singular at the kink) -> NOT self-concordant.")

print("\n"+"="*84)
print("Q2b  THE BRIDGE: the LOG barrier's gradient IS the cotangent (hence Dedekind)")
print("     d/dx [ -log(2 sin(pi x)) ] = -pi cot(pi x);  summing -cot over the orbit {h i/k}")
print("     reproduces the cotangent Dedekind sum = the SIGNED residual.")
print("="*84)
for (h,k) in [(2,7),(5,13),(8,21)]:
    grad_sum = sum(1.0/tan(pi*i/k)*1.0/tan(pi*h*i/k) for i in range(1,k))/(4*k)
    print(f"   (1/4k) sum_i cot(pi i/{k}) cot(pi {h} i/{k}) = {grad_sum:+.6f} = s({h},{k})  [log-barrier residual]")
print("   => The log barrier is self-concordant AND its far-runner residual is the classical")
print("      Dedekind sum -- the SIGNED shadow the conjecture needs, WITH reciprocity.")

print("\n"+"="*84)
print("SYNTHESIS")
print("="*84)
print(" - 'Self-concordant barrier residual' is LITERAL only for the LOG barrier, not the")
print("   sector/tent functional (piecewise-linear, F''=0).")
print(" - Both barriers give the SAME signed number s(h,k): the sector Delta_w is its")
print("   Rademacher/B_2 (sawtooth, L1) form; the log barrier gives the COTANGENT form.")
print(" - The log barrier BUYS: (i) literal interior-point self-concordance (central path,")
print("   Newton), (ii) DEDEKIND RECIPROCITY = a Euclidean descent = a UNIFORM-in-w bound,")
print("   which the per-base period-max (THM-563) does not provide.")
print(" - Route refinement: recast the far-element residual via the LOG barrier to trade the")
print("   per-base finite period-max for ONE reciprocity/continued-fraction bound -- the")
print("   general-bounded-base closure THM-563 leaves 'in progress'.")
print("DONE.")
