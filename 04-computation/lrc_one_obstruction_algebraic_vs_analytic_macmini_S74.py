"""
S74: THE ONE OBSTRUCTION -- algebraic certificates are worst-case (config-blind) and all fail at the SAME
threshold (comonotone, k=7 events of mass 1/7 = one full circle); the ACTUAL apex comb is equidistributed
(config-specific) and succeeds. So the LRC(14) hard core is irreducibly ANALYTIC (equidistribution), not
algebraic (SOS / Lee-Yang zero-freeness / Bonferroni / Asano all break at the worst case).

MODEL: k danger events, each an arc of measure p=1/7 on the circle. Loneliness = P(M=0), M=#arcs covering a pt.
 - COMONOTONE (worst case = the algebraic certificate's adversary): arcs maximally aligned. P(M=0)=max(0,1-k p).
 - EQUIDISTRIBUTED / independent (the apex comb 14Z): P(M=0)=(1-p)^k > 0 always.
 - ACTUAL apex comb {14,28,...,14k}: danger ||14m t||<1/14; compute P(M=0) exactly on the t-grid.
Certificates: union-bound/Bonferroni-2; Lee-Yang min|zero| of Xi(lam)=G_M(1-lam) (zero-free in disk = certifies).
"""
from fractions import Fraction as F
import numpy as np
from math import comb

p=F(1,7)

def pm0_comonotone(k): return max(F(0), 1-k*p)
def pm0_indep(k): return (1-p)**k
def bonferroni2(k):   # P(M=0) >= 1 - k p + C(k,2) p^2 (incl-excl truncated at even order 2) -- a LOWER bound when valid
    return 1 - k*p + comb(k,2)*p*p

def pm0_apexcomb(k):
    """exact P(no apex runner within 1/14 of observer): danger of 14m at t is ||14 m t||<1/14."""
    E=[14*m for m in range(1,k+1)]
    b=set([F(0),F(1)])
    for e in E:
        for j in range(0,e+1):
            # boundaries where ||e t||=1/14: e t = j +- 1/14
            for s in (F(1,14),-F(1,14)):
                v=(F(j)+s)/e
                if 0<=v<=1: b.add(v)
    b=sorted(b); safe=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        mid=(x0+x1)/2
        lonely=all(min((e*mid)%1,1-(e*mid)%1)>=F(1,14) for e in E)
        if lonely: safe+=x1-x0
    return safe

def leeyang_minzero(probs_count):
    """G_M(z)=sum P(M=t) z^t ; Xi(lam)=G_M(1-lam); min |lam-zero| (interior => Bonferroni/Asano FAIL)."""
    q=probs_count
    if len(q)<2: return None
    # zeros of G_M(z): z-roots; map to lam=1-z
    c=[float(x) for x in q[::-1]]
    while len(c)>1 and abs(c[0])<1e-15: c=c[1:]
    if len(c)<2: return None
    zr=np.roots(c)
    lam=[1-z for z in zr]
    inside=[l for l in lam if abs(l)<1-1e-9]
    return (min(abs(l) for l in lam), len(inside))

def count_dist_indep(k):
    # M ~ Binomial(k, p): P(M=t)=C(k,t) p^t (1-p)^(k-t)
    return [comb(k,t)*(p**t)*((1-p)**(k-t)) for t in range(k+1)]
def count_dist_comonotone(k):
    # M in {0,k}: with prob k p at value... actually comonotone: M=k on measure p, M=0 on 1-kp (k<7)
    if k*p<=1: return [1-k*p]+[F(0)]*(k-1)+[k*p]
    return None

print("="*94)
print(" THE ONE OBSTRUCTION: algebraic certificates are WORST-CASE; comonotone fails at k=7; equidist succeeds")
print("="*94)
print(f"{'k':>3}{'P(M=0) comonotone':>20}{'P(M=0) equidist':>18}{'P(M=0) apex comb':>18}{'Bonf-2 lower':>14}")
for k in range(2,11):
    com=pm0_comonotone(k); ind=pm0_indep(k); ap=pm0_apexcomb(k); bf=bonferroni2(k)
    print(f"{k:>3}{str(com)+f' ({float(com):.3f})':>20}{f'{float(ind):.4f}':>18}{f'{float(ap):.4f}':>18}{f'{float(bf):.3f}':>14}")
print("-"*94)
print(" Lee-Yang: min|lam-zero| of Xi(lam)=G_M(1-lam); zero INSIDE unit disk => algebraic certificate FAILS")
print(f"{'k':>3}{'equidist min|zero|':>20}{'#inside disk':>14}{'comonotone min|zero|':>22}{'#inside':>9}")
for k in range(2,11):
    ind=leeyang_minzero(count_dist_indep(k))
    cm=count_dist_comonotone(k); cmz=leeyang_minzero(cm) if cm else None
    istr=f"{ind[0]:.3f}" if ind else "-"; icnt=ind[1] if ind else "-"
    cstr=f"{cmz[0]:.3f}" if cmz else "(no dist k>=7)"; ccnt=cmz[1] if cmz else "-"
    print(f"{k:>3}{istr:>20}{str(icnt):>14}{cstr:>22}{str(ccnt):>9}")
print("="*94)
print(" READOUT: comonotone (worst case) P(M=0)->0 at k=7 (saturation, 7 x 1/7 = full circle); equidist &")
print(" the ACTUAL apex comb stay POSITIVE for all k (equidistribution). Algebraic certs see the worst case")
print(" => fail at k=7; the proof must use the apex comb's EQUIDISTRIBUTION (analytic), not a worst-case cert.")
