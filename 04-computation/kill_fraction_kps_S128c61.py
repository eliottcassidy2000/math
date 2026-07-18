#!/usr/bin/env python3
"""kill_fraction_kps_S128c61.py -- kind-pasteur S128 cont.61.
REPLACING THE FINITE HORN BY A COUNTING LEMMA, uniformly in r.
A family fails the small-modulus criterion iff, for EVERY (q,a) safe for the core, some
killer is unsafe.  Write bits(P) = {(q,a) safe for the core} and kill(k) = bits(P) minus
k's safe set.  Then the family is UNCERTIFIED iff the kill-sets COVER bits(P).  Hence

      sum_i |kill(k_i)|  <  |bits(P)|      ==>   CERTIFIED  (no enumeration needed).

Each killer is unsafe at a (q,a) with density about 2*ceil(q/14)/q ~ 1/7, so r killers
should cover at most ~r/7 -- but a killer DIVISIBLE by q is unsafe at that q for EVERY a,
which is what inflates |kill|.  Measure the worst kill-fraction.  PRINT DATA ONLY."""
import sys, itertools
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
def corebits(P):
    return [i for i,(q,a) in enumerate(QS) if all(la(p*a,q)>=-(-q//14) for p in P)]
def killfrac(bits,k):
    bad=0
    for i in bits:
        q,a=QS[i]
        if la(k*a,q)<-(-q//14): bad+=1
    return bad/len(bits) if bits else 1.0
print("### worst kill-fraction over killers, by core size ###")
print("  |P|  #cores  min|bits|  worst kill-frac   worst killer   r_max = floor(1/frac)")
for size,KB in [(11,900),(10,900),(9,900)]:
    CS=[sorted(c) for c in itertools.combinations(range(1,13),size)]
    worst=0.0; argk=None; minbits=10**9
    for P in CS:
        bits=corebits(P)
        minbits=min(minbits,len(bits))
        if not bits: continue
        M=max(P)
        for k in range(13*M+1,KB):
            f=killfrac(bits,k)
            if f>worst: worst=f; argk=(tuple(P),k)
    print("  %-4d %-7d %-10d %-18.5f %-14s %d"%(
        size,len(CS),minbits,worst,str(argk[1]) if argk else "-",int(1/worst) if worst>0 else -1))
print()
print("### which killers have the biggest kill-fraction?  (divisibility is the driver) ###")
P=[x for x in range(1,13) if x!=1]; bits=corebits(P)
rows=[]
for k in range(157,900):
    rows.append((killfrac(bits,k),k,[q for q in range(15,41) if k%q==0]))
rows.sort(reverse=True)
print("  core {2..12}, |bits| = %d"%len(bits))
print("  frac     killer  divisors in [15,40]")
for f,k,dv in rows[:10]:
    print("  %-8.5f %-7d %s"%(f,k,dv))
print("  ...")
for f,k,dv in rows[-3:]:
    print("  %-8.5f %-7d %s"%(f,k,dv))
print()
print("### the counting lemma verdict ###")
for size in [11,10,9]:
    CS=[sorted(c) for c in itertools.combinations(range(1,13),size)]
    worst=0.0
    for P in CS:
        bits=corebits(P)
        if not bits: continue
        M=max(P)
        for k in range(13*M+1,900):
            f=killfrac(bits,k)
            worst=max(worst,f)
    r=13-size
    print("  r=%d killers, core size %d: worst frac %.5f ; r*frac = %.4f  -> %s"%(
        r,size,worst,r*worst,"CERTIFIED by counting" if r*worst<1 else "counting insufficient"))
print("DONE")
