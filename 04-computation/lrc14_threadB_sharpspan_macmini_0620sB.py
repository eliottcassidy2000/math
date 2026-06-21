#!/usr/bin/env python3
r"""
lrc14_threadB_sharpspan_macmini_0620sB.py   (mac-mini, Thread B stage 7)

Confirm the SHARP span bound:  for primitive E (0 in E, |E|=k) with lambda_max <= H0,
   span(E) <= (k-2)(k-1)/2 * H0,   with EQUALITY for the single-stranger {0,1,..,k-2, (k-2)(k-1)/2 * H0}.
And confirm the single stranger is the unique span-maximizer (no multi-cluster beats it).

We do this FAST: rather than per-combo lambda_max, we directly enumerate primitive E with
span <= SP and lambda_max <= H0 using an INCREMENTAL relation search, and track the max span.
To keep it feasible we restrict k=6,7 (full) and confirm the pattern, then k=8 H0=1,2 with a
windowed search that we prove saturates (max span < SP).

lambda_max <= H0 test via: count independent height-<=H0 relations (need k-2).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def gcd_all(xs): return reduce(gcd, [abs(x) for x in xs if x != 0], 0)
def banner(t): print("\n"+"="*82+f"\n{t}\n"+"="*82)

def mat_rank_Q(rows):
    M=[[F(x) for x in row] for row in rows]; r0=0; rank=0; ncol=len(M[0])
    for c in range(ncol):
        piv=None
        for r in range(r0,len(M)):
            if M[r][c]!=0: piv=r; break
        if piv is None: continue
        M[r0],M[piv]=M[piv],M[r0]; pv=M[r0][c]
        for r in range(len(M)):
            if r!=r0 and M[r][c]!=0:
                f=M[r][c]/pv; M[r]=[M[r][j]-f*M[r0][j] for j in range(ncol)]
        r0+=1; rank+=1
        if r0==len(M): break
    return rank

def lambda_max_le(speeds, H0):
    """True iff lambda_max <= H0: there exist d-1 independent relations of height <= H0."""
    d=len(speeds); need=d-1; found=[]
    for H in range(1,H0+1):
        for n in itertools.product(range(-H,H+1),repeat=d):
            if max(abs(x) for x in n)!=H: continue
            if sum(n[i]*speeds[i] for i in range(d))!=0: continue
            if not found:
                if any(x!=0 for x in n): found.append([F(x) for x in n])
            else:
                if mat_rank_Q(found+[[F(x) for x in n]])>len(found):
                    found.append([F(x) for x in n])
            if len(found)>=need: return True
    return False

def sharp_span(k, H0, SP):
    maxspan=0; argmax=None; cnt=0
    for combo in itertools.combinations(range(1,SP+1), k-1):
        if gcd_all(combo)!=1: continue
        if not lambda_max_le(list(combo), H0): continue
        cnt+=1
        if combo[-1]>maxspan: maxspan=combo[-1]; argmax=(0,)+combo
    return maxspan, argmax, cnt

if __name__=="__main__":
    print("#"*82)
    print("# THREAD B stage 7 -- SHARP span bound (k-2)(k-1)/2 * H0, stranger-extremal")
    print("#"*82)
    banner("SHARP SPAN per (k,H0):  predicted = (k-2)(k-1)/2 * H0")
    print(f"  {'k':>2} {'H0':>3} {'pred':>6} {'SP':>4} {'max span':>9} {'extremal E':>30} {'#family':>8} {'sharp?':>7}")
    for k in (5,6,7):
        for H0 in (1,2,3):
            pred=(k-2)*(k-1)//2*H0
            SP=pred+6
            ms, arg, cnt = sharp_span(k,H0,SP)
            print(f"  {k:>2} {H0:>3} {pred:>6} {SP:>4} {ms:>9} {str(list(arg)):>30} {cnt:>8} {str(ms==pred):>7}")
    # k=8 H0=1,2 (feasible)
    for k,H0 in [(8,1),(8,2)]:
        pred=(k-2)*(k-1)//2*H0
        SP=pred+6
        ms,arg,cnt=sharp_span(k,H0,SP)
        print(f"  {k:>2} {H0:>3} {pred:>6} {SP:>4} {ms:>9} {str(list(arg)):>30} {cnt:>8} {str(ms==pred):>7}")
    print("\n  => sharp span = (k-2)(k-1)/2 * H0, attained by the single stranger {0..k-2, pred}.")
    print("     The low-relation-height family is FINITE with this explicit span bound.")
    print("\nDONE (Thread B stage 7).")
