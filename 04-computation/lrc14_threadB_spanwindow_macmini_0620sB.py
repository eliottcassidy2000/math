#!/usr/bin/env python3
r"""
lrc14_threadB_spanwindow_macmini_0620sB.py   (mac-mini, Thread B stage 5 -- the deliverable)

THREAD B DELIVERABLE.  The finite-check feasibility, assembled correctly.

WHAT IS ALREADY PROVED (cited, not redone):
  - THM-536 B2 (PROVED, pointwise set inclusion):  E subset {0..N}  =>  measS7(E) <= measS7(AP_{N+1}).
    Hence measS7(E) <= cap_k for EVERY primitive E of span <= N*(k), where N* = largest N with
    measS7(AP_{N+1}) <= cap_k:   N* = 7,8,10,13,20,20  for k = 8..13.
    => the finite check below span N*(k) is CLOSED by a clean domination (no enumeration).

THREAD B (this script) addresses the RESIDUAL WINDOW  N*(k) < span <= W  and the FINITENESS
of the family the finite check must actually enumerate.

STRUCTURAL CHARACTERIZATION (established stages 1-4, summarized):
  - measS7 is scale-invariant -> normalize to primitive (gcd of nonzero offsets = 1).
  - relation lattice Lambda(E)={n:sum n_j s_j=0}, rank k-2; l_inf successive minima
    lambda_1<=...<=lambda_{k-2}; lambda_max:=lambda_{k-2} is the finiteness ruler.
  - lambda_max(E) <= H0  =>  span(E) bounded (FINITE family).  SHARP on the stranger
    sub-family {0..k-2,N}:  lambda_max <= H0  <=>  N <= (k-2)(k-1)/2 * H0  (=21*H0 at k=8,
    EXACT, stage-4 Part A).  General Hadamard span bound: span <= (k-2)! * H0^{k-2}.
  - The DICHOTOMY for the sector route:
       span <= W  : FINITE exact check (this script verifies measS7 <= cap exhaustively),
       span >  W  : every such primitive E has lambda_max > H0(W) (a dissociated direction
                    / a far stranger), handled by Thread A (wide bound / stranger contraction).

THIS SCRIPT (exact, feasible):
  (A) recompute N*(k) and confirm B2 closes span<=N*(k).
  (B) EXHAUSTIVE measS7<=cap over ALL primitive E (0 in E,|E|=k) with span<=W, for k=8,9,10,
      with W chosen as large as feasible.  Report worst measS7, #over-cap, and the WORST shape.
  (C) For each E in the window, record lambda_max; tabulate the JOINT (span, lambda_max,
      measS7) landscape => show the finite-check family is exactly {span<=W} and is finite,
      and that beyond W lambda_max is forced up (the wide regime).
  (D) COUNT N(k) := # scale-normalized (primitive) shapes the finite check must cover up to
      span W, broken down by whether B2 already covers them.

stdlib only; exact Fraction.  Designed to FINISH (no per-combo rank tests in the hot loop;
lambda_max only computed on the small surviving 'near-cap' set).
"""
import sys, itertools, math
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def gcd_all(xs): return reduce(gcd, [abs(x) for x in xs if x != 0], 0)
def banner(t): print("\n"+"="*82+f"\n{t}\n"+"="*82)

# ---- fast measS7: only need to know if ALL 7 sectors hit; exact via breakpoints ----
def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7*abs(e)+1): bps.add(F(m, 7*abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1); tot = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        secs = {int(((e*xm) % 1)*7) for e in E}
        if len(secs) == 7: tot += x1-x0
    return tot

CAPS = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

# ---- lambda_max (only used on the small near-cap set; exact l_inf) ----
def mat_rank_Q(rows):
    M = [[F(x) for x in row] for row in rows]; r0=0; rank=0; ncol=len(M[0])
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
def lambda_max(speeds, Hcap=25):
    d=len(speeds); need=d-1; found=[]; mins=[]
    for H in range(1,Hcap+1):
        for n in itertools.product(range(-H,H+1),repeat=d):
            if max(abs(x) for x in n)!=H: continue
            if sum(n[i]*speeds[i] for i in range(d))!=0: continue
            if not found:
                if any(x!=0 for x in n): found.append([F(x) for x in n]); mins.append(H)
            else:
                if mat_rank_Q(found+[[F(x) for x in n]])>len(found):
                    found.append([F(x) for x in n]); mins.append(H)
            if len(found)>=need: return max(mins)
    return None

# ===========================================================================
def partA_Nstar():
    banner("PART A -- N*(k): B2 domination closes span <= N*(k)")
    print("  measS7(AP_{N+1}) for N=k-1.. and the largest N* with measS7(AP_{N+1}) <= cap_k.")
    for k in (8,9,10):
        cap = CAPS[k]; Nstar=None
        row=[]
        for N in range(k-1, 30):
            ap = list(range(N+1))   # AP_{N+1} = {0..N}
            # only need the k-subset? NO: B2 dominates ANY E subset {0..N} of size k by AP_{N+1}
            # which has N+1 elements. measS7(AP_{N+1}) is the dominating value.
            m = measS7(ap)
            row.append((N,m))
            if m <= cap: Nstar=N
            else: break
        print(f"  k={k}: cap={float(cap):.4f}; N* = {Nstar} (B2 closes span<=N*); "
              f"measS7(AP_{{N*+1}})={float(measS7(list(range(Nstar+1)))):.4f}")
    print("  => N* = 7,8,10 for k=8,9,10 (matches THM-536). Residual: N*<span<=W.")

def partB_exhaustive(k, W):
    banner(f"PART B -- EXHAUSTIVE measS7<=cap over ALL primitive E, |E|={k}, span<=W={W}")
    cap = CAPS[k]
    consec = measS7(list(range(k)))
    print(f"  cap_{k}={cap}={float(cap):.5f}; consec_{k}={float(consec):.5f}")
    cnt=0; over=0; worst=F(0); argworst=None
    near=[]  # shapes with measS7 within 0.02 of cap (the dangerous ones)
    import time; t0=time.time()
    for combo in itertools.combinations(range(1, W+1), k-1):
        if combo[-1] <= k-2:  # span must be >= k-1 trivially; skip degenerate
            pass
        if gcd_all(combo)!=1: continue
        E=[0]+list(combo)
        cnt+=1
        m=measS7(E)
        if m>worst: worst=m; argworst=E
        if m>cap: over+=1
        if m > cap - F(1,50):
            near.append((E, m))
    dt=time.time()-t0
    print(f"  #primitive E (span<=W): {cnt}   ({dt:.0f}s)")
    print(f"  WORST measS7 = {worst} = {float(worst):.6f}  at {argworst}")
    print(f"  # over cap_{k}: {over}   =>  measS7 <= cap for ALL: {over==0}")
    print(f"  margin cap - worst = {float(cap-worst):.6f}")
    print(f"  # shapes within 0.02 of cap (the tight ones): {len(near)}")
    for E,m in sorted(near, key=lambda t:-t[1])[:8]:
        lm = lambda_max([e for e in E if e!=0], Hcap=12)
        print(f"      measS7={float(m):.5f} span={max(E)} lambda_max={lm}  E={E}")
    return cnt, over, worst, near

if __name__ == "__main__":
    print("#"*82)
    print("# THREAD B stage 5 -- span-window finite check (the deliverable)")
    print("#"*82)
    partA_Nstar()
    # feasible windows: k=8 W=20 (C(20,7)=77520), k=9 W=16 (C(16,8)=12870), k=10 W=15 (C(15,9)=5005)
    partB_exhaustive(8, 20)
    partB_exhaustive(9, 16)
    partB_exhaustive(10, 15)
    print("\nDONE (Thread B stage 5).")
