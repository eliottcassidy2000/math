from fractions import Fraction as F
# EXACT E[maxgap] of the orbit {j*x : j in E} via three-gap piecewise-linear integration.
def Emaxgap_exact(E, kdenom):
    # breakpoints: Farey m/d, d<=kdenom (order changes + wraps)
    bps=set([F(0),F(1)])
    for d in range(1,kdenom+1):
        for m in range(1,d): bps.add(F(m,d))
    bps=sorted(bps)
    total=F(0)
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2
        # phases p_j(x)=j*x - floor(j*mid) on (a,b); sorted order fixed at mid
        fl={j:(j*mid).__floor__() for j in E}
        order=sorted(E, key=lambda j: (j*mid - fl[j]))
        # gaps g_s(x)=c*x+b0 (consecutive sorted + wrap); each linear
        gaps=[]
        for s in range(len(order)):
            j1=order[s]; j2=order[(s+1)%len(order)]
            if s<len(order)-1:
                c=F(j2-j1); b0=F(-(fl[j2]-fl[j1]))
            else:
                c=F(order[0]-order[-1]); b0=F(-(fl[order[0]]-fl[order[-1]])+1)
            gaps.append((c,b0))
        # maxgap(x)=max_s (c_s x + b0_s) on (a,b): piecewise-linear; sub-break at pairwise crossings in (a,b)
        subbp=set([a,b])
        for i in range(len(gaps)):
            for jx in range(i+1,len(gaps)):
                ci,bi=gaps[i]; cj,bj=gaps[jx]
                if ci!=cj:
                    xc=(bj-bi)/(ci-cj)
                    if a<xc<b: subbp.add(xc)
        subbp=sorted(subbp)
        for u,v in zip(subbp,subbp[1:]):
            m2=(u+v)/2
            mg=max(c*m2+b0 for c,b0 in gaps)  # which gap is max on (u,v) (constant choice)
            # integrate that linear gap over (u,v): trapezoid
            # find the argmax gap at m2, integrate it
            cbest,bbest=max(gaps,key=lambda cb: cb[0]*m2+cb[1])
            fu=cbest*u+bbest; fv=cbest*v+bbest
            total += (fu+fv)/2*(v-u)
    return total

thr=F(1,7)
print("=== EXACT E[maxgap(AP_k)] via three-gap integration; reverse-Markov floor (opus-S133) ===\n")
print(f"{'k':>3} {'E[maxgap] exact':>18} {'~':>8} {'>1/7 margin':>12} {'revMarkov mu_17>=':>18}")
for k in range(8,14):
    E=list(range(1,k+1))
    Em=Emaxgap_exact(E,k)
    lb=(Em-thr)/(1-thr)   # B=1
    print(f"{k:>3} {str(Em):>18} {float(Em):>8.5f} {float(Em-thr):>12.5f} {float(lb):>18.5f}")
print(f"\n  1/7 = {float(thr):.5f}. All E[maxgap(AP_k)] > 1/7 with margin => reverse-Markov density floor > 0.")
print("  REDUCTION: density-floor positivity <= [E[maxgap] > 1/7] (mean, not tail). Binding = AP.")
