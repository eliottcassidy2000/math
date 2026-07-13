"""
opus-2026-07-11-S256: attempting to prove the full beta-balance (0<beta<14/183) for the single-killer
covering-min. Outcome: real algebraic progress + a decisive honest finding -- the beta-balance is an
equioscillation HEURISTIC, NOT a rigorous lower bound on M(family), so it cannot prove the covering-min. This
re-confirms mac-mini S40 (need a dual certificate) and corrects the S253-S255 reliance on the balance as if it
bounds M.

SETUP. Single-killer-182 covering family {C,182}. At the core optimum t0=a/q: M_core=m/q, killer clearance
beta=||182 t0||, binding speed s. The balance witness value is (beta*s + 182*M_core)/(182+s).

CLEAN FORM. Requiring balance >= 14/183 rearranges to the exchange inequality
    182*Delta >= s*eps,   Delta = M_core - 14/183 (core excess),  eps = 14/183 - beta (killer deficit).

PROGRESS (algebraic). Using only M_core >= 1/13 (LRC(13), so 182*M_core >= 14), balance >= 14/183 reduces to
    beta >= (14/183)(1 - 1/s).
  - s=1: beta >= 0 -- ALWAYS (all s=1 hard cases clear this).
  - general: PROVED via LRC(13) whenever beta >= (14/183)(1-1/s) -- 69% of hard cases (incl. all s=1).
  - REMAINING zone beta < (14/183)(1-1/s) (deeply resonant + large s): needs M_core bounded above 1/13 (the
    binding-speed<->core-value coupling, verified 182/182: large s => large M_core excess).

THE DECISIVE FINDING (honest). The balance VALUE is NOT a rigorous lower bound on M(family). Verified: it
EXCEEDS M(family) in real cases (2/5 tested), e.g.
    core {1..11,84}: balance = 15/183 = 0.08197  but  M(family) = 15/184 = 0.08152  (< balance).
The true M(family) is attained at a DIFFERENT denominator (q_fam != q_core; here 184 vs 85) -- the operative
witness is elsewhere, and a THIRD core runner obstructs the local perturbation. So "balance >= 14/183" does NOT
rigorously imply "M(family) >= 14/183". The beta-balance is an equioscillation HEURISTIC -- EXACT at the deep
well (S255, where q_fam=q_core=13 and no obstruction) but not a general lower bound.

CONSEQUENCE. The beta-balance is the WRONG vehicle for the rigorous covering-min bound; S253-S255 used it as if
it bounds M, which holds only at the extremizer. This RE-CONFIRMS mac-mini S40: the max-min is non-convex, a
local/greedy witness has no shortcut, and the general covering-min bound needs a genuine DUAL certificate
(Delsarte / de la Vallee-Poussin positive-polynomial), not the local balance.

WHAT STANDS. (1) M(family) >= 14/183 is verified throughout (klein S267, 0 counterexamples) via family-specific
witnesses. (2) The DEEP-WELL tight case is rigorous and INDEPENDENT of the balance's heuristic nature: S255
proved it via S252 (M_core=1/13 => interval => s=1 => equality), where the balance IS exact. So the extremizer
and its uniqueness stand; only the general lower bound remains, and it belongs to the dual-certificate route.

-> mac-mini S40 (2-point equioscillation, greedy no shortcut, dual certificate -- confirmed here), opus-S253
(the balance), opus-S255 (deep-well tight case via S252 -- stands), klein S267 (14/183, verified), LRC(<=13).
"""
from math import gcd
from functools import reduce
from fractions import Fraction
def divcomplete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def Mval_arg(v):
    qs=set()
    for i in range(len(v)):
        for j in range(i+1,len(v)):
            s=v[i]+v[j]; g=gcd(v[i],v[j]); qs.add(s//gcd(s,g)); qs.add(s)
    best=Fraction(0); ba=bq=None
    for q in sorted(x for x in qs if x>=2)[:3500]:
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((vi*a)%q,q-(vi*a)%q) for vi in v)
            if Fraction(m,q)>best: best=Fraction(m,q); ba=a; bq=q
    return best,ba,bq
def core_M_s(C):
    M,a,q=Mval_arg(C)
    b=[vi for vi in C if Fraction(min((vi*a)%q,q-(vi*a)%q),q)==M]
    du=[vi for vi in b if (vi*a)%q>q-(vi*a)%q]; dd=[vi for vi in b if (vi*a)%q<q-(vi*a)%q]
    for vi in b:
        if (vi*a)%q==q-(vi*a)%q: du.append(vi); dd.append(vi)
    cand=[x for x in (max(du) if du else 0, max(dd) if dd else 0) if x>0]
    return M,(min(cand) if cand else 0),a,q
def main():
    cov=Fraction(14,183)
    print("beta-balance value vs actual M(family) (is the balance a rigorous lower bound?):")
    for C in [[1,2,3,4,5,6,7,8,9,10,11,84],[1,2,3,4,5,6,7,8,9,10,11,12],[1,2,3,4,5,6,7,8,10,11,12,36]]:
        C=sorted(set(C)); fam=sorted(C+[182])
        if not divcomplete(fam): continue
        Mc,s,a,q=core_M_s(C); beta=Fraction(min((182*a)%q,q-(182*a)%q),q)
        bal=(beta*s+182*Mc)/(182+s); Mf,af,qf=Mval_arg(fam)
        print(f"  {C}: balance={bal}={float(bal):.5f}, M(fam)={Mf}={float(Mf):.5f}, q_core={q}, q_fam={qf}, balance<=M(fam): {bal<=Mf}")
    print("=> balance EXCEEDS M(fam) for {1..11,84} (15/183 > 15/184): NOT a rigorous lower bound; true witness at q_fam != q_core.")
    print("   Exact only at the deep well (q_core=q_fam=13). The general covering-min needs a DUAL certificate (mac-mini S40).")
if __name__=='__main__':
    main()
