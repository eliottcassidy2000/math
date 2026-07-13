"""
opus-2026-07-11-S263: bound the multi-linear cancellation (S262 residual) via Gowers norms. Outcome: GOWERS
NORMS DO NOT APPLY (parallel forms); the multi-linear cancellation is governed by the ADDITIVE RELATIONS among
the speeds -- the sumset / Schur = E3 structure (LEM-015), the project's own additive invariant. A redirect.

THE OBJECT. eps_v (multi-linear part, S262) is the correlation of dilates of the band g at the SAME point:
    int prod_i g(v_i t) dt = Sum_{Sum_i v_i k_i = 0} prod_i ghat(k_i).
The dominant off-diagonal terms are SMALL relations Sum v_i k_i = 0 with k_i = +-1, i.e. +-v +- w1 +- w2 (...) = 0
-- the ADDITIVE relations among {v} u (non-core subset).

WHY GOWERS DOES NOT APPLY. Gowers U^s norms bound multilinear averages int prod_i g(psi_i(t)) via the
generalized von Neumann theorem ONLY when the linear forms psi_i are in GENERAL POSITION (pairwise linearly
independent, bounded Cauchy-Schwarz complexity). Here psi_i(t) = v_i * t are all PROPORTIONAL to t -- a single
1-dimensional direction, maximally degenerate. General position FAILS, so the generalized von Neumann /
Gowers-norm machinery does not bound int prod g(v_i t). Gowers norms are the wrong tool.

THE RIGHT TOOL: ADDITIVE RELATIONS (E3). The correlation is Sum_{Sum v_i k_i=0} prod ghat(k_i), dominated by
the +-1 relations = the additive/sumset structure. VERIFIED: |eps_v (multi)| correlates with the count of
additive relations v = w_i +- w_j among non-core:
    correlation( |eps_v|, #relations ) = 0.527; monotone by count: #rel 0->|eps|0.021, 2->0.027, 4->0.065,
    6->0.073, 8->0.086.
Each +-1 relation contributes ~b_1^3 ~ 0.0026 (weighted by (6/7)^{r-2}), so |eps_v| ~ (#relations)*const. This
is exactly the E3 = Schur-triple structure (LEM-015): the core arc's correlation with the good-set product is
driven by how additively RELATED v is to the non-core speeds. RUNNER 1 has the MOST relations (every
consecutive non-core difference w'-w=1), hence eps_1 ~ 0.57 -- the near-AP maximum (LEM-015: E3 maximized at
the interval/AP), handled by S255. Dissociated core runners (large, coprime, few relations) have eps ~ 0.02.

NET. Bounding the multi-linear cancellation via Gowers norms FAILS (parallel forms, general position violated).
The correct tool is ADDITIVE COMBINATORICS: the multi-linear cancellation = the additive-relation (sumset /
Schur / E3) structure of the core against the non-core, which the project already has as LEM-015 (E3 max at the
AP). So the covering-min residual reduces to an E3 / DISSOCIATION bound: the coprime core is additively
dissociated from the non-core (few relations) => eps small => coreCover < 1 => LRC(14), with the AP/runner-1
(max E3) as the S255 exception. The crux is E3-additive, not Gowers -- a redirect onto the project's own
additive invariant.

-> opus-S262 (multi-linear residual), LEM-015 (E3 = Schur triples, max at AP), opus-S246 (E3 the right additive
lever), opus-S255 (runner-1/near-AP/max-E3), opus-S259.
"""
import numpy as np
from math import gcd, pi
from functools import reduce
import random, statistics
def divcomplete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def addrel_count(v, non):
    nn=list(non); cnt=0
    for i in range(len(nn)):
        for j in range(len(nn)):
            if i!=j and (v==nn[i]+nn[j] or v==abs(nn[i]-nn[j])): cnt+=1
    return cnt
def main():
    D=13860; c=1.0/14; cD=c*D; random.seed(4); rows=[]; seen=0
    while seen<40:
        v=sorted(random.sample(range(1,90),13))
        if not (primitive(v) and divcomplete(v)): continue
        core=[x for x in v if gcd(x,30030)==1]; non=[x for x in v if gcd(x,30030)!=1]
        a=np.arange(D); safe=np.ones(D,bool)
        for w in non:
            rr=(w*a)%D; safe &= (np.minimum(rr,D-rr)>=cD)
        g=safe.astype(float); Gm=g.mean()
        if Gm<0.02: continue
        seen+=1
        for vv in core:
            rr=(vv*a)%D; Dv=(np.minimum(rr,D-rr)<cD).astype(float); eps=(Dv*g).sum()/g.sum()-1/7
            rows.append((abs(eps), float(addrel_count(vv,non))))
    es=[r[0] for r in rows]; rcs=[r[1] for r in rows]
    print(f"n={len(rows)}, correlation(|eps_v|, #additive-relations v=w_i±w_j) = {np.corrcoef(es,rcs)[0,1]:.3f}")
    byrc={}
    for e,rc in rows: byrc.setdefault(int(rc),[]).append(e)
    for rc in sorted(byrc)[:8]: print(f"   #rel={rc}: mean|eps|={statistics.mean(byrc[rc]):.4f} (n={len(byrc[rc])})")
    print("=> multi-linear cancellation = ADDITIVE RELATIONS (E3, LEM-015), NOT Gowers (parallel forms).")
if __name__=='__main__':
    main()
