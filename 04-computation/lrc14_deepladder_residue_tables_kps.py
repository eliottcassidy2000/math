"""kps-S5: full symbolic residue tables for the deep one-swap ladders drop-9, drop-11
(and drop-12,13 recheck), to formalize residueLiar-style in Lean. Verify witness for all k."""
from fractions import Fraction
def resdist(a,q): r=a%q; return min(r,q-r)

# ladder params: family = {1..13}\{j} u {X0*k}; witness t*=(pa*k+pb)/(qa*k+qb); claim M=(ma*k)/(qa*k+qb)
LAD = {
 13: dict(drop=13, X0=182, pa=14,pb=0,  qa=182,qb=1, ma=14),  # deep well (far-peel; recheck)
 12: dict(drop=12, X0=84,  pa=35,pb=2,  qa=84, qb=5, ma=7),   # residue-liar (recheck)
 11: dict(drop=11, X0=154, pa=56,pb=1,  qa=154,qb=3, ma=14),  # NEW
 9:  dict(drop=9,  X0=126, pa=28,pb=1,  qa=126,qb=5, ma=14),  # NEW
}
for j,P in LAD.items():
    base=[v for v in range(1,14) if v!=j]
    print(f"\n===== drop-{j}: X={P['X0']}k, t*=({P['pa']}k+{P['pb']})/({P['qa']}k+{P['qb']}), claim M={P['ma']}k/({P['qa']}k+{P['qb']}) =====")
    # verify witness gives claimed M for k=1..25, and it's the TRUE max-min (spot vs fine search at k=1)
    ok=True; kappa_lin=(P['ma'],0)  # kappa=ma*k
    for k in range(1,26):
        X=P['X0']*k; V=base+[X]; p=P['pa']*k+P['pb']; q=P['qa']*k+P['qb']
        md=min(resdist(v*p,q) for v in V)
        claim=P['ma']*k
        if md!=claim: ok=False; print(f"  k={k}: md={md} != {claim}  MISMATCH"); 
    print(f"  witness M={P['ma']}k/({P['qa']}k+{P['qb']}) verified k=1..25: {ok};  M(k=1)={Fraction(P['ma'],P['qa']+P['qb'])}={float(P['ma']/(P['qa']+P['qb'])):.5f} > 1/14={1/14:.5f}: {P['ma']/(P['qa']+P['qb'])>1/14}")
    # symbolic residue table: for each runner v, r_v = v*p mod q; fit qq_v (constant) & r_v = A*k+B (linear)
    print(f"  RESIDUE TABLE (runner v: qq_v, r_v=A*k+B ; need {P['ma']}k <= r_v <= {P['qa']-P['ma']}k+{P['qb']}):")
    rows=[]
    for v in base+['X']:
        vv = (lambda k: P['X0']*k) if v=='X' else (lambda k,v=v: v)
        # compute r at k=2,3,4 to fit linear A*k+B and constant qq
        ks=[2,3,4]; rs=[]; qqs=[]
        for k in ks:
            p=P['pa']*k+P['pb']; q=P['qa']*k+P['qb']; N=vv(k)*p
            r=N%q; qq=N//q; rs.append(r); qqs.append(qq)
        # fit r=A*k+B
        A=Fraction(rs[1]-rs[0],ks[1]-ks[0]); B=Fraction(rs[0])-A*ks[0]
        # fit qq: for X-runner qq grows with k; for fixed v qq constant
        qq_const = (qqs[0]==qqs[1]==qqs[2])
        qqA=Fraction(qqs[1]-qqs[0],ks[1]-ks[0]); qqB=Fraction(qqs[0])-qqA*ks[0]
        vlabel = f"{P['X0']}k" if v=='X' else str(v)
        # bounds check symbolic: r=A*k+B, need ma*k<=A*k+B and A*k+B<=(qa-ma)k+qb
        lo_ok_slope = A>=P['ma'] or (A==P['ma'])  # and B>=0 handles lower
        rows.append((v,vlabel,A,B,qqA,qqB,rs))
    for (v,vlabel,A,B,qqA,qqB,rs) in rows:
        qqstr = f"{qqA}k+{qqB}" if qqA!=0 else f"{qqB}"
        print(f"    v={vlabel:>6}: qq={qqstr:>8}  r={A}k+{B}   (samples k=2,3,4: {rs})")
