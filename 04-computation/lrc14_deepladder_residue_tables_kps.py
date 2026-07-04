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

# --- drop-8 (X0=56=lcm(8,14)) and drop-10 (X0=70=lcm(10,14)): complete the unique-coverer hexad ---
print("\n\n########## drop-8, drop-10 (complete deep hexad) ##########")
def resd(a,q): r=a%q; return min(r,q-r)
for j,X0 in [(8,56),(10,70)]:
    base=[v for v in range(1,14) if v!=j]
    print(f"\n===== drop-{j}: X={X0}k =====")
    # find witness at k=1,2,3 by search, fit p,q,md linear
    from math import gcd
    data=[]
    for k in range(1,5):
        X=X0*k; V=base+[X]; best=(0,0,1,0)
        for q in range(2,200*k+80):
            for p in range(1,q):
                if gcd(p,q)!=1: continue
                md=min(resd(v*p,q) for v in V)
                if md*best[2]>best[0]*q: best=(md,p,q,md)
        md,p,q,_=best
        data.append((k,X,p,q,md))
        print(f"  k={k}: X={X} t*={p}/{q} M={md}/{q}={md/q:.5f} md={md}")
    # fit q=qa*k+qb, p=pa*k+pb, md=ma*k
    (k1,_,p1,q1,m1),(k2,_,p2,q2,m2)=data[0],data[1]
    qa=q2-q1; qb=q1-qa; pa=p2-p1; pb=p1-pa; ma=m2-m1
    print(f"  FIT: t*=({pa}k+{pb})/({qa}k+{qb}) M={ma}k/({qa}k+{qb}); kappa={ma}k, q-kappa={qa-ma}k+{qb}")
    print(f"  verify k=1..20:", end=" ")
    okv=True
    for k in range(1,21):
        X=X0*k; V=base+[X]; p=pa*k+pb; q=qa*k+qb
        if min(resd(v*p,q) for v in V)!=ma*k: okv=False
    print("OK" if okv else "FAIL")
    # residue table
    print(f"  residue table (v: qq, r=Ak+B; need {ma}k<=r<={qa-ma}k+{qb}):")
    for v in base+['X']:
        vv=(lambda k:X0*k) if v=='X' else (lambda k,v=v:v)
        rs=[];qqs=[]
        for k in [2,3,4]:
            p=pa*k+pb;q=qa*k+qb;N=vv(k)*p; rs.append(N%q);qqs.append(N//q)
        from fractions import Fraction
        A=Fraction(rs[1]-rs[0],1);B=rs[0]-A*2
        qqA=Fraction(qqs[1]-qqs[0],1);qqB=qqs[0]-qqA*2
        lbl=f"{X0}k" if v=='X' else str(v)
        qstr=f"{qqA}k+{qqB}" if qqA!=0 else f"{qqB}"
        print(f"    v={lbl:>5}: qq={qstr:>7} r={A}k+{B}")
