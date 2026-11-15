import os
HERE=os.path.dirname(os.path.abspath(__file__))
RIB=os.path.join(HERE,'..','..','..','..','04-computation','experiments','ribbon_spectrum_20261104')
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'w29_lp_reachable_sectors.py')).read().rsplit("def feasible(rho):",1)[0].replace("print('invariance violations (numeric):',viol)",""))
from fractions import Fraction as Fr
from sympy import Rational
def fr_ray(a):
    # rational ray close to angle a
    return (Fr(math.cos(a)).limit_denominator(10**7), Fr(math.sin(a)).limit_denominator(10**7))
def cross(u,v): return u[0]*v[1]-u[1]*v[0]
def apply(X,v): return (Fr(int(X[0,0].p),int(X[0,0].q))*v[0]+Fr(int(X[0,1].p),int(X[0,1].q))*v[1], Fr(int(X[1,0].p),int(X[1,0].q))*v[0]+Fr(int(X[1,1].p),int(X[1,1].q))*v[1])
def inside(v,lo,hi):  # v within closed sector from lo to hi (counterclockwise, sector < pi)
    return cross(lo,v)>=0 and cross(v,hi)>=0
Kf={s:[K[s][0],K[s][1]] for s in rec}
for it in range(400):
    changed=False
    for (s0,j,t),X in T.items():
        Xf=np.array(X.tolist(),dtype=float)
        for a in Kf[s0]:
            bb=ang(Xf@ray(a))
            if bb<Kf[t][0]-1e-15: Kf[t][0]=bb; changed=True
            if bb>Kf[t][1]+1e-15: Kf[t][1]=bb; changed=True
    if not changed: break
print('hull iterations',it,{s:(round(Kf[s][0],9),round(Kf[s][1],9)) for s in rec})
K={s:tuple(Kf[s]) for s in rec}
for eps in [1e-9,1e-8,1e-7,1e-6,1e-5,1e-4]:
    def snap(a,sgn):
        for p in range(-6,7):
            for q in range(-6,7):
                if (p,q)!=(0,0) and abs(math.atan2(q,p)-a)<1e-8: return (Fr(p),Fr(q))
        return fr_ray(a+sgn*eps)
    KR={s:(snap(K[s][0],-1),snap(K[s][1],1)) for s in rec}
    bad=0
    for (s,j,t),X in T.items():
        for v in KR[s]:
            w=apply(X,v)
            if not inside(w,*KR[t]): bad+=1
    print('eps',eps,'exact invariance failures',bad)
    if bad==0: break
# LP at rho0 with these rational rays, then exact check
rho0=Fr(17699,10000)
idx={s:i for i,s in enumerate(rec)}; nv=2*len(rec)
A=[];b=[]
for (s,j,t),X in T.items():
    Xf=np.array(X.tolist(),dtype=float)
    for v in KR[s]:
        e=np.array([float(v[0]),float(v[1])]); row=np.zeros(nv)
        row[2*idx[t]:2*idx[t]+2]-=Xf@e; row[2*idx[s]:2*idx[s]+2]+=float(rho0)**Lw[j]*e
        A.append(row); b.append(-1e-9)
for s in rec:
    for v in KR[s]:
        e=np.array([float(v[0]),float(v[1])]); row=np.zeros(nv); row[2*idx[s]:2*idx[s]+2]=-e; A.append(row); b.append(-1.0)
res=linprog(np.zeros(nv),A_ub=np.array(A),b_ub=np.array(b),bounds=[(None,None)]*nv,method='highs')
print('LP status',res.status)
phi={s:(Fr(res.x[2*idx[s]]).limit_denominator(10**9),Fr(res.x[2*idx[s]+1]).limit_denominator(10**9)) for s in rec}
ok=True; worst=None
for (s,j,t),X in T.items():
    for v in KR[s]:
        w=apply(X,v)
        lhs=phi[t][0]*w[0]+phi[t][1]*w[1]; rhs=rho0**Lw[j]*(phi[s][0]*v[0]+phi[s][1]*v[1])
        if lhs<rhs: ok=False
        r=lhs/rhs if rhs>0 else None
        worst=r if worst is None or (r is not None and r<worst) else worst
for s in rec:
    for v in KR[s]:
        if phi[s][0]*v[0]+phi[s][1]*v[1]<=0: ok=False
print('exact certificate at rho0 = 1.7699:',ok,' min slack ratio',float(worst))
# H positivity: fin.B_s x > 0 on sector rays
def F(x): x=Rational(x); return Fr(int(x.p),int(x.q))
for s0 in rec:
    Bm=Bs[s0]
    for v in KR[s0]:
        Hv=sum(F(fin[0,k])*(F(Bm[k,0])*v[0]+F(Bm[k,1])*v[1]) for k in range(Bm.shape[0]))
        print('  state',s0,'H on sector ray >0:',Hv>0)
print('done')
