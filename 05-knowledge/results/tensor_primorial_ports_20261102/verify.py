"""Exact audit of the order42 certificate and its scalar stopping witness."""
from fractions import Fraction as F
from math import factorial
from itertools import combinations, product
from pathlib import Path
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def run():
    T=sum(F(5)**j/factorial(j) for j in range(13))
    c0=F(2483,10000); c=F(99,100); boot=[]; endpoints=[]
    need((1-1/T)/4>c0,'initial strict reserve')
    for R in (21,22,23):
        K=F(4,R)*(1-F(3,4)**R); E=F(1,4)**R; N=2*R+1
        need((1-1/T)/5>K and K<F(1,5),'V>5/4 and local lemma domain')
        mass=lambda z:z/K-N*K/8
        steps=[]; z=c0
        for i in range(30):
            m=mass(z); raw=1-(1-m/R)**R-E
            nxt=F(raw.numerator*10**6//raw.denominator,10**6)
            need(0<m<F(R,4) and nxt>z,'strict bootstrap progress')
            steps.append({'from':str(z),'mass':str(m),'to':str(nxt),'floor_slack':str(raw-nxt)})
            z=nxt
            if z>=c:break
        need(z>=c,'bootstrap terminates at required reserve')
        boot.append({'R':R,'steps':steps})
        m=mass(c)
        for n in (2*R,2*R+1):
            eps=int(n%2==0); k=K-F(4,R)*eps*E
            phi=lambda x:1-(1-x/R)**R-eps*E-k*k*R-k*(1-4*k)*x
            need(phi(F(R,4))==0,'capacity endpoint identity')
            need(phi(m)>F(1,25),'positive lower endpoint')
            endpoints.append({'n':n,'mass':str(m),'slack':str(phi(m))})
    # Full-domain algebra is in the proof; this grid catches implementation/sign errors.
    G=lambda w,s:max(w*s,w*s/2+s*s/4)
    local=0
    for k in [F(1,20),F(1,10),F(19,100),F(1,5)]:
        for i in range(101):
            w=F(i,400)
            for j in range(101):
                s=F(j,50)
                need(G(w,s)>=k*s-k*k*(1-4*w),'capacity potential')
                local+=1
    # Literal small graphs, independent of the scalar minimization.
    edges=list(combinations(range(4),2)); graphs=0
    for values in product([F(0),F(1,8),F(1,4)],repeat=6):
        q=dict(zip(edges,values)); todo={e for e in edges if q[e]}; groups=[]
        while todo:
            e=max(todo,key=lambda e:(q[e],e)); own={f for f in todo if set(e)&set(f)}
            todo-=own; groups.append((q[e],sum(q[f] for f in own-{e})))
        W=sum(q[e]*q[f] for e,f in combinations(edges,2) if set(e)&set(f))
        M=sum(w for w,s in groups);D=sum(s for w,s in groups);k=F(1,5)
        need(W>=k*D-k*k*(2-4*M),'literal padded-group bound')
        graphs+=1
    R=20;n=41;k=F(4,R)*(1-F(3,4)**R);x=F(1,70)
    M=R*x;D=F(2);V=M+D;B=F(8,41)
    groups=[F(11,105)]*18+[F(3,35),F(1,35)]
    need(sum(groups)==D,'hostile mass')
    need(all(s<=2*(n-2*i-2)*x for i,s in enumerate(groups)),'hostile group capacity')
    need(B>=sum(G(x,s) for s in groups),'hostile disjoint invoices')
    need(B>=2*D*D/n and B>=2*V*V/n-V/4,'hostile degree bounds')
    need(B<k*D and B+F(1,4)<k*V,'hostile linear and exponential reserves')
    slack=k*V-B-(1-(1-x)**R)
    need(slack>0,'hostile tensor reserve')
    improved=2*V*V/n-x*V
    need(improved>B and improved+1-(1-x)**R>k*V,'max-weight coordinate rejects hostile')
    m=c0/k-n*k/8; nxt=1-(1-m/R)**R-F(1,4)**R
    need(nxt<c0,'order41 bootstrap fails')
    return {'status':'PASS','new_orders':list(range(42,48)),'remaining_open':[10,41],
        'bootstrap':boot,'endpoints':endpoints,'local_grid_checks':local,'literal_graphs':graphs,
        'order41_relaxation_hostile':{'w':str(x),'M':str(M),'D':str(D),'V':str(V),'B':str(B),
            'tensor_slack':str(slack),'new_max_weight_B':str(improved),
            'meaning':'Feasible only for listed scalar inequalities; not an actual tournament. Keeping max q rejects it.'}}

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
