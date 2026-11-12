"""First marked jet with factorial and dyadic denominator compilation."""
from backbone import *

def integer_jet(n,B=(3,5)):
    B=validate(n,B);nf=factorial(n);scale=1<<n
    ansU=ansC=0;stats={'backbone_configs':0,'transitions':0}
    fact=[factorial(i) for i in range(n+1)]
    # Insert j cycles into m labelled places: always integer.
    placements={(b,m,j):fact[m+b*j]//(fact[m]*b**j*fact[j])
                for b in range(3,n+1,2) if b not in B
                for m in range(n+1) for j in range((n-m)//b+1)}
    for ms in configs(n,B):
        stats['backbone_configs']+=1
        s=sum(a*m for a,m in zip(B,ms));N=n-s;r=sum(ms)
        v=[1<<(comb(f,2)+f*r) for f in range(N+1)]
        c=[0]*(N+1);mom=[];prev=[]
        for b in range(3,N+1,2):
            if b in B:continue
            R=sum(m*gcd(a,b) for a,m in zip(B,ms))
            ov=[0]*(N+1);oc=[0]*(N+1);om=[[0]*(N+1) for _ in range(len(prev)+1)]
            qs=[scale-(1<<(n-(a-gcd(a,b)))) for a in prev]
            for m in range(N+1):
                if not v[m]:continue
                Q=sum(q*arr[m] for q,arr in zip(qs,mom))
                for j in range((N-m)//b+1):
                    stats['transitions']+=1;t=m+b*j
                    A=placements[b,m,j]<<(j*m+j*(b-1)//2+b*comb(j,2)+j*R)
                    ov[t]+=A*v[m];oc[t]+=A*(c[m]+j*Q)
                    for z,arr in enumerate(mom):om[z][t]+=A*arr[m]
                    om[-1][t]+=A*j*v[m]
            v,c,mom=ov,oc,om;prev.append(b)
        eb=sum(m*(a-1)//2+a*comb(m,2) for a,m in zip(B,ms))
        eb+=sum(ms[i]*ms[j]*gcd(a,b) for i,a in enumerate(B) for j,b in enumerate(B) if i<j)
        zb=1
        for a,m in zip(B,ms):zb*=a**m*fact[m]
        den=zb*fact[N];assert nf%den==0
        factor=(nf//den)<<eb
        ansU+=factor*v[N];ansC+=factor*c[N]
    return F(ansU,nf),F(ansC,nf*scale),stats

if __name__=='__main__':
    rows=[];start=time.time()
    for n in list(range(21))+[32,48,64,80,96]:
        for B in [(),(3,),(3,5),(3,7)] if n<=20 else [(3,5)]:
            t=time.perf_counter();u,c,stats=integer_jet(n,B);secs=time.perf_counter()-t
            if n<=80:
                ur=F(0);cr=F(0);pc=0;t=time.perf_counter()
                for p in parts(n):
                    w,qs=sector(p,B);ur+=w;cr+=w*sum(qs,F(0));pc+=1
                refs=time.perf_counter()-t
                assert (u,c)==(ur,cr)
            else:
                t=time.perf_counter();ur,cr,_=envelope(n,B,True);refs=time.perf_counter()-t
                assert (u,c)==(ur,cr);pc=None
            rows.append({'n':n,'backbone':list(B),'integer_jet_seconds':secs,
                         'reference_seconds':refs,'reference':'partition' if n<=80 else 'rational_jet',
                         'partition_count':pc,'C0':str(u),'C1':str(c),**stats})
    print(json.dumps({'status':'PASS','seconds':time.time()-start,'rows':rows},indent=2))
