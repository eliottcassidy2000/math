"""First-jet response compressed into primitive-word divisor moments."""
from integer_jet import *

def divisor_jet(n,B=(3,5)):
    B=validate(n,B);nf=factorial(n);scale=1<<n
    fact=[factorial(i) for i in range(n+1)]
    ds=list(range(1,n//3+1,2));primitive={}
    for d in ds:primitive[d]=(1<<d)-sum(primitive[e] for e in ds if e<d and d%e==0)
    lengths=[b for b in range(3,n+1,2) if b not in B]
    div={b:[i for i,d in enumerate(ds) if b%d==0] for b in lengths}
    placements={(b,m,j):fact[m+b*j]//(fact[m]*b**j*fact[j]) for b in lengths
                for m in range(n+1) for j in range((n-m)//b+1)}
    outU=outC=0;stats={'backbone_configs':0,'transitions':0,'divisor_moments':len(ds)}
    for ms in configs(n,B):
        stats['backbone_configs']+=1
        N=n-sum(a*m for a,m in zip(B,ms));r=sum(ms)
        v=[1<<(comb(f,2)+f*r) for f in range(N+1)]
        c=[0]*(N+1);total=[0]*(N+1);res=[[0]*(N+1) for _ in ds]
        for b in lengths:
            if b>N:break
            R=sum(m*gcd(a,b) for a,m in zip(B,ms))
            ov=[0]*(N+1);oc=[0]*(N+1);ot=[0]*(N+1);orr=[[0]*(N+1) for _ in ds]
            for m in range(N+1):
                if not v[m]:continue
                Q=scale*total[m]-sum(primitive[ds[i]]*res[i][m] for i in div[b])
                for j in range((N-m)//b+1):
                    stats['transitions']+=1;t=m+b*j
                    A=placements[b,m,j]<<(j*m+j*(b-1)//2+b*comb(j,2)+j*R)
                    ov[t]+=A*v[m];oc[t]+=A*(c[m]+j*Q);ot[t]+=A*(total[m]+j*v[m])
                    for i,arr in enumerate(res):orr[i][t]+=A*arr[m]
                    if j:
                        val=(A*j*v[m])<<(n-b)
                        for i in div[b]:orr[i][t]+=val
            v,c,total,res=ov,oc,ot,orr
        eb=sum(m*(a-1)//2+a*comb(m,2) for a,m in zip(B,ms))
        eb+=sum(ms[i]*ms[j]*gcd(a,b) for i,a in enumerate(B) for j,b in enumerate(B) if i<j)
        zb=1
        for a,m in zip(B,ms):zb*=a**m*fact[m]
        den=zb*fact[N];assert nf%den==0
        factor=(nf//den)<<eb;outU+=factor*v[N];outC+=factor*c[N]
    return F(outU,nf),F(outC,nf*scale),stats

if __name__=='__main__':
    start=time.time();rows=[]
    for n in list(range(21))+[32,48,64,80,96,128]:
        for B in [(),(3,5),(3,7)] if n<=20 else [(3,5)]:
            t=time.perf_counter();u,c,st=divisor_jet(n,B);td=time.perf_counter()-t
            t=time.perf_counter();ur,cr,_=integer_jet(n,B);ti=time.perf_counter()-t
            assert (u,c)==(ur,cr)
            rows.append({'n':n,'backbone':list(B),'divisor_seconds':td,'integer_length_seconds':ti,
                         'C0':str(u),'C1':str(c),**st})
    print(json.dumps({'status':'PASS','seconds':time.time()-start,'rows':rows},indent=2))
