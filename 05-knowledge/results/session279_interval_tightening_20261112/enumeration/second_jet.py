"""Degree-two marked transfer: C0,C1,C2 without unrestricted partitions."""
from backbone import *

def second_jet(n,B=(3,5)):
    B=validate(n,B);ans=[F(0)]*3;stats={'backbone_configs':0,'transitions':0}
    for ms in configs(n,B):
        stats['backbone_configs']+=1
        N=n-sum(a*m for a,m in zip(B,ms));r=sum(ms)
        v=[F(1<<(comb(f,2)+f*r),factorial(f)) for f in range(N+1)]
        cs=[[F(0)]*(N+1) for _ in range(2)]
        mom=[];mixed=[];second={};prev=[]
        for b in range(3,N+1,2):
            if b in B:continue
            tlen=len(prev);R=sum(m*gcd(a,b) for a,m in zip(B,ms))
            zero=lambda:[F(0)]*(N+1)
            ov=zero();oc=[zero(),zero()]
            om=[zero() for _ in range(tlen+1)];on=[zero() for _ in range(tlen+1)]
            os={(a,c):zero() for a in range(tlen+1) for c in range(a,tlen+1)}
            qs=[1-F(1,1<<(min(a,b)-gcd(a,b))) for a in prev]
            for m in range(N+1):
                if not v[m]:continue
                Q=sum((q*arr[m] for q,arr in zip(qs,mom)),F(0))
                Q2=sum((q*q*arr[m] for q,arr in zip(qs,mom)),F(0))
                QQ=sum((qs[a]*qs[c]*arr[m]*(1 if a==c else 2) for (a,c),arr in second.items()),F(0))
                NQ=sum((q*arr[m] for q,arr in zip(qs,mixed)),F(0))
                MQ=[sum((qs[a]*second[tuple(sorted((a,c)))][m] for a in range(tlen)),F(0)) for c in range(tlen)]
                for j in range((N-m)//b+1):
                    stats['transitions']+=1;ix=m+b*j
                    w=F(1<<(j*m+j*(b-1)//2+b*comb(j,2)+j*R),b**j*factorial(j))
                    ov[ix]+=w*v[m];oc[0][ix]+=w*(cs[0][m]+j*Q)
                    oc[1][ix]+=w*(cs[1][m]+j*NQ+(j*j*QQ-j*Q2)/2)
                    for a in range(tlen):
                        om[a][ix]+=w*mom[a][m];on[a][ix]+=w*(mixed[a][m]+j*MQ[a])
                        os[a,tlen][ix]+=w*j*mom[a][m]
                    om[-1][ix]+=w*j*v[m];on[-1][ix]+=w*(j*cs[0][m]+j*j*Q)
                    for ac,arr in second.items():os[ac][ix]+=w*arr[m]
                    os[tlen,tlen][ix]+=w*j*j*v[m]
            v,cs,mom,mixed,second=ov,oc,om,on,os;prev.append(b)
        w=bweight(B,ms)
        for i,a in enumerate([v[N],cs[0][N],cs[1][N]]):ans[i]+=w*a
    return ans,stats

if __name__=='__main__':
    rows=[];start=time.time()
    for n in list(range(19))+[24,32,40]:
        for B in [(),(3,5)] if n<=18 else [(3,5)]:
            t=time.perf_counter();cs,stats=second_jet(n,B);seconds=time.perf_counter()-t
            reference=[F(0)]*3
            for p in parts(n):
                w,qs=sector(p,B)
                for i,a in enumerate(elementary(qs,2)):reference[i]+=w*a
            assert cs==reference
            rows.append({'n':n,'backbone':list(B),'seconds':seconds,'C':[str(v) for v in cs],**stats})
    print(json.dumps({'status':'PASS','total_seconds':time.time()-start,'rows':rows},indent=2))
