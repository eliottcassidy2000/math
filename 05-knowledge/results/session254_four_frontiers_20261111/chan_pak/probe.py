"""Exact small controls for cost-zero padding normal forms; stdlib only."""
import math, json, random, hashlib, subprocess

def profile(n,mask):
    out=[0]*n; z=0
    for i in range(n):
        for j in range(i+1,n):
            a,b=(i,j) if (mask>>z)&1 else (j,i)
            out[a]|=1<<b; z+=1
    size=1<<n; end=[[0]*n for _ in range(size)]; h=[0]*size
    for s in range(1,size):
        for v in range(n):
            if s>>v&1:
                t=s^(1<<v)
                end[s][v]=1 if not t else sum(end[t][u] for u in range(n) if t>>u&1 and out[u]>>v&1)
        h[s]=sum(end[s])
    p=[[0]*(n+1) for _ in range(size)]; p[0][0]=1
    for s in range(1,size):
        bit=s&-s; t=s
        while t:
            if t&bit:
                for k in range(1,n+1): p[s][k]+=h[t]*p[s^t][k-1]
            t=(t-1)&s
    c=[0]*(n+1); fall=[1]
    for k in range(1,n+1):
        f=[0]*(len(fall)+1)
        for j,a in enumerate(fall): f[j]-=(k-1)*a; f[j+1]+=a
        fall=f
        for j,a in enumerate(fall): c[j]+=p[-1][k]*a
    return c,out

def data(c):
    H=sum(c); D=sum(i*x for i,x in enumerate(c)); E=sum(math.comb(i,2)*x for i,x in enumerate(c))
    return H,D,E,(2*H*E-D*D+H*D)%3

def val(c,x): return sum(a*x**i for i,a in enumerate(c))
def jet(c,m,k): return [sum((-1)**(i-j)*math.comb(i,j)*val(c,j) for j in range(i+1))%m for i in range(1,k)]
def require(t):
    if not t: raise RuntimeError('check failed')

def main():
    rng=random.Random(254); cases=[(n,m) for n in range(1,6) for m in range(1<<(n*(n-1)//2))]
    cases += [(6,rng.randrange(1<<15)) for _ in range(500)]
    tests=0; fibres={}; hostile=None; states={}
    for n,m in cases:
        c,rows=profile(n,m); H,D,E,K=data(c)
        if H%3==0: continue
        fibres.setdefault(H,{}) .setdefault(K,{'n':n,'mask':m,'rows':rows,'Q':c,'H':H,'K':K})
        for parity in (0,1):
            for mod,lo,period,kap in ((9,2,6,6),(27,3,18,9)):
                r=next(r for r in range(lo,lo+period) if (n+r)%2==parity and (D+r*H)%(mod//3)==0)
                q=[0]*r+c; hh,dd,ee,kk=data(q)
                require(hh==H and kk==K and dd%(mod//3)==0)
                expected=[]
                for x in range(mod):
                    if x%3==0: want=0
                    else:
                        sign=1 if x%3==1 else -1
                        u=((sign*x)-1)//3
                        want=(H+(9*u*u*ee if mod==27 else 0))*(sign**parity)
                    require(val(q,x)%mod==want%mod)
                key=(mod,H%mod,parity,K if mod==27 else None)
                j=jet(q,mod,kap)
                if key in states: require(states[key]==j)
                states[key]=j; tests+=1
    for H,ks in sorted(fibres.items()):
        if len(ks)>1: hostile=list(ks.values()); break
    a,arows=profile(6,11356); b,brows=profile(6,11356^(1<<13))
    require(data(a)[0]==data(b)[0]==37 and data(a)[3]==2 and data(b)[3]==0)
    require([i for i in range(6) if arows[i]!=brows[i]]==[3,5])
    surgery={'edge':[3,5], 'before_rows':arows,'after_rows':brows,'before_Q':a,'after_Q':b,'H':37,'before_K':2,'after_K':0}
    for c0 in ([0,1],[0,0,4,0,1],[0,0,24,0,12,0,1]):
        def mul(c,d):
            ans=[0]*(len(c)+len(d)-1)
            for i,x in enumerate(c):
                for j,y in enumerate(d): ans[i+j]+=x*y
            return ans
        aa=mul(a,c0);bb=mul(b,c0)
        require(data(aa)[0]==data(bb)[0] and (data(bb)[3]-data(aa)[3])%3==1)
    print(json.dumps({'neutral_edge_surgery':surgery,'universe':'All labelled tournaments through order5 plus500 deterministic order6 samples','cases':len(cases),'normalizations_checked':tests,'distinct_normal_form_states':len(states),'same_exact_H_distinct_curvature':hostile,'curvatures_seen':sorted({k for ks in fibres.values() for k in ks})},indent=2,sort_keys=True))
if __name__=='__main__': main()
