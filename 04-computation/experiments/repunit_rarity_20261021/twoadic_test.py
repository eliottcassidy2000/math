# Test: does the local log-residual of P_d(1) (strong classes, all orders) depend on v2(d+1) (closeness to -1) and
# v2(d-1) (closeness to +1)?  Also per fixed order j (sigma(j,d)), and a permutation test for significance.
import json, math, random
from collections import defaultdict
tab={int(k):{int(a):b for a,b in v.items()} for k,v in json.load(open('strong_table_lt375.json')).items()}
def v2(x):
    c=0
    while x%2==0: x//=2; c+=1
    return c
def residuals(count, ds, frac=0.15):
    res={}
    for d in ds:
        w=[count[e] for e in ds if abs(e-d)<=frac*d and e!=d]
        if len(w)<4: continue
        res[d]=math.log((count[d]+0.5)/(sum(w)/len(w)+0.5))
    return res
def report(name,count,ds):
    res=residuals(count,ds)
    for key,f in [('v2(d+1)',lambda d:v2(d+1)),('v2(d-1)',lambda d:v2(d-1))]:
        g=defaultdict(list)
        for d,r in res.items(): g[min(f(d),5)].append(r)
        s=' '.join(f'{k}:{sum(v)/len(v):+.2f}(n={len(v)})' for k,v in sorted(g.items()))
        # permutation test of the correlation between f(d) and residual
        xs=[min(f(d),5) for d in res]; ys=list(res.values())
        def corr(a,b):
            ma=sum(a)/len(a); mb=sum(b)/len(b)
            num=sum((x-ma)*(y-mb) for x,y in zip(a,b)); den=math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
            return num/den
        c0=corr(xs,ys); random.seed(1); cnt=0; R=4000
        for _ in range(R):
            sh=ys[:]; random.shuffle(sh)
            if abs(corr(xs,sh))>=abs(c0): cnt+=1
        print(f'{name:22s} by {key}: {s}   corr={c0:+.3f} perm-p={cnt/R:.4f}')
ds=[d for d in range(25,375,2)]
report('ALL orders P_d(1)',{d:sum(tab[d].values()) for d in ds},ds)
for j in [7,8,9]:
    dsj=[d for d in range(1,375,2) if tab[d].get(j,0)>0 or any(tab[e].get(j,0)>0 for e in range(max(1,d-30),d+31,2))]
    dsj=[d for d in dsj if d>=25]
    report(f'order {j} sigma(j,d)',{d:tab[d].get(j,0) for d in dsj},dsj)
