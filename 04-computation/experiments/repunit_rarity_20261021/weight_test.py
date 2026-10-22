# Compare predictors of the rarity residual of P_d(1) (strong classes, all orders, odd 25 <= d < 375):
# s2(Q) (binary weight of Q=(d-1)/2), s2 of the low 5 bits of Q, v2(d+1), v2(d-1); partial check controlling for v2(d+1).
import json, math, random
tab={int(k):{int(a):b for a,b in v.items()} for k,v in json.load(open('strong_table_lt375.json')).items()}
tot={d:sum(v.values()) for d,v in tab.items()}
ds=[d for d in range(25,375,2)]
res={}
for d in ds:
    w=[tot[e] for e in ds if abs(e-d)<=0.15*d and e!=d]
    res[d]=math.log((tot[d]+0.5)/(sum(w)/len(w)+0.5))
def v2(x):
    c=0
    while x%2==0: x//=2; c+=1
    return c
s2=lambda x: bin(x).count('1')
feats={'s2(Q)':lambda d:s2((d-1)//2), 's2(Q mod 32)':lambda d:s2(((d-1)//2)%32), 'v2(d+1)':lambda d:min(v2(d+1),6),
       'v2(d-1)':lambda d:min(v2(d-1),6), 's2(d)':lambda d:s2(d), 'trailing ones of Q':lambda d:min(v2((d-1)//2+1),6)}
def corr(a,b):
    ma=sum(a)/len(a); mb=sum(b)/len(b)
    num=sum((x-ma)*(y-mb) for x,y in zip(a,b)); den=math.sqrt(sum((x-ma)**2 for x in a)*sum((y-mb)**2 for y in b))
    return num/den
ys=[res[d] for d in ds]
random.seed(7)
for name,f in feats.items():
    xs=[f(d) for d in ds]; c0=corr(xs,ys); cnt=0; R=5000
    for _ in range(R):
        sh=ys[:]; random.shuffle(sh)
        if abs(corr(xs,sh))>=abs(c0): cnt+=1
    print(f'{name:20s} corr={c0:+.3f}  perm-p={cnt/R:.4f}')
# partial: residualize on v2(d+1) group means, then correlate with s2(Q)
from collections import defaultdict
g=defaultdict(list)
for d in ds: g[min(v2(d+1),6)].append(res[d])
gm={k:sum(v)/len(v) for k,v in g.items()}
rr=[res[d]-gm[min(v2(d+1),6)] for d in ds]
xs=[s2((d-1)//2) for d in ds]; c0=corr(xs,rr); cnt=0
for _ in range(5000):
    sh=rr[:]; random.shuffle(sh)
    if abs(corr(xs,sh))>=abs(c0): cnt+=1
print(f's2(Q) after removing v2(d+1) group means: corr={c0:+.3f} perm-p={cnt/5000:.4f}')
