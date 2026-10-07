"""For each class r mod 2^K (K small) find: descent depth; translation-join (n-d, any d) depth; minimal backward depth i of a
branch certificate (ratio<1) at each s.  Report classes certified ONLY by branch words of backward depth >= 2 (i.e. not by
descent, joins, or the depth-1 predecessor (2x_s-1)/3 of Angeltveit 2026), with an explicit certificate m(n)."""
import itertools, random, sys
from fractions import Fraction
def T(x): return x//2 if x%2==0 else (3*x+1)//2
random.seed(11)
K=int(sys.argv[1])
def analyze(r):
    reps=[r+(random.getrandbits(90)<<K) for _ in range(3)]
    orbs=[]
    for n in reps:
        xs=[n];a=[0]
        for s in range(K): a.append(a[-1]+(xs[-1]&1)); xs.append(T(xs[-1]))
        orbs.append((n,xs))
    a_=a
    res={'desc':None,'join':None,'pred1':None,'branch':None}
    for s in range(1,K+1):
        if res['desc'] is None and all(xs[s]<n for n,xs in orbs): res['desc']=s
        # joins: n-d merging at equal time s, uniform: test d up to 2^(s - a_s)
        if res['join'] is None:
            for d in range(1, 2**max(0,s-a_[s])+1):
                if all(T_iter(n-d,s)==xs[s] for n,xs in orbs): res['join']=(s,d);break
        # branch words
        for i in range(0,s):
            found=None
            for u in itertools.product((0,1),repeat=i):
                if sum(u)>a_[s]: continue
                ok=True; vals=[]
                for n,xs in orbs:
                    v=xs[s]
                    for bit in u:
                        if bit:
                            if v%3!=2: ok=False;break
                            v=(2*v-1)//3
                        else: v=2*v
                    if not ok or not v<n: ok=False;break
                    vals.append((n,v))
                if ok: found=(s,i,u,vals);break
            if found:
                if res['branch'] is None or i<res['branch'][1]: res['branch']=found
                if i==1 and res['pred1'] is None: res['pred1']=s
                break
    return res
def T_iter(x,s):
    for _ in range(s): x=T(x)
    return x
newonly=[]
cnt={'desc':0,'join':0,'pred1':0,'deep':0,'none':0}
for r in range(2**K):
    res=analyze(r)
    if res['desc']: cnt['desc']+=1; continue
    if res['join']: cnt['join']+=1
    if res['pred1']: cnt['pred1']+=1
    if not res['join'] and not res['pred1']:
        if res['branch']: cnt['deep']+=1; newonly.append((r,res['branch']))
        else: cnt['none']+=1
print("K",K,"classes",2**K,cnt)
for r,(s,i,u,vals) in newonly[:12]:
    (n,v)=vals[0]; (n2,v2)=vals[1]
    rho=Fraction(v2-v, n2-n); kap=v-rho*n
    print(f"  r={r} mod 2^{K}: from x_{s}, backward word (reading backward, 1=odd) {''.join(map(str,u))} (depth {i}), m = {rho}*n + {kap}")
