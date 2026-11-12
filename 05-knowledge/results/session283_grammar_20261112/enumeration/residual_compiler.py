"""Exact complementary elimination for stratified Burnside defect moments.

Only off-backbone multisets are listed; all backbone completions are summed by
an external-field recursion. No exact tournament count is used by produce().
"""
from fractions import Fraction as F
from collections import Counter,defaultdict
from functools import lru_cache
from math import gcd,comb,factorial
from pathlib import Path
import sys,json,time,argparse
sys.set_int_max_str_digits(0)
OLD=Path(__file__).resolve().parents[2]/'session279_interval_tightening_20261112'/'cross'
sys.path.insert(0,str(OLD))
from ground_moments import certify,produce as old_produce

def residuals(n,allowed,lo=0):
    yield ()
    for i in range(lo,len(allowed)):
        a=allowed[i]
        if a>n:break
        for p in residuals(n-a,allowed,i):yield (a,)+p

def elementary(ds,K):
    e=[F(1)]+[F(0)]*K
    for d in ds:
        q=1-F(1,1<<d)
        for j in range(K,0,-1):e[j]+=q*e[j-1]
    return e

def core(p):
    ds=[];e=sum((a-1)//2 for a in p)
    for i,a in enumerate(p):
        for b in p[i+1:]:
            e+=min(a,b);d=min(a,b)-gcd(a,b)
            if d:ds.append(d)
    z=1
    for a,m in Counter(p).items():z*=a**m*factorial(m)
    return F(1<<e,z),ds

def produce(n,B=(3,5,7,9,11),K=4,projective=False,integer=True):
    if n<0 or K!=4 or len(set(B))!=len(B) or any(b<3 or b%2==0 for b in B):raise ValueError('domain')
    if projective and integer:raise ValueError('signed projective gauge needs rational mode')
    B=tuple(sorted(b for b in B if b<=n));short=(1,)+B
    allowed=tuple(a for a in range(3,n+1,2) if a not in B)
    count=0;groups={};transitions=0
    # Every core has one continuation key; its defect information is separate.
    for p in residuals(n,allowed):
        count+=1;w,ds=core(p);delta=min(ds,default=0)
        cell=(len(ds),delta,sum(d>delta for d in ds))
        N=n-sum(p);q=tuple(sum(gcd(a,b) for a in p) for b in short)
        if projective:
            w*=1<<(N*q[0])
            q=(0,)+tuple((b*q[0]-q[j])//2 for j,b in enumerate(short) if j)
        key=(N,q)
        cs=groups.setdefault(key,defaultdict(lambda:[F(0)]*(K+1)))[cell]
        for j,v in enumerate(elementary(ds,K)):cs[j]+=w*v
    facts=[factorial(i) for i in range(n+1)]
    @lru_cache(None)
    def complete(i,R,q):
        nonlocal transitions
        if i==0:
            value=1<<(comb(R,2)+R*q[0])
            return value if integer else F(value,facts[R])
        b=short[i];ans=0 if integer else F(0)
        for m in range(R//b+1):
            transitions+=1
            exp=m*(b-1)//2+b*comb(m,2)+(m*(R-b*m)-2*m*q[i] if projective else m*q[i])
            if integer:
                denom=b**m*facts[m]*facts[R-b*m]
                factor,rem=divmod(facts[R],denom)
                if rem:raise ArithmeticError('nonintegral labelled cycle-choice coefficient')
                w=factor<<exp
            else:w=(F(1<<exp) if exp>=0 else F(1,1<<(-exp)))/(b**m*facts[m])
            qq=tuple(q[j]+m*((short[j]-gcd(b,short[j]))//2 if projective else gcd(b,short[j])) for j in range(i))
            ans+=w*complete(i-1,R-b*m,qq)
        return ans
    cells=defaultdict(lambda:[F(0)]*(K+1))
    for (R,q),data in groups.items():
        factor=complete(len(short)-1,R,q)
        if integer:factor=F(factor,facts[R])
        for cell,cs in data.items():
            for j,v in enumerate(cs):cells[cell][j]+=factor*v
    packet={'schema':'ground-moments-v1','n':n,'backbone':list(B),'degree':K,
            'producer':'complementary-external-field-v1',
            'cells':[{'M':m,'delta':d,'R':r,'C':list(map(str,cs))} for (m,d,r),cs in sorted(cells.items())]}
    stats={'integer_short':integer,'projective':projective,'residual_cores':count,'continuation_keys':len(groups),'recursion_states':complete.cache_info().currsize,'recursion_cache_hits':complete.cache_info().hits,'transitions':transitions,'moment_cells':len(cells)}
    return packet,stats

def normalized(packet):return packet['cells']
def run(large=(64,70,80)):
    checks=[]
    for n in range(25):
        for B in [(),(3,),(3,5,7),(3,7)]:
            p,s=produce(n,B);raw,_=produce(n,B,projective=False,integer=False);proj,_=produce(n,B,projective=True,integer=False);o,_=old_produce(n,B)
            if normalized(proj)!=normalized(p):raise AssertionError((n,B,'integer/projective'))
            if normalized(raw)!=normalized(p):raise AssertionError((n,B,'raw/projective'))
            if normalized(p)!=normalized(o):raise AssertionError((n,B,'packet disagreement'))
            checks.append([n,list(B)])
    rows=[]
    for n in large:
        t=time.perf_counter();p,s=produce(n);secs=time.perf_counter()-t
        out=certify(p)['ground_with_R']
        row={'n':n,**s,'seconds':round(secs,6),'interval':out}
        if n<=70:
            t=time.perf_counter();o,_=old_produce(n,(3,5,7,9,11))
            row['old_seconds']=round(time.perf_counter()-t,6);row['old_partitions']=o['odd_partitions']
            if normalized(p)!=normalized(o):raise AssertionError((n,'large packet disagreement'))
        rows.append(row)
    return {'status':'PASS','small_packet_controls':checks,'rows':rows}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int);ap.add_argument('--packet');ap.add_argument('--large',default='64,70,80');args=ap.parse_args()
    if args.n is not None:
        t=time.perf_counter();p,s=produce(args.n);result={'n':args.n,**s,'seconds':round(time.perf_counter()-t,6),'certificate':certify(p)}
        if args.packet:Path(args.packet).write_text(json.dumps(p,indent=2)+'\n')
    else:result=run(tuple(map(int,args.large.split(','))))
    print(json.dumps(result,indent=2))
