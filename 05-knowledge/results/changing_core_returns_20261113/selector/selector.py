"""Exact physical 881 preparations and finite exponent congruence certificates."""
from pathlib import Path
from functools import lru_cache
from itertools import permutations
from math import gcd,lcm
import argparse,importlib.util,json,sys
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
sys.path.insert(0,str(RESULTS/'berggren_target_return_20261111/reset'))
import funnel_six_state as funnel
spec=importlib.util.spec_from_file_location('old_escape',RESULTS/'eight_nine_twentyseven_20261112/constructors/escape.py')
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)
def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def valid(t):
    n=len(t)
    require(all(0<=x<1<<n and not(x>>i&1) for i,x in enumerate(t)),'diagonal/range')
    require(all(((t[i]>>j&1)+(t[j]>>i&1))==1 for i in range(n) for j in range(i)),'not tournament')
def physical_swap(t,a,b):
    p=list(range(len(t)));p[a],p[b]=p[b],p[a]
    return tuple(sum(1<<p[j] for j in range(len(t)) if t[i]>>j&1) for i in sorted(range(len(t)),key=lambda x:p[x]))
def profiles(t):
    valid(t);full=(1<<len(t))-1
    ends=e.ends(t);starts=e.ends(e.converse(t));h=[sum(x) for x in ends];h[0]=1
    return [(h[full],starts[full][r],h[full^(1<<r)],sum(ends[m][r]*h[full^m] for m in range(1,full+1))) for r in range(len(t))]
def direct_h(t):
    @lru_cache(None)
    def walk(last,left):
        if not left:return 1
        return sum(walk(v,left^(1<<v)) for v in range(len(t)) if left>>v&1 and t[last]>>v&1)
    full=(1<<len(t))-1
    return sum(walk(v,full^(1<<v)) for v in range(len(t)))
def factors(n):
    ans=[];p=2
    while p*p<=n:
        if n%p==0:
            ans.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:ans.append(n)
    return ans
def wheel(z):
    h,s,q,F=z;d=F-q;delta=d-h+s
    if delta%3==0:return dict(delta=delta,d=d,period=1,good=[],reason='delta divisible by three')
    if (q-delta)%2==0:return dict(delta=delta,d=d,period=1,good=[],reason='return parameter even')
    bad=[];period=2
    for p in factors(abs(delta)):
        if p==2:continue
        if d%p==0:
            if q%p==0:return dict(delta=delta,d=d,period=1,good=[],reason='permanent common prime '+str(p))
            continue
        x=(-q*pow(d,-1,p))%p
        vals=[];a=1
        while not vals or a!=1:
            vals.append(a);a=2*a%p
        if x in vals:
            order=len(vals);bad.append(dict(prime=p,modulus=order,residue=vals.index(x)));period=lcm(period,order)
    require(period<=1000000,'wheel too large for explicit residue listing')
    good=[k for k in range(period) if (d*pow(2,k,3)+q-delta)%3==0 and all(k%x['modulus']!=x['residue'] for x in bad)]
    return dict(delta=delta,d=d,period=period,good=good,bad=bad)
def returned(z,k):
    h,s,q,F=z;d=F-q;delta=d-h+s;t=d*(1<<k)+q-delta;b=3*t+2*delta
    return b,t,t,b
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='receipt.json');args=parser.parse_args()
    controls=0
    for n in range(1,5):
        pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
        for mask in range(1<<len(pairs)):
            t=[0]*n
            for bit,(i,j) in enumerate(pairs):
                if mask>>bit&1:t[i]|=1<<j
                else:t[j]|=1<<i
            t=tuple(t);p=profiles(t)
            require(p[0][0]==sum(all(t[a]>>b&1 for a,b in zip(order,order[1:])) for order in permutations(range(n))),'literal H control')
            for z in p:
                wh=wheel(z)
                for k in range(1,101):
                    b,a,_,_=returned(z,k);guard=(gcd(a,b)==1 and a%3==0 and b%3!=0)
                    require(guard==(k%wh['period'] in wh['good']),'wheel vs direct gcd')
                controls+=1
    sources=[]
    for parent_order,rank,expected in [(9,168,(881,169,195,977)),(10,174,(881,145,227,1139))]:
        t,root=funnel.a.birth(funnel.a.unrank(parent_order,rank),0)
        before=profiles(t);require(e.strong(t) and before[root][0]==881,'881 source')
        prepared=physical_swap(t,root,3);after=profiles(prepared)
        require(after[root]==expected and after[root]==before[3],'physical named-root transposition')
        require(e.strong(prepared) and direct_h(prepared)==881,'independent full count')
        flips=[(i,j) for i in range(len(t)) for j in range(i+1,len(t)) if ((t[i]^prepared[i])>>j)&1]
        wh=wheel(expected);literal=[]
        for k in (1,2,3):
            u,r=prepared,root
            for _ in range(k):u,r=e.child(u,r)
            u,r=e.child(u,r,True);u,r=e.child(u,r)
            actual=e.root_state(u,r)
            require(actual==returned(expected,k) and e.strong(u),'actual strong return')
            literal.append(dict(k=k,rows=u,root=r,profile=actual))
        for z in before:
            test=wheel(z)
            for k in range(1,201):
                b,a,_,_=returned(z,k)
                require((gcd(a,b)==1 and a%3==0 and b%3!=0)==(k%test['period'] in test['good']),'881 root wheel')
        sources.append(dict(parent_order=parent_order,parent_rank=rank,rows=t,root=root,root_profiles=before,prepared_rows=prepared,arc_flips=flips,prepared_profile=expected,ready_defect=expected[2]-expected[1],wheel=wh,literal_returns=literal))
    require(sources[0]['wheel']['good']==[1] and sources[0]['wheel']['period']==2,'A branch guard')
    require(sources[1]['wheel']['good']==[0,1,2,3,4,5,6,8,9] and sources[1]['wheel']['period']==10,'B branch guard')
    data=dict(status='PASS',scope='Two strong 881 sources; all labelled rooted controls through order four; all-positive-k congruence proof in PROOF.md',rooted_controls=controls,sources=sources)
    (HERE/args.output).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='PASS',rooted_controls=controls,source_orders=[len(x['rows']) for x in sources],first_returns=[x['literal_returns'][0]['profile'] for x in sources])))
if __name__=='__main__':main()
