#!/usr/bin/env python3
from automata import *
from itertools import product
import json,random

def require(x,detail):
    if not x:raise RuntimeError(detail)

def core(k,bits):
    a=[0]*k;q=0
    for i in range(k):
        for j in range(i+1,k):
            if bits>>q&1:a[j]|=1<<i
            else:a[i]|=1<<j
            q+=1
    return tuple(a)

def joint_dp(out,k):
    n=len(out);ans=[]
    for start in range(k):
        dp=[[0]*n for _ in range(1<<n)];dp[1<<start][start]=1
        for mask in range(1,1<<n):
            for v,c in enumerate(dp[mask]):
                if not c:continue
                rest=out[v]&~mask
                while rest:
                    bit=rest&-rest;rest-=bit;u=bit.bit_length()-1
                    dp[mask|bit][u]+=c
        ans.append(tuple(dp[-1][:k]))
    return tuple(ans)

report={'status':'PASS','literal_cases':0,'random_cases':0,'boundary_cases':{}}

def check(apex,w):
    tf=Transfer(apex);out=tournament(apex,w)
    h,starts,ends=subset_dp(out);v,(st,en)=tf.evaluate(w)
    require((v,st,en)==(h,starts[:len(apex)],ends[:len(apex)]),('transfer',apex,w,h,v,st,en))
    nw,deleted=normalize(apex,w)
    actual=tuple(sorted(c[0]-len(apex) for c in strong_components(out) if len(c)==1 and c[0]>=len(apex)))
    require(deleted==actual,('delete',apex,w,deleted,actual))
    require(subset_dp(tournament(apex,nw))[0]==h,('H normalization',apex,w))
    require(normalize(apex,nw)==(nw,()),('idempotence',apex,w))
    require(len(nw)<=length_cap(h),('Busch cap',apex,w,h,nw))
    ew,_=normalize_endpoints(apex,w)
    eh,es,ee=subset_dp(tournament(apex,ew))
    require((eh,es[:len(apex)],ee[:len(apex)])==(h,starts[:len(apex)],ends[:len(apex)]),('endpoint normalization',apex,w,ew))
    require(len(ew)<=length_cap(h)+2,('endpoint cap',apex,w,ew))
    state=tf.initial
    for mask in w:state=tf.step(state,mask)
    joint=joint_dp(out,len(apex))
    require(tf.joint_endpoints(state)==joint,('joint transfer',apex,w))
    require(joint_dp(tournament(apex,ew),len(apex))==joint,('joint normalization',apex,w,ew))

# Full stated labeled universe, including nontransitive three-apex cores.
for k,maxlen in [(0,4),(1,7),(2,5),(3,3)]:
    for bits in range(1<<(k*(k-1)//2)):
        apex=core(k,bits)
        for n in range(maxlen+1):
            for w in product(range(1<<k),repeat=n):check(apex,w);report['literal_cases']+=1
rng=random.Random(221)
for k in (2,3,4):
    for _ in range(30):
        apex=core(k,rng.randrange(1<<(k*(k-1)//2)))
        w=tuple(rng.randrange(1<<k) for _ in range(rng.randrange(4,8)))
        check(apex,w);report['random_cases']+=1
# Normalization is NOT a continuation equivalence; this test must stay.
a=(2,0)
require(normalize(a,(3,))[0]==(), 'bad hostile normalization')
left=Transfer(a).evaluate((3,0))[0];right=Transfer(a).evaluate((0,))[0]
require((left,right)==(5,1),('continuation hostile',left,right))
report['boundary_cases']['continuation_hostile']={'word':[3],'suffix':[0],'original_H':left,'normalized_H':right}
# Saturated numerical bound and ordinary scalar edge cases.
for m in range(1,1000):
    z=length_cap(m);require(5**z<=m**3,('coarse integer bound',m,z))
report['boundary_cases']['caps']={str(m):length_cap(m) for m in [1,3,5,7,9,15,21,25,63,1000000]}
# Exact all-length inverse census; independently check every returned graph.
r=census(a,63)
for h,w in r['witnesses'].items():require(subset_dp(tournament(a,w))[0]==int(h),('census witness',h,w))
missing=[m for m in range(1,64,2) if str(m) not in r['witnesses']]
require(missing==[7,21],('exact all-length census',missing))
r['missing_odd']=missing;report['two_apex_census']=r
# Sharpness of the core-independent constant: ordinal C3/S4 blocks, one apex each.
for sizes in [(2,),(3,),(2,3),(3,3),(2,2,3)]:
    k=len(sizes);apex=core(k,0);w=[];expected=1
    for i,n in enumerate(sizes):
        # One-apex C3 word10 or S4 word100 has H=3 or5.
        expected*=3 if n==2 else 5
        for j in range(n):w.append((1<<i)-1 | ((1<<i) if j==0 else 0))
    check(apex,tuple(w));h=Transfer(apex).evaluate(w)[0]
    require(h==expected,('sharp product',sizes,h,expected))
    require(len(w)==sum(sizes),('cost',sizes))
    report['boundary_cases']['product_'+','.join(map(str,sizes))]={'H':h,'chain':len(w),'cap':length_cap(h),'word':w}
report['target_compiler']=[]
for target in [1,7,21,23,63]:
    r=compile_target((2,0),target)
    require((r['status']=='FOUND')==(target not in [7,21]),('target status',r))
    if r['status']=='FOUND':require(subset_dp(tuple(r['outgoing']))[0]==target,('target graph',r))
    report['target_compiler'].append(r)
report['endpoint_compiler']=[]
for w in [(),(0,),(3,),(0,3),(3,0),(1,2,2)]:
    a=(2,0);tf=Transfer(a);state=tf.initial
    for mask in w:state=tf.step(state,mask)
    st,en=tf.endpoints(state)
    wanted={'start':st,'end':en,'joint':tf.joint_endpoints(state)}
    r=compile_target(a,tf.output(state),wanted)
    require(r['status']=='FOUND',('endpoint compiler status',w,r))
    out=tuple(r['outgoing']);h,st2,en2=subset_dp(out)
    require(h==r['target'] and st2[:2]==st and en2[:2]==en and joint_dp(out,2)==wanted['joint'],('endpoint compiler graph',w,r))
    report['endpoint_compiler'].append(r)
for r in report['endpoint_compiler']:
    wanted=json.loads(json.dumps(r['endpoint_constraints']))
    reread=compile_target((2,0),r['target'],wanted)
    require(reread==r,('JSON endpoint coercion',r,reread))
require(compile_target((2,0),1,{'start':[2,0]})['status']=='ABSENT','oversized endpoint')
require(compile_target((2,0),1,{'end':[-1,0]})['status']=='ABSENT','negative endpoint')
report['endpoint_input_controls']=8
print(json.dumps(report,indent=2,sort_keys=True))
