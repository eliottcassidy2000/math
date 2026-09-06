"""Full thirteen-row decoder-span equality classification in the Q^2 box.

The target is W=V_dec, not safety. Three or more graph components do not
determine rank W; a literal rank-eleven hostile retains this distinction.
"""
from itertools import combinations
from functools import reduce
from math import gcd,isqrt
import hashlib,json,sys

sys.stdout.reconfigure(newline='\n')
Q=91**6
GATES=0

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(label)

def cgcd(v):return reduce(gcd,v,0)
def ceildiv(a,b):return -((-a)//b)

def admissible_sums():
    out=set()
    for n in range(2,357):
        residual=n;good=True
        for p in range(2,isqrt(n)+1):
            if residual%p:continue
            e=0
            while residual%p==0:residual//=p;e+=1
            if p%3!=2 or e>2:good=False;break
        if good and (residual==1 or residual%3==2):out.add(n)
    return out

SUMS=admissible_sums()

def graph(row):
    n=len(row);neighbors=[set() for _ in row];edges=[];vectors=[]
    for i,j in combinations(range(n),2):
        d=gcd(row[i],row[j]);a=row[i]//d;b=row[j]//d
        if a+b in SUMS:
            neighbors[i].add(j);neighbors[j].add(i)
            edges.append((i,j,a,b))
            w=[0]*n;w[i]=b;w[j]=-a;vectors.append(w)
    unseen=set(range(n));components=[]
    while unseen:
        todo=[min(unseen)];seen=set()
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i);todo.extend(neighbors[i]-seen)
        unseen-=seen;components.append(sorted(seen))
    return sorted(components,key=lambda x:(len(x),x)),edges,vectors

def rank(rows):
    basis={}
    for v in rows:
        v=list(v)
        for p in sorted(basis):
            if v[p]:
                a=v[p];b=basis[p][p]
                v=[b*x-a*y for x,y in zip(v,basis[p])]
                d=cgcd(abs(x) for x in v)
                if d:v=[x//d for x in v]
        p=next((i for i,x in enumerate(v) if x),None)
        if p is not None:
            d=cgcd(abs(x) for x in v)
            if v[p]<0:d=-d
            basis[p]=[x//d for x in v]
    return len(basis)

def support(row,pair,k):
    i,j=sorted(pair,key=lambda i:row[i])
    D=gcd(row[i],row[j]);a=row[i]//D;b=row[j]//D
    need(b<=Q,'minimal-coefficient lemma is used only after height gate')
    delta=gcd(D,row[k]);c=D//delta;x=row[k]//delta
    r0=pow(a,-1,b)*x%b;s0=(x-a*r0)//b
    lower=max(ceildiv(-Q-r0,b),ceildiv(s0-Q,a))
    upper=min((Q-r0)//b,(s0+Q)//a)
    yes=c<=Q and lower<=upper
    answer=dict(pair=[i,j],distinguished=k,D=D,a=a,b=b,delta=delta,c=c,x=x,
                lower=lower,upper=upper,crossing=yes)
    if yes:
        r=r0+lower*b;s=s0-lower*a;w=[0]*13
        w[i]=-r;w[j]=-s;w[k]=c;answer['witness']=w
    return answer

def classify(row):
    base=dict(row=list(row))
    if len(row)!=13 or any(type(x) is not int or x<=0 for x in row):
        return dict(base,domain=False,reason='thirteen_positive_integers_required')
    if len(set(row))!=13 or cgcd(row)!=1:
        return dict(base,domain=False,reason='distinct_primitive_row_required')
    if sum(row)>Q*Q:
        return dict(base,domain=False,reason='outside_physical_box')
    comps,edges,vecs=graph(row);c=len(comps)
    base.update(domain=True,components=comps,decoder_rank=13-c)
    if c==1:
        return dict(base,equality=True,reason='connected_rank_twelve',W_rank=12)
    if c>=3:
        return dict(base,equality=False,reason='at_least_three_components',W_rank_interval=[11,12])
    bad=[]
    for part in comps:
        for i,j in combinations(part,2):
            h=max(row[i],row[j])//gcd(row[i],row[j])
            if h>Q:bad.append(dict(pair=[i,j],height=h))
    if bad:
        return dict(base,equality=False,reason='internal_height_failure',W_rank=12,failed_heights=bad)
    I,J=comps;tests=[]
    for first,second in ((I,J),(J,I)):
        for pair in combinations(first,2):
            for k in second:tests.append(support(row,pair,k))
    yes=any(t['crossing'] for t in tests)
    return dict(base,equality=not yes,reason='complete_mixed_supports',W_rank=12 if yes else 11,
                tests=tests,mixed_count=len(tests))

def complement_member(a,b,x):
    if x>Q*(a+b):return False
    M=Q*(a+b)-x
    v0=0 if a==1 else pow(b,-1,a)*M%a
    u0=(M-b*v0)//a
    return max(0,ceildiv(u0-2*Q,b))<=min((2*Q-v0)//a,u0//b)

def validate(row,answer):
    need(answer['domain'],'control lies in actual domain')
    comps,edges,vecs=graph(row)
    need(rank(vecs)==13-len(comps),'literal decoder weighted-edge rank')
    need(all(max(abs(x) for x in v)<=355 for v in vecs),'all decoder rows lie in W')
    for w in vecs:need(sum(x*y for x,y in zip(w,row))==0,'literal edge relation')
    if 'tests' in answer:
        I,J=comps
        required={(tuple(sorted(pair,key=lambda i:row[i])),k)
                  for first,second in ((I,J),(J,I)) for pair in combinations(first,2) for k in second}
        seen=set()
        for test in answer['tests']:
            key=(tuple(test['pair']),test['distinguished'])
            need(key in required and key not in seen,'complete mixed supports with no duplicates')
            seen.add(key)
            alt=test['c']<=Q and complement_member(test['a'],test['b'],test['x'])
            need(alt==test['crossing'],'independent complement-semigroup membership')
            if alt:
                w=test['witness']
                need(sum(x*y for x,y in zip(w,row))==0 and max(map(abs,w))<=Q and sum(x!=0 for x in w)<=3,
                     'returned literal mixed relation')
                need(sum(row[i]*w[i] for i in I)!=0,'mixed relation escapes decoder span')
        need(seen==required,'all mixed supports retained')
        need(len(seen)==11*len(I)*len(J)//2<=231,'exact sharp mixed-triple count')
    return vecs

def summary(name,answer):
    return dict(name=name,component_sizes=list(map(len,answer['components'])),
                decoder_rank=answer['decoder_rank'],reason=answer['reason'],equality=answer['equality'],
                W_rank=answer.get('W_rank'),W_rank_interval=answer.get('W_rank_interval'),
                mixed_count=answer.get('mixed_count'))

def main():
    atlas=[(a,b) for s in SUMS for a in range(1,(s+1)//2) if gcd(a,s)==1 for b in (s-a,)]
    need(len(atlas)==5855 and max(b for a,b in atlas)==355,'actual all-scale atlas')
    records=[]
    row=list(range(1,14));answer=classify(row);vecs=validate(row,answer)
    need(answer['equality'] and answer['W_rank']==12 and len(answer['components'])==1,'connected automatic equality')
    records.append(summary('connected',answer))
    # All six possible unordered two-component size types, with a separate
    # dominance proof of equality rather than only the membership routine.
    total_tests=0
    for a in range(1,7):
        b=13-a;A=[3**j for j in range(a)];K=max(A);g=2*Q*K+1
        row=A+[g*3**j for j in range(b)]
        answer=classify(row);vecs=validate(row,answer)
        need(sorted(map(len,answer['components']))==[a,b],'actual prescribed component sizes')
        need(answer['equality'] and answer['W_rank']==11,'all six types have actual equality controls')
        need(g>2*Q*K,'nonzero large-component partial sum dominates at most two small labels')
        need(answer['mixed_count']==11*a*b//2,'exact count for each size type')
        need(sum(row)<=Q*3**12+3**13<Q*Q,'uniform physical-box bound for all six controls')
        # Repeat after a label rotation/reversal; no component ordering gauge.
        moved=list(reversed(row[3:]+row[:3]));other=classify(moved)
        need(other['domain'] and other['equality'] and other['mixed_count']==answer['mixed_count'],
             'label-order equivariance including singleton orientation')
        total_tests+=answer['mixed_count'];records.append(summary('two_components_'+str(a)+'_'+str(b),answer))
    need(total_tests==1001,'complete six-size-type support count')
    U=[1,4,6,8,10,12,14,15,16,18,22];t=3*Q+1
    row=[t*x for x in U]+[1,3];answer=classify(row);vecs=validate(row,answer)
    positive=[test for test in answer['tests'] if test['crossing']]
    need(not answer['equality'] and len(positive)==1 and rank(vecs+[positive[0]['witness']])==12,
         'mixed-support rejection is a literal rank-twelve obstruction')
    records.append(summary('opposite_orientation_hostile',answer))
    A=[355**j for j in range(6)]+[3,9,27,81,243];K=max(A);g=1000*K+1
    row=A+[g,3*g];answer=classify(row);vecs=validate(row,answer)
    w=[0]*13;w[0]=-1;w[5]=-1000;w[11]=1
    need(answer['reason']=='internal_height_failure' and answer['W_rank']==12,'in-box internal-height branch')
    need(sum(x*y for x,y in zip(w,row))==0 and rank(vecs+[w])==12,'explicit internal-height crossing witness')
    records.append(summary('height_gate_hostile',answer))
    # Three actual decoder components, but exactly rank eleven in W.
    A=[3**j for j in range(10)]+[1000];K=max(A);g=2*Q*K+1
    row=A+[g,3*g];answer=classify(row);vecs=validate(row,answer)
    need(sorted(map(len,answer['components']))==[1,2,10],'actual three-component hostile')
    w=[0]*13;w[0]=-1000;w[10]=1
    need(sum(x*y for x,y in zip(w,row))==0 and max(map(abs,w))<=Q,'new internal bounded row omitted by decoder atlas')
    need(rank(vecs+[w])==11,'lower bound eleven attained by literal vectors')
    need(g>2*Q*K,'all crossings to large pair excluded, giving W rank at most eleven')
    need(not answer['equality'] and answer['W_rank_interval']==[11,12] and 'W_rank' not in answer,
         'three-component branch does not falsely assert terminal rank twelve')
    records.append(dict(summary('three_components_W_rank_eleven',answer),independently_proved_W_rank=11))
    # Three graph components can also have rank twelve in W.
    A=[3**j for j in range(11)];K=max(A);g=Q*K+1;h=2*g+1
    row=A+[g,h];answer=classify(row);vecs=validate(row,answer)
    need(sorted(map(len,answer['components']))==[1,1,11],'second actual three-component control')
    w1=[0]*13;w1[0]=-1;w1[10]=-Q;w1[11]=1
    w2=[0]*13;w2[0]=-1;w2[11]=-2;w2[12]=1
    need(all(sum(x*y for x,y in zip(w,row))==0 and max(map(abs,w))<=Q for w in (w1,w2)),
         'two literal mixed rows are within support and height budget')
    need(rank(vecs+[w1,w2])==12,'upper endpoint twelve is also attained in the multi-component branch')
    records.append(dict(summary('three_components_W_rank_twelve',answer),independently_proved_W_rank=12))
    # Its necessity relies on the physical box, just as in the 11+2 theorem.
    A=[355**j for j in range(11)];g=2*Q*max(A)+1
    out=classify(A+[g,3*g])
    need(not out['domain'] and 'equality' not in out,'outside-box is not an equality verdict')
    print('STATUS: PASS; general thirteen-row decoder equality classification')
    print('TWO-COMPONENT COUNTS:',[11*a*(13-a)//2 for a in range(1,7)],'; maximum 231; six controls total',total_tests)
    for record in records:print(json.dumps(record,sort_keys=True))
    print('SCOPE: W=V_dec only; no safety conclusion; >=3 components permit W ranks 11 and 12')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':main()
