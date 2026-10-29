"""Exact relaxed-threshold volume germs and owner-sensitive continuation."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from functools import reduce
from math import gcd
import json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,m):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(m)
def components(V,deltas):
    result=[(F(0),F(1))]
    for v,d in zip(V,deltas):
        allowed=[((F(k)+d)/v,(F(k+1)-d)/v) for k in range(v)]
        new=[]
        for a,b in result:
            for c,d in allowed:
                lo,hi=max(a,c),min(b,d)
                if lo<=hi:new.append((lo,hi))
        result=[]
        for a,b in sorted(new):
            if result and a<=result[-1][1]:result[-1]=(result[-1][0],max(b,result[-1][1]))
            else:result.append((a,b))
    return result
def measure(pieces):return sum((b-a for a,b in pieces),F(0))
def distance(x):
    r=x-x.numerator//x.denominator
    return min(r,1-r)
def owners(V,t,delta):
    pos=[];neg=[]
    for i,v in enumerate(V):
        r=v*t-(v*t).numerator//(v*t).denominator
        if r==delta:pos.append(i)
        if r==1-delta:neg.append(i)
    return pos,neg
def germ(V,delta,rates=None):
    if rates is None:rates=[F(1)]*len(V)
    pieces=components(V,[delta]*len(V));D=F(0);ends=[]
    for a,b in pieces:
        P,_=owners(V,a,delta);_,N=owners(V,b,delta)
        need(bool(P) and bool(N),'both outward endpoint owner sets present')
        left=min(F(rates[i],V[i]) for i in P)
        right=min(F(rates[i],V[i]) for i in N)
        D+=left+right
        ends.append({'left':str(a),'right':str(b),'positive_owners':[V[i] for i in P],
                     'negative_owners':[V[i] for i in N],
                     'left_rate':str(left),'right_rate':str(right)})
    return measure(pieces),D,pieces,ends
def tight_update(records,w,delta):
    result=[]
    for t,P,N in records:
        d=distance(w*t)
        if d<delta:continue
        a,b=owners((w,),t,delta)
        result.append((t,max(P,w) if a else P,max(N,w) if b else N))
    return result
def tight_records(V,delta):
    pieces=components(V,[delta]*len(V))
    if any(a!=b for a,b in pieces):return None
    return [(a,max(V[i] for i in owners(V,a,delta)[0]),max(V[i] for i in owners(V,a,delta)[1])) for a,b in pieces]
def depth_histogram(V,delta):
    events={F(0),F(1)}
    for v in V:
        for k in range(v):events.update(((F(k)+delta)/v,(F(k+1)-delta)/v))
    events=sorted(events);hist=[F(0)]*(len(V)+1)
    for a,b in zip(events,events[1:]):
        t=(a+b)/2;depth=sum(distance(v*t)<delta for v in V);hist[depth]+=b-a
    need(sum(hist)==1,'independent event depth mass')
    return hist
def gibbs_control(V,delta):
    q=delta.denominator;M=max(V);epsilon=F(1,4*q*M);z=F(1,4*q*M*M)
    hist=depth_histogram(V,delta-epsilon);Z=sum(a*z**j for j,a in enumerate(hist))
    good=bool(components(V,[delta]*len(V)))
    need(hist[0]==measure(components(V,[delta-epsilon]*len(V))),'event depth versus interval measure')
    need((Z>z)==good,'single height-adapted Gibbs certificate for weak witnesses')
    need(Z>=2*z if good else Z<=z,'Gibbs quantitative separation margin')
    return {'epsilon':str(epsilon),'fugacity':str(z),'Z':str(Z),'good':good,
            'depth_masses':[str(v) for v in hist]}
def main():
    cases=0;tight_cases=0;direction_cases=0;gibbs_cases=0;seen={};collision=None
    for n in range(1,5):
        for V in combinations(range(1,11),n):
            if reduce(gcd,V)!=1:continue
            for delta in sorted({F(1,n+1),F(1,3),F(1,4),F(1,2)}):
                L,D,pieces,ends=germ(V,delta);M=max(V);q=delta.denominator
                epsilon=F(1,4*q*M)
                relaxed=measure(components(V,[delta-epsilon]*n))
                need(relaxed==L+epsilon*D,'explicit height-dependent exact first jet')
                need(bool(pieces)==(relaxed>0),'weak witness iff positive relaxed volume')
                need(not pieces or relaxed>=F(1,2*q*M*M),'quantitative positive-volume gap')
                if n<=3:
                    gibbs_control(V,delta);gibbs_cases+=1
                    small=F(1,4*q*M*M)
                    observed=[]
                    for S in range(1<<n):
                        rates=[V[i]*((S>>i)&1) for i in range(n)]
                        _,E,_,_=germ(V,delta,rates)
                        direct=measure(components(V,[delta-small*x for x in rates]))
                        need(direct==L+small*E,'directional owner minimum law')
                        expected=sum(all(S>>i&1 for i in owners(V,a,delta)[0])+
                                     all(S>>i&1 for i in owners(V,b,delta)[1]) for a,b in pieces)
                        need(E==expected,'Boolean directional germ counts contained owner sets')
                        observed.append(int(E));direction_cases+=1
                    # Boolean-lattice inversion reconstructs endpoint-owner multiplicities.
                    inverse=observed[:]
                    for bit in range(n):
                        for S in range(1<<n):
                            if S>>bit&1:inverse[S]-=inverse[S^(1<<bit)]
                    direct=[0]*(1<<n)
                    for a,b in pieces:
                        for O in (owners(V,a,delta)[0],owners(V,b,delta)[1]):
                            direct[sum(1<<i for i in O)]+=1
                    need(inverse==direct,'Mobius inversion recovers side-owner multiset')
                R=tight_records(V,delta)
                if R is not None:
                    for w in range(1,16):
                        updated=tight_update(R,w,delta)
                        literal=tight_records(tuple(V)+(w,),delta)
                        need(updated==literal,'tight-stratum attachment closure')
                        tight_cases+=1
                key=(delta,L,D)
                if collision is None and key in seen:
                    old=seen[key]
                    if components(old,[delta]*len(old))!=pieces:
                        for w in range(1,21):
                            a=components(old+(w,),[delta]*(len(old)+1))
                            b=components(V+(w,),[delta]*(n+1))
                            if bool(a)!=bool(b):
                                collision={'V':old,'W':V,'delta':str(delta),'jet':[str(L),str(D)],
                                           'added_speed':w,'future_nonempty':[bool(a),bool(b)],
                                           'future_components':[[[str(x),str(y)] for x,y in a],[[str(x),str(y)] for x,y in b]]}
                                break
                seen.setdefault(key,V);cases+=1
    examples=[]
    for V,delta in [((1,2),F(1,3)),((1,2),F(1,4)),((1,2),F(2,5)),
                    (tuple(range(1,14)),F(1,14)),((2,6,10),F(1,2)),
                    (tuple(26*j for j in range(1,13))+(339,),F(1,14))]:
        L,D,pieces,ends=germ(V,delta)
        epsilon=F(1,4*delta.denominator*max(V))
        need(measure(components(V,[delta-epsilon]*len(V)))==L+epsilon*D,'large explicit control jet')
        examples.append({'speeds':V,'delta':str(delta),'bulk':str(L),'slope':str(D),
                         'epsilon':str(epsilon),'components':[[str(a),str(b)] for a,b in pieces],
                         'owners':ends,'gibbs':gibbs_control(V,delta)})
        gibbs_cases+=1
    need(collision is not None,'scalar-germ continuation hostile found')
    # Even the complete local directional germ loses the phase of its atoms.
    left,right=(1,3,6),(2,3,6);delta=F(1,3)
    A=germ(left,delta);B=germ(right,delta)
    need(A[:2]==B[:2]==(F(0),F(1)),'tight phase-blind germ coincidence')
    for mask in range(8):
        rates=[F((mask>>i)&1) for i in range(3)]
        need(germ(left,delta,rates)[:2]==germ(right,delta,rates)[:2],'all Boolean directional germs agree')
    for rates in product((-1,0,1),repeat=3):
        h=F(1,10000);thresholds=[delta-h*x for x in rates]
        expected=2*max(F(0),h*F(rates[1],3)+h*F(rates[2],6))
        need(measure(components(left,thresholds))==measure(components(right,thresholds))==expected,
             'full signed local threshold-volume germ formula')
    need(not components(left+(7,),[delta]*4) and bool(components(right+(7,),[delta]*4)),
         'new speed seven separates phase-blind tight owners')
    phase_hostile={'V':left,'W':right,'delta':str(delta),'jet':['0','1'],'added_speed':7,
                   'future_nonempty':[False,True],
                   'components':[[[str(a),str(b)] for a,b in x[2]] for x in (A,B)],
                   'owner_records':[x[3] for x in (A,B)]}
    for V,delta in [((1,2),F(1,3)),((1,3),F(1,3)),((2,3,6),F(1,3))]:
        reference=germ(V,delta)[:2]
        for scale in range(2,6):
            need(germ(tuple(scale*v for v in V),delta)[:2]==reference,'uniform germ dilation invariance')
    out={'status':'FINITE-EXACT','primitive_speed_bound':10,'sizes':[1,2,3,4],
         'base_cases':cases,'direction_cases':direction_cases,'tight_attachment_cases':tight_cases,
         'gibbs_cases':gibbs_cases,
         'scalar_germ_collision':collision,'phase_blind_owner_collision':phase_hostile,
         'examples':examples,'checks':CHECKS}
    (HERE/'lrc_germs.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='examples'},indent=2))
if __name__=='__main__':main()
