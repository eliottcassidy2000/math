"""Bounded rational lassos, corrected zero-drift scope, and trace/carry loss.

No lasso search is assumed to terminate. A pending frontier is not periodicity.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product

CHECKS=0


def need(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(message)


def integer(n,minimum=None):
    need(type(n) is int and (minimum is None or n>=minimum),'exact integer domain')
    return n


def rational(x):
    need(type(x) in (int,F),'exact rational input')
    x=F(x)
    need(x.denominator%2==1,'odd denominator required')
    return x


def parity(x):
    return rational(x).numerator%2


def step(x):
    x=rational(x)
    return (3*x+1)/2 if x.numerator%2 else x/2


def power3(k):
    integer(k)
    return F(3**k) if k>=0 else F(1,3**(-k))


@dataclass(frozen=True)
class Lasso:
    source: F
    child: F
    debt: int
    prefix: tuple
    cycle: tuple


def audit_lasso(witness):
    need(type(witness) is Lasso,'exact Lasso')
    need(type(witness.source) is F and type(witness.child) is F,'normalized rational roots')
    rational(witness.source);rational(witness.child);integer(witness.debt)
    need(type(witness.prefix) is tuple and type(witness.cycle) is tuple and
         bool(witness.cycle),'finite prefix and nonempty cycle')
    states=witness.prefix+witness.cycle
    for state in states:
        need(type(state) is tuple and len(state)==2,'exact pair state')
        for x in state:
            need(type(x) is F,'normalized exact pair coordinate');rational(x)
    need(states[0]==(witness.source,witness.child),'same supplied limit pair')
    need(len(set(states))==len(states),'canonical first repetition')
    for i,(u,v) in enumerate(states):
        target=states[i+1] if i+1<len(states) else witness.cycle[0]
        need((step(u),step(v))==target,'exact pair edge including cycle closure')
    debt=witness.debt
    source_count=child_count=0
    absorption=None
    for i,(u,v) in enumerate(states):
        if u==v and debt==0 and absorption is None:absorption=i
        if i==len(witness.prefix):entry_debt=debt
        pu,pv=parity(u),parity(v)
        if i>=len(witness.prefix):source_count+=pu;child_count+=pv
        debt+=pu-pv
    delta=source_count-child_count
    if absorption is not None:kind='absorb'
    elif delta:kind='nonzero-drift'
    elif witness.cycle[0][0]==witness.cycle[0][1]:kind='anchored'
    else:kind='zero-drift-pair'
    return dict(kind=kind,entry=len(witness.prefix),period=len(witness.cycle),
                entry_debt=entry_debt,cycle_source_odd=source_count,
                cycle_child_odd=child_count,delta=delta,
                drift=F(delta,len(witness.cycle)),absorption=absorption)


def read_lasso(source,child,debt,budget):
    """At most budget exact transitions, retaining a pending frontier."""
    source=rational(source);child=rational(child)
    integer(debt);integer(budget,0)
    u,v=source,child;k=debt
    seen={};states=[]
    for time in range(budget+1):
        if u==v and k==0:
            return dict(status='absorbed',time=time,states=tuple(states)+((u,v),),debt=k)
        pair=(u,v)
        if pair in seen:
            entry=seen[pair]
            witness=Lasso(source,child,debt,tuple(states[:entry]),tuple(states[entry:]))
            result=audit_lasso(witness)
            return dict(status='lasso',witness=witness,**result)
        if time==budget:
            return dict(status='pending',time=time,states=tuple(states)+(pair,),debt=k)
        seen[pair]=time;states.append(pair)
        k+=parity(u)-parity(v);u,v=step(u),step(v)
    raise ArithmeticError('unreachable bounded-reader branch')


def cycle_values(seed,budget=1000):
    """Bounded single-orbit reader used only by the declared finite census."""
    seed=rational(seed);integer(budget,1)
    states=[];seen={};x=seed
    for _ in range(budget):
        if x in seen:return tuple(states[seen[x]:])
        seen[x]=len(states);states.append(x);x=step(x)
    raise ValueError('single-orbit cycle not certified within budget')


def carrier(word):
    need(type(word) is tuple,'exact valuation tuple')
    p,q,b=1,1,0
    for a in word:
        integer(a,1);p,b,q=3*p,3*b+q,q*(1<<a)
    return p,q,b


def finite_cycles():
    seen=set();cycles=[]
    for length in range(1,5):
        for word in product(range(1,9),repeat=length):
            if sum(word)>8:continue
            if any(length%d==0 and word==word[:d]*(length//d) for d in range(1,length)):continue
            p,q,b=carrier(word);c=F(b,q-p)
            cyc=cycle_values(c)
            key=frozenset(cyc)
            if key not in seen:
                seen.add(key);cycles.append((word,c,cyc))
    return tuple(cycles)+(((0,),F(0),(F(0),)),)


def reject(fn,*args):
    try:fn(*args)
    except (ValueError,TypeError):need(True,'hostile rejected')
    else:need(False,'hostile accepted')


def main():
    # The false finiteness implication is not replaced by assumed termination.
    for n in range(1,64,2):
        need(F(-1)+n+1==n,'arbitrary integer appears as the other limit')
        x=F(n)
        for _ in range(12):
            y=step(x);need(x.denominator%y.denominator==0,'denominators divide the initial denominator')
            x=y
    pending=read_lasso(27,-1,0,10)
    need(pending['status']=='pending' and len(pending['states'])==11,'budget expiry is honestly pending')

    simple=read_lasso(F(19,37),F(23,37),0,10)
    need(simple['status']=='lasso' and simple['entry']==0 and simple['period']==6 and
         simple['delta']==0 and simple['kind']=='zero-drift-pair','different-cycle zero drift')
    need(set(cycle_values(F(19,37))).isdisjoint(cycle_values(F(23,37))),
         'two cycles really are disjoint')
    actual=read_lasso(27*F(7,55)-26,F(7,55),3,100)
    need((actual['entry'],actual['entry_debt'],actual['period'],actual['delta'])==(24,8,12,0),
         'incoming universal-state table row')
    need(actual['witness'].cycle[0]==(F(2,55),F(7,55)),'exact lasso entry coordinates')
    need(len(cycle_values(F(2,55)))==12 and len(cycle_values(F(7,55)))==6,
         'same drift does not identify a cycle')

    # Every finite census outcome is established by its own bounded exact witness.
    cycles=finite_cycles();need(len(cycles)==55,'declared primitive-word cycle universe')
    states=((3,F(-26)),(4,F(10)),(3,F(13)),(-5,F(1,243)-1),(1,F(1)))
    counts={};distinct_zero=0;largest=0
    for k,e in states:
        for _,c,_ in cycles:
            for child_runs in (True,False):
                u,v=(power3(k)*c+e,c) if child_runs else (c,(c-e)/power3(k))
                result=read_lasso(u,v,k,2000)
                need(result['status']!='pending','finite census has an authenticated outcome')
                name='absorbed' if result['status']=='absorbed' else result['kind']
                counts[name]=counts.get(name,0)+1
                used=result['time'] if name=='absorbed' else result['entry']+result['period']
                largest=max(largest,used)
                if name=='zero-drift-pair':
                    a,b=result['witness'].cycle[0]
                    distinct_zero+=set(cycle_values(a))!=set(cycle_values(b))
    need(sum(counts.values())==550,'five states times55cycles times two run sides')

    p,q,b=carrier((1,2));pp,qq,bb=carrier((2,1))
    need((p,q,b)==(9,8,5) and (pp,qq,bb)==(9,8,7),'trace-equal carry-distinct words')
    need(p+q==pp+qq and F(p,q)==F(pp,qq),'trace and multiplier forget order')
    need(((q-b)*pow(p,-1,2*q)%(2*q),(qq-bb)*pow(pp,-1,2*qq)%(2*qq))==(11,9),
         'different forward native cylinders')
    need(b*pow(q,-1,p)%p==4 and bb*pow(qq,-1,pp)%pp==2,'different inverse endpoint guards')
    for E in range(1,61):
        n=2**E-1
        need(((q*n-b)%p==0)==(E%6==5),'only word12 Mersenne inverse phase')
        need((qq*n-bb)%pp!=0,'word21 has no Mersenne inverse phase')
    for bad in (True,1.0,F(1,2)):reject(read_lasso,bad,-1,0,10)
    reject(read_lasso,3,-1,True,10);reject(read_lasso,3,-1,0,-1)
    w=simple['witness']
    reject(audit_lasso,replace(w,debt=True))
    reject(audit_lasso,replace(w,cycle=w.cycle[:-1]))
    reject(audit_lasso,replace(w,source=F(21,37)))
    print('Universal periodicity is NOT deduced from bounded rational denominators.')
    print('Bounded n27/child-1 reader: PENDING after10 transitions; no claim about eventual fate.')
    print('Simple disjoint zero-drift cycles:19/37 versus23/37; period6, equal odd count3.')
    print('Universal childrun(2,4): entry24,debt8,pair(2/55,7/55),jointperiod12,drift0; distinct cycles12/6.')
    print('FINITE-EXACT lasso census:55cycles,5states,2sides; outcomes',sorted(counts.items()))
    print('Zero-drift outcomes on distinct cycles:',distinct_zero,'largest verified transition count:',largest)
    print('Trace12/21:P9,Q8,B5/7; forward native11/9mod16; only12 inverse Mersenne phase E5mod6.')
    print('No universal lasso termination, no broad ladder-completeness audit, no trace-only guard inference.')
    print('Exact checks:',CHECKS)


if __name__=='__main__':main()
