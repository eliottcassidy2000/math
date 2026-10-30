"""Independent raw B2 circuit enumeration, actual gate constructions, row tests.
Does not import either producer's implementation. All checks survive python -O.
"""
from fractions import Fraction as Q
from itertools import product
import json

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def gate(a,b,table,width):
    # Deliberately pointwise: no bitwise truth-table synthesis formula.
    return sum(((table >> (((a>>x)&1) + 2*((b>>x)&1)))&1)<<x for x in range(width))

def base_certificate():
    ins=tuple(sum(((x>>i)&1)<<x for x in range(8)) for i in range(3))
    states={ins}; first={a:0 for a in ins}; counts=[]
    for depth in range(1,4):
        following=set()
        for signals in states:
            values={gate(a,b,t,8) for a,b in product(signals, repeat=2) for t in range(16)}
            for value in values:
                first.setdefault(value,depth)
                if depth<3:
                    following.add(signals+(value,))
        counts.append(len(following) if depth<3 else None)
        states=following
    costs=[]
    for r in range(8):
        expected=2 if r.bit_count()<=1 else r.bit_count()+1
        for v in range(8):
            for c in range(2):
                f=sum((c^((r&x).bit_count()%2)^(x==v))<<x for x in range(8))
                check(first.get(f,4)==expected, ('base',r,v,c,first.get(f)))
        if r in (0,1,3,7): costs.append(first.get(sum((((r&x).bit_count()%2)^(x==7))<<x for x in range(8)),4))
    check(len(first)==232,'232 raw reachable functions')
    return {'ordered_states':counts,'reachable_truth_tables':len(first),'canonical_costs':costs}

def constructor(n,k):
    # Gate truth tables follow output bit index a+2b; no free negations.
    gs=[]
    def add(a,b,t):
        gs.append((a,b,t));return n+len(gs)-1
    def chain(wires,t):
        out=wires[0]
        for b in wires[1:]:out=add(out,b,t)
        return out
    if k==0:
        out=chain(list(range(n)),8)
    else:
        rest=chain(list(range(1,n)),8)
        out=add(0,rest,2) # a AND NOT b
        if k>=2:
            parity=chain(list(range(1,k)),6)
            out=add(parity,out,6)
    return gs,out

def construction_checks():
    evaluations=0
    for n in range(3,11):
        for k in range(n+1):
            gs,out=constructor(n,k)
            check(len(gs)==(n-1 if k<=1 else n+k-2),'gate count')
            for x in range(1<<n):
                wires=[(x>>i)&1 for i in range(n)]
                for a,b,t in gs:wires.append((t>>(wires[a]+2*wires[b]))&1)
                check(wires[out]==(((x&((1<<k)-1)).bit_count()%2)^(x==(1<<n)-1)),'constructor')
                evaluations+=1
    # Independent three-response obstruction; additive constants do not matter.
    for s,t in product(range(2),repeat=2):
        responses={(a^b^s,a^b^(a&b)^t) for a,b in product(range(2),repeat=2)}
        check(len(responses)==3,'branching obstruction')
    return evaluations

def potential_checks():
    # All row types, in order (ds,dI,dq,da), including GP k=1..9.
    rows=[(-1,0,0,0),(0,0,0,0),(0,0,-1,0),(-2,1,-1,0),(-1,0,-1,0),(-3,-1,-1,-1),(-3,0,-2,-1)]
    rows += [(-(1+k),-2,1,-1) for k in range(1,10)]
    for j in range(251):
        c=Q(5,2)+Q(3*j,500);i=3*c-7;p=2*c-5;t=4*c-7
        check(i>=0 and p>=0,'nonnegative weights')
        check(t-i==c,'coefficient')
        for ds,di,dq,da in rows:check(ds+i*di+p*dq<=t*da,'row inequality')
        check(2+2*i-p==t and 3+2*p==t,'binding rows')
    # Integer feasible weights test the derived necessary bounds too.
    feasible=0
    for i,p,t in product(range(13),range(13),range(21)):
        if all(ds+i*di+p*dq<=t*da for ds,di,dq,da in rows):
            c=t-i
            check(c<=4 and t>=4*c-7,'necessary bounds')
            feasible+=1
    return {'rational_frontier_points':251,'integer_feasible_weight_triples':feasible,'scope':'archived local row inequalities only'}

if __name__=='__main__':
    print(json.dumps({'base':base_certificate(),'explicit_gate_evaluations':construction_checks(),'potential':potential_checks()},sort_keys=True,indent=2))
