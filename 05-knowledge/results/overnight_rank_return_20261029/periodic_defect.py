"""A finite bulk transfer and exact recurrence for actual (2,2,...,2) targets."""
from collections import defaultdict,deque
import json
from sympy import Matrix,Poly,symbols,rem,gcd,invert,factorint
from collections import Counter
from math import prod
from complement_kernel import HERE,need
from defect_transfer import query,interval_target


def step(states,v,back,last,start=None,end=None):
    out=defaultdict(int)
    for (done,components),weight in states.items():
        lefts=[None] if start is not None else [None]+[i for i,(h,t) in enumerate(components) if t==-1 or t>=0 and t not in back]
        rights=[None] if end is not None else [None]+[i for i,(h,t) in enumerate(components) if h>=0 and h in back]
        for a in lefts:
            for b in rights:
                if a is not None and a==b:continue
                h=components[a][0] if a is not None else -start-2 if start is not None else v
                t=components[b][1] if b is not None else -end-2 if end is not None else v
                pairs=[pair for i,pair in enumerate(components) if i not in (a,b)]
                finished=done
                if h<0 and t<-1:
                    if h!=t or done>>(-h-2)&1:continue
                    finished|=1<<(-h-2)
                else:pairs.append((h,t))
                if any(h>=0 and last[h]<=v for h,t in pairs):continue
                pairs=[(h,-1 if t>=0 and last[t]<=v else t) for h,t in pairs]
                out[finished,tuple(sorted(pairs))]+=weight
    return dict(out)


def shift(state,offset):
    done,components=state
    return done,tuple(sorted((h-offset if h>=0 else h,t-offset if t>=0 else t) for h,t in components))


def initial():
    states={(0,()):1};last={0:3,1:2,2:5,3:4}
    for v in range(4):
        back=set() if v<2 else {1} if v==2 else {0}
        states=step(states,v,back,last,start=v if v<3 else None)
    return {shift(state,2):value for state,value in states.items()}


def bulk(state):
    states={state:1};last={0:3,1:2,2:5,3:4}
    states=step(states,2,{1},last)
    states=step(states,3,{0},last)
    return {shift(s,2):v for s,v in states.items()}


def terminal(state):
    states={state:1};last={0:3,1:2,2:5,3:4,4:-1,5:-1}
    for v,back in ((2,{1}),(3,{0}),(4,{3}),(5,{2})):
        states=step(states,v,back,last,end=v-3 if v>=3 else None)
    return states.get((7,()),0)


def run():
    init=initial();states=list(init);known=set(states);pending=deque(states);edges={}
    while pending:
        state=pending.popleft();edges[state]=bulk(state)
        for target in edges[state]:
            if target not in known:known.add(target);states.append(target);pending.append(target)
    d=len(states);index={s:i for i,s in enumerate(states)}
    B=[[0]*d for _ in states]
    for s,targets in edges.items():
        for t,value in targets.items():B[index[s]][index[t]]=value
    alpha=[init.get(s,0) for s in states];beta=[terminal(s) for s in states]
    vector=alpha[:];values=[]
    for r in range(2*d+2):
        values.append(sum(x*y for x,y in zip(vector,beta)))
        vector=[sum(vector[i]*B[i][j] for i in range(d)) for j in range(d)]
    coefficients=None
    for order in range(1,d+1):
        H=Matrix([[values[i+j] for j in range(order+1)] for i in range(d)])
        kernel=H.nullspace()
        if kernel:
            v=next((v for v in kernel if v[-1]),None)
            if v is not None:
                coefficients=[x/v[-1] for x in v];break
    need(coefficients is not None and all(x.q==1 for x in coefficients),'monic integral scalar recurrence')
    coefficients=list(map(int,coefficients));order=len(coefficients)-1
    # Independent vector proof: beta Krylov polynomial is annihilated by
    # every initial row alpha*B^j for j=0,...,d-1; Cayley-Hamilton closes all j.
    power=beta[:];residual=[0]*d
    for coefficient in coefficients:
        residual=[x+coefficient*y for x,y in zip(residual,power)]
        power=[sum(B[i][j]*power[j] for j in range(d)) for i in range(d)]
    vector=alpha[:]
    for j in range(d):
        need(sum(x*y for x,y in zip(vector,residual))==0,'finite Krylov annihilation certificate')
        vector=[sum(vector[i]*B[i][j] for i in range(d)) for j in range(d)]
    minor=int(Matrix([[values[i+j] for j in range(order)] for i in range(order)]).det())
    need(minor!=0,'nonzero Hankel minor proves minimal order from r=0')
    for r in range(20):
        n,P,D=interval_target([2]*(r+4))
        actual,_=query(n,D,(0,1,2),(n-3,n-2,n-1))
        need(actual==values[r] if r<len(values) else False,'independent boundary query')
    generated=values[:order]
    for r in range(order,497):
        generated.append(-sum(coefficients[j]*generated[r-order+j] for j in range(order)))
    long=json.loads((HERE/'defect_transfer.json').read_text())['long_targets']
    for record in long:need(generated[(record['n']-8)//2]==record['value'],'recurrence agrees at 20,100,1000')
    need(not any(residual),'stronger exact identity p(B)*beta=0')
    x=symbols('x');p=sum(c*x**i for i,c in enumerate(coefficients))
    f6=x**6-x**5-6*x**4+4*x**3-7*x**2+5*x-4
    f11=x**11-13*x**10-16*x**9-268*x**8-382*x**7-1170*x**6-1048*x**5+2072*x**4+3197*x**3+459*x**2-1232*x-576
    need(Poly(p,x)==Poly((x+1)*f6*f11,x),'exact characteristic factorization')
    g=Poly(x**5+x**4-x**3-x**2-1,x,modulus=3)
    need(Poly(p,x,modulus=2)==Poly(x**3*(x+1)**15,x,modulus=2),'binary characteristic layers')
    need(Poly(p,x,modulus=3)==Poly(x*(x+1)**2*g.as_expr()**3,x,modulus=3),'ternary characteristic layers')
    need(rem(Poly(x**243-x,x,modulus=3),g).is_zero and gcd(g,Poly(x**3-x,x,modulus=3)).degree()==0,
         'prime-degree irreducibility criterion for the actual F243 chart')
    need(rem(Poly(x**121,x,modulus=3),g)==Poly(1,x,modulus=3)
         and rem(Poly(x**11,x,modulus=3),g)!=Poly(1,x,modulus=3),'multiplicative order exactly121')
    modulus=g*g
    eps=rem(g*invert(g.diff(),modulus),modulus)
    field_alpha=rem(Poly(x,x,modulus=3)-eps,modulus)
    def ring_power(value,exponent):
        answer=Poly(1,x,modulus=3)
        while exponent:
            if exponent&1:answer=rem(answer*value,modulus)
            value=rem(value*value,modulus);exponent//=2
        return answer
    need(not eps.is_zero and rem(eps*eps,modulus).is_zero,'visible nonzero nilpotent direction')
    need(ring_power(field_alpha,121)==Poly(1,x,modulus=3),'Hensel field lift has order121')
    need(ring_power(Poly(x,x,modulus=3),363)==Poly(1,x,modulus=3)
         and ring_power(Poly(x,x,modulus=3),121)!=Poly(1,x,modulus=3),'clock retains the nilpotent factor3')
    need(ring_power(eps,3).is_zero,'cubing is not an automorphism of the actual local algebra')
    q=Poly((x+1)**2*g.as_expr()**2,x,modulus=3)
    qcoeff=[int(q.nth(j))%3 for j in range(13)]
    need(all(sum(c*values[i+j] for j,c in enumerate(qcoeff))%3==0 for i in range(d)),
         'projected ternary recurrence; Cayley-Hamilton closes all later coefficients')
    det12=int(Matrix([[values[i+j]%3 for j in range(12)] for i in range(12)]).det())%3
    need(det12==1,'ternary scalar recurrence has minimal order12')
    orbits=[]
    for modulus in (2,3,6):
        state=tuple(v%modulus for v in values[:order]);seen={};r=0
        while state not in seen:
            seen[state]=r
            state=state[1:]+((-sum(coefficients[j]*state[j] for j in range(order))%modulus),)
            r+=1
        need(seen[state]==0 and r==(1 if modulus==2 else 726),'exact actual initial-state period')
        orbits.append(dict(modulus=modulus,preperiod=0,period=r))
    from paley_pencil import diagonal,local_layers,valuation
    H=[[values[i+j] for j in range(order)] for i in range(order)]
    ds=diagonal(H)
    need(prod(ds)==abs(minor),'Hankel lattice index independently matches determinant')
    observable=[];power=beta[:]
    for j in range(order):
        observable.append(power)
        power=[sum(B[i][j]*power[j] for j in range(d)) for i in range(d)]
    obs_smith=diagonal(observable)
    obs_layers={}
    for prime in sorted({int(p) for value in obs_smith for p in factorint(value)}):
        exponents=[valuation(value,prime) for value in obs_smith]
        expected=Counter(exponents);precision=max(exponents)+1
        actual=local_layers(observable,prime,precision)
        need(actual==[expected[i] for i in range(precision)],'independent observer precision layers')
        obs_layers[str(prime)]=actual
    print('periodic dimension',d,'minimal order',order,'polynomial',coefficients,'Hankel',minor,flush=True)
    return dict(states=states,transfer=B,initial=alpha,terminal=beta,dimension=d,
                recurrence_low_first=coefficients,values=values,hankel_minor=minor,
                residual=residual,long_indices=[6,46,496],hankel_smith=ds,
                observable=observable,observable_smith=obs_smith,observable_layers=obs_layers,
                ternary_recurrence_low_first=qcoeff,ternary_hankel_minor=det12,orbits=orbits,
                field_polynomial_low_first=[int(g.nth(j))%3 for j in range(6)],
                field_degree=5,multiplicative_order=121,frobenius_order=5,
                hensel_field_lift=[int(field_alpha.nth(j))%3 for j in range(10)],
                nilpotent=[int(eps.nth(j))%3 for j in range(10)])


if __name__=='__main__':
    import sys
    receipt=json.loads(json.dumps(run()))
    path=HERE/'periodic_defect.json'
    if '--write' in sys.argv:path.write_text(json.dumps(receipt,separators=(',',':'))+'\n')
    else:need(receipt==json.loads(path.read_text()),'unchanged finite transfer certificate')
    print('PASS actual-complement periodic recurrence',flush=True)
