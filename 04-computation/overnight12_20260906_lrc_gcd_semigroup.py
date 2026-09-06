"""Exact signed-box and actual-entry LRC gcd addendum controls.

Standard library only, always-active gates, no producer imports. Compare
exact interval membership with complete coefficient boxes, and the minimal
external-coefficient test with all bounded three-coefficient relations.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, isqrt, prod
import sys

sys.stdout.reconfigure(newline="\n")
Q = 91**6
CUTOFF = Q//(42*177)
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(label)


def ceildiv(n,d):
    return -((-n)//d)


def radius(a,b,B):
    return B*(a+b)-(a-1)*(b-1)


def membership(x,a,b,B):
    r0 = (pow(a,-1,b)*x)%b
    s0 = (x-a*r0)//b
    lo = max(ceildiv(-B-r0,b),ceildiv(s0-B,a))
    hi = min((B-r0)//b,(s0+B)//a)
    return lo<=hi


def witness(x,a,b,B):
    r0 = (pow(a,-1,b)*x)%b
    s0 = (x-a*r0)//b
    lo = max(ceildiv(-B-r0,b),ceildiv(s0-B,a))
    hi = min((B-r0)//b,(s0+B)//a)
    if lo>hi:
        return None
    return r0+lo*b,s0-lo*a


def native(u,v,t,g,s,B):
    D = gcd(u,v)
    a,b = u//D,v//D
    delta = gcd(t*D,g*s)
    c,x = t*D//delta,g*s//delta
    R = radius(a,b,B)
    z = gcd(delta,s)
    d,e = delta//z,s//z
    lower = d*(R//e+1)
    crosses = c<=B and membership(x,a,b,B)
    return dict(D=D,a=a,b=b,delta=delta,c=c,x=x,R=R,lower=lower,crosses=crosses)


def clearance(v,x):
    y=v*x
    r=y-y.numerator//y.denominator
    return min(r,1-r)


def pair_phase(p,q):
    s=p+q
    return F(pow(p,-1,s)*(s//2)%s,s)


def safe_lift(V,t,g,p,q,y):
    K=max(V)
    need(min(clearance(v,y) for v in V)>=F(1,12),"literal eleven-core supplier phase")
    need(gcd(t,g)==1 and g>=42*K,"complete grid meets global core arc")
    z=pair_phase(p,q)
    need(min(clearance(p,z),clearance(q,z))>=F(1,3),"literal pair supplier phase")
    left=y-F(1,84*K)
    a=g*left-t*z
    k=ceildiv(a.numerator,a.denominator)
    j=(pow(t,-1,g)*k)%g
    x=(z+j)/g
    target=(t*z+k)/g
    need(left<=target<=y+F(1,84*K),"exact global core-clock arc hit")
    need((t*x-target).denominator==1,"clock residue transported exactly")
    physical=tuple(t*v for v in V)+(g*p,g*q)
    need(min(clearance(w,x) for w in physical)>=F(1,14),"literal full physical clearance")
    return x


def decoder_components(speeds):
    primes=[p for p in range(2,357) if p%3==2 and all(p%d for d in range(2,isqrt(p)+1))]
    allowed={1}
    for p in primes:
        allowed|={s*p**e for s in tuple(allowed) for e in (1,2) if s*p**e<=356}
    adj=[set() for _ in speeds]
    for i,j in combinations(range(len(speeds)),2):
        if (speeds[i]+speeds[j])//gcd(speeds[i],speeds[j]) in allowed:
            adj[i].add(j)
            adj[j].add(i)
    unseen=set(range(len(speeds)))
    out=[]
    while unseen:
        todo=[min(unseen)]
        reached=set()
        while todo:
            i=todo.pop()
            if i in reached:
                continue
            reached.add(i)
            todo.extend(adj[i]-reached)
        unseen-=reached
        out.append(reached)
    return out


def main():
    box_cases=point_cases=0
    for B in range(2,18):
        for b in range(2,B+1):
            for a in range(1,b):
                if gcd(a,b)!=1:
                    continue
                values={r*a+s*b for r in range(-B,B+1) for s in range(-B,B+1)}
                R=radius(a,b,B)
                L=B*(a+b)
                need(2*R>L,"central interval exceeds half of box support")
                need(all(x in values for x in range(-R,R+1)),"complete central interval by literal coefficient set")
                need(R+1 not in values and -(R+1) not in values,"sharp first holes in both directions")
                need(R>=B*b and R<=L,"radius dominates old residue strip and respects support")
                for x in range(-L-2,L+3):
                    actual=x in values
                    need(membership(x,a,b,B)==actual,"integer interval membership versus literal box")
                    w=witness(x,a,b,B)
                    need((w is not None)==actual,"constructive exact witness iff membership")
                    if w is not None:
                        need(w[0]*a+w[1]*b==x and max(map(abs,w))<=B,"literal witness identity and height")
                    point_cases+=1
                for x in range(1,L+2):
                    if x in values:
                        continue
                    need(all(j*x not in values for j in range(1,B+1)),"a missing positive ray point has no bounded positive-multiple rescue")
                box_cases+=1

    relation_cases=0
    for B in range(2,7):
        for b in range(2,B+1):
            for a in range(1,b):
                if gcd(a,b)!=1:
                    continue
                values={r*a+s*b for r in range(-B,B+1) for s in range(-B,B+1)}
                for D in (1,2,5):
                    for t in (1,2,3):
                        A,Bphys=t*D*a,t*D*b
                        rhs={t*D*x for x in values}
                        for outside in range(1,31):
                            literal=any(c*outside in rhs for c in range(1,B+1))
                            delta=gcd(t*D,outside)
                            c=t*D//delta
                            x=outside//delta
                            predicted=c<=B and membership(x,a,b,B)
                            need(predicted==literal,"minimal external coefficient iff full three-label bounded crossing")
                            relation_cases+=1
    need(not any(c in {3*r+6*s for r in range(-2,3) for s in range(-2,3)} for c in (1,2)),
         "external coefficient gate cannot be omitted")
    need(membership(5,2,3,2) and 2*15==12+18,"nontrivial delta-clearing positive control")
    need(radius(2,3,3)==13 and not membership(14,2,3,3) and membership(15,2,3,3),
         "central first hole does not exclude all larger points")
    need(not membership(3,1,6,2) and membership(6,1,6,2),
         "omitting b<=Q invalidates the minimal external-coefficient theorem")

    lattice_cases=0
    for R in range(1,31):
        for delta in range(1,13):
            for s in range(1,13):
                z=gcd(delta,s)
                d,e=delta//z,s//z
                lower=d*(R//e+1)
                literal=next(g for g in range(1,lower+1) if (g*s)%delta==0 and g*s//delta>R)
                need(lower==literal,"strongest integer scale lower bound from central radius and divisibility")
                lattice_cases+=1

    need(Q==567869252041 and CUTOFF==76388115,"exact inherited height and uniform gcd cutoff")
    need(677*CUTOFF==51714753855<Q,"small-t cleared coefficient budget")
    need(Q>=42*177*CUTOFF and Q<42*177*(CUTOFF+1),"uniform cutoff division is exact")
    for p in range(1,178):
        need(F(Q,CUTOFF*p)>=42,"all-parameter maximal-p gluing comparison")
    for q in range(2,356):
        need((21*q-1)//11<=677,"large-t complement is bounded")

    examples=[]
    for V,g,y in [((1,4,6,8,10,12,14,15,16,18,22),2**45,F(9,19)),
                  ((2,3,4,5,6,10,12,14,15,20,30),60*Q+1,F(1,17))]:
        t,p,q=1,1,3
        K=max(V)
        speeds=V+(g,3*g)
        need(len(V)==11 and len(set(speeds))==13 and reduce(gcd,speeds)==1,"distinct primitive physical positive control")
        need(sum(speeds)<=Q**2 and g>2*Q*K,"actual box and complete crossing exclusion by dominance")
        need(decoder_components(speeds)==[set(range(11)),{11,12}],"literal connected 11+2 decoder graph")
        scores=[]
        for u,v in combinations(V,2):
            for s in (p,q):
                data=native(u,v,t,g,s,Q)
                need(data['b']<=Q and not data['crosses'],"all 110 exact selected crossing tests in actual-entry control")
            data=native(u,v,t,g,p,Q)
            scores.append((data['R'],u,v,data['lower']))
        phase=safe_lift(V,t,g,p,q,y)
        examples.append((V,g,str(phase),max(scores)))
        if min(V)>1:
            maximum_pair=max(x for x in scores if x[2]==K)
            interior=max(x for x in scores if x[2]<K)
            need(maximum_pair[:3]==(22*Q-84,14,30),"best maximum-coordinate pair score")
            need(interior[:3]==(29*Q-182,14,15),"strictly stronger arbitrary-pair score")
            need(interior[0]>maximum_pair[0],"arbitrary core pairs retain quantitatively useful information")
            need(all(gcd(g,v)==1 for v in V),"score control has delta one at outside p=1")

    # A pair can improve a relation score, but its own maximum is not the
    # Lipschitz constant for the complete eleven-speed safe arc.
    V=tuple(range(1,11))+(85,)
    y=F(1,12)
    lost=F(7,85)
    need(min(clearance(v,y) for v in V)>=F(1,12),"full core is safe at the proposed center")
    need(abs(lost-y)<F(1,84*10) and clearance(85,lost)==0,"substituting selected-pair maximum loses the full-core predicate")

    primes=(37,43,61,67,73,79,97,103,127)
    P=15*prod(primes)
    primorial=tuple(2*P//r for r in primes)+(P//3,P//5)
    need(reduce(gcd,primorial)==1 and min(primorial)>1,"primitive nonunit normalization hostile")
    need(max(primorial)==237907127334685115>Q,"primitive maximum can exceed Q")
    need(max(max(a,b)//gcd(a,b) for a,b in combinations(primorial,2))==127,"all primitive internal pair heights still small")
    hierarchy=[]
    for size,subset,bound,wanted in [(7,primorial[:7],90,392430),(8,primorial[:8],30,3810),
                                   (9,primorial[:9],9,30),(10,primorial[:9]+(P//3,),4,5)]:
        value=reduce(gcd,subset)
        need(value==wanted>bound,"primorial hostile is already excluded as a strict failure by inherited hierarchy")
        hierarchy.append((size,value,bound))

    print("STATUS: PASS; exact signed coefficient boxes and native actual-entry gcd criteria")
    print("BOX UNIVERSE:",box_cases,"primitive pairs Qtoy2..17;",point_cases,"literal membership points")
    print("FULL CROSSING UNIVERSE:",relation_cases,"three-label coefficient boxes Qtoy2..6")
    print("DISCRETE SCALE UNIVERSE:",lattice_cases,"exact divisibility lower-bound controls")
    print("UNIFORM GCD CUTOFF:",CUTOFF,"at Q=",Q,"with residual coefficient cap",677*CUTOFF)
    for V,g,phase,score in examples:
        print("ACTUAL ENTRY CONTROL:",V,"g=",g,"phase=",phase,"best radius pair=",score)
    print("ARBITRARY PAIRS: interior(14,15) radius29Q-182 > maximum-pair(14,30) radius22Q-84")
    print("GLOBAL ARC HOSTILE: V=(1,...,10,85), center1/12, selected max10 loses safety at7/85")
    print("PRIMORIAL HIERARCHY EXCLUSIONS (subset size,gcd,inherited cap):",hierarchy)
    print("SCOPE: actual W=Vdec and box retained; hierarchy constraints apply only under strict failure")
    print("ACTIVE GATES:",GATES)


if __name__=="__main__":
    main()
