"""Exact, dependency-free controls. Run from this directory with python [-O]."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod
import json


def need(ok, detail):
    if not ok:
        raise ValueError(detail)


def bounded(c):
    if c < -4:
        return []
    M = 2 + isqrt(c + 4)
    M -= (M + 1) % 2
    S = set(range(-M, M + 1, 2))
    while True:
        R = {a for a in S if (a*a-c)//4 in S}
        if R == S:
            return sorted(S)
        S = R


SMALL = {-3: [-3,-1,1,3], 5: [-5,-3,-1,1,3,5],
         13: [-5,-3,-1,1,3,5], 21: list(range(-7,8,2)),
         29: list(range(-7,8,2)), 37: [-7,-5,-3,3,5,7],
         45: [-9,-5,-3,3,5,9]}


def classified(c):
    if c < -3:
        return []
    if c in SMALL:
        return SMALL[c]
    for d, shift in [(c+4, 2), (c-12, -2)]:
        s = isqrt(d)
        if s*s == d:
            return sorted({v for a in [shift-s, shift+s] for v in [a,-a]})
    return []


def cycles(D):
    out = set()
    for x in D:
        path = []
        while x not in path:
            path.append(x)
            x = D[x]
        cyc = path[path.index(x):]
        out.add(min(tuple(cyc[i:]+cyc[:i]) for i in range(len(cyc))))
    return sorted(out)


def layer_data(m):
    V = [(i,j) for i in range(1,m+1) for j in range(i)]
    E = [((i,u),(i+1,v)) for i in range(1,m)
         for u in range(i) for v in range(i+1)]
    ix = {e:k for k,e in enumerate(E)}
    def vec(terms):
        d = [0]*len(E)
        for e,s in terms:
            d[ix[e]] += s
        return d
    R=[]; H=[]
    for i in range(1,m):
        for u in range(1,i):
            for v in range(1,i+1):
                R.append(vec([(((i,u),(i+1,v)),1),
                    (((i,u),(i+1,0)),-1),(((i,0),(i+1,v)),-1),
                    (((i,0),(i+1,0)),1)]))
    for i in range(2,m):
        for u in range(1,i):
            H.append(vec([(((i-1,0),(i,u)),1),
                (((i,u),(i+1,0)),1),(((i-1,0),(i,0)),-1),
                (((i,0),(i+1,0)),-1)]))
    return V,E,R,H


def rank2(rows):
    piv={}
    for row in rows:
        v=sum((x%2)<<i for i,x in enumerate(row))
        while v:
            p=v.bit_length()-1
            if p in piv: v ^= piv[p]
            else:
                piv[p]=v
                break
    return len(piv)


def tournament(j):
    free = {(0,3):0,(0,2):1,(1,3):2}
    return {(a,b): (1 if b==a+1 else (j>>free[a,b])&1)
            for a,b in combinations(range(4),2)}


def arc(T,a,b):
    return T[a,b] if a<b else 1-T[b,a]


def canonical(T):
    return min(tuple(arc(T,p[a],p[b]) for a,b in combinations(range(4),2))
               for p in permutations(range(4)))


def primes_to(n):
    return [p for p in range(2,n+1) if all(p%d for d in range(2,isqrt(p)+1))]


def prime_factors(m):
    out=[];p=2
    while p*p<=m:
        if m%p==0:
            out.append(p)
            while m%p==0:m//=p
        p+=1
    if m>1:out.append(m)
    return out


def wheel(m):
    U={a for a in range(1,m,2) if gcd(a,m)==1}
    P=[p for p in prime_factors(m) if p!=2]
    if not P:
        need(len(U)==m//2,('power two',m))
        return {'cycle':m//2}
    L=m//(2*prod(P));small=min(P)
    def S(r):return L*prod(p-r for p in P) if r<small else 0
    predicted={r:S(r)-2*S(r+1)+S(r+2) for r in range(1,small)}
    predicted={r:v for r,v in predicted.items() if v}
    observed=Counter()
    for a in U:
        if (a-2)%m in U:continue
        r=0;b=a
        while b in U:
            r+=1;b=(b+2)%m
        observed[r]+=1
    need(dict(observed)==predicted,('wheel',m,observed,predicted))
    need(sum(r*v for r,v in predicted.items())==len(U),('wheel coverage',m))
    return predicted


def main():
    supported=0
    for k in range(12502):
        tri=lambda v:v>=0 and isqrt(8*v+1)**2==8*v+1
        total=4*(tri(k)+tri(k-2))+2*(k in [1,2,5,6])+8*(k==4)
        need(total==len(classified(8*k-3)),('generating function',k))
        supported+=total>0
        N=lambda v:0 if v<0 else (isqrt(8*v+1)-1)//2+1
        need(supported==N(k)+N(k-2)-(k>=3)+(k>=4),('support count',k))
    for c in range(-99,100006,8):
        need(bounded(c)==classified(c),('classification',c))
    small={c:{'bounded':bounded(c),
        'cycles':cycles({a:(a*a-c)//4 for a in bounded(c)})} for c in SMALL}
    layer=[]
    for m in range(2,13):
        V,E,R,H=layer_data(m)
        need(len(R)+len(H)==len(E)-len(V)+1,('dimension',m))
        need(rank2(R+H)==len(R)+len(H),('basis rank',m))
        # A basis of non-tree columns has triangular +/-1 matrix.
        non_tree=[i for i,(u,v) in enumerate(E) if u[1]!=0]
        pivots=[]
        for row in H+R:
            support=[i for i in non_tree if row[i]]
            # Eliminate each previously certified unit pivot, over Z.
            row=row[:]
            for pivot,prior in pivots:
                coefficient=row[pivot]//prior[pivot]
                if coefficient:
                    row=[x-coefficient*y for x,y in zip(row,prior)]
            support=[i for i in non_tree if row[i]]
            need(len(support)==1 and abs(row[support[0]])==1,('integral pivot',m))
            pivots.append((support[0],row))
        layer.append([m,len(V),len(E),len(R),len(H)])
    V,E,R,H=layer_data(3)
    gradients={tuple((x[V.index(v)]-x[V.index(u)])%2 for u,v in E)
               for x in product(range(2),repeat=len(V))}
    checked=0
    for L in product(range(2),repeat=len(E)):
        zero=all(sum(x*y for x,y in zip(row,L))%2==0 for row in R+H)
        need(zero==(L in gradients),('syndrome',L));checked+=1
    hostile=[0]*len(E);hostile[E.index(((1,0),(2,1)))]=1
    need(all(sum(a*b for a,b in zip(r,hostile))==0 for r in R), 'local control')
    need(any(sum(a*b for a,b in zip(h,hostile))!=0 for h in H),'hourglass hostile')
    F=[(j-3)*(j-4)//2 for j in range(8)]
    q=lambda j: ((j&1)^((j>>1)&1),(j&1)^((j>>2)&1))
    Q={};trows=[]
    for j in range(8):
        u,v=q(j)
        output=(u,1^u^v,(1^u)&(1^v))
        need(F[j]==sum(b<<i for i,b in enumerate(output)),('Boolean rule',j))
        need(F[7-j]==F[j],('fold',j))
        Q[q(j)]=q(F[j])
        T=tournament(j);C=canonical(T)
        conv=canonical({k:1-v for k,v in T.items()})
        hp=sum(all(arc(T,p[i],p[i+1]) for i in range(3)) for p in permutations(range(4)))
        aut=sum(all(arc(T,p[a],p[b])==T[a,b] for a,b in T) for p in permutations(range(4)))
        trows.append({'j':j,'a':2*j-7,'next':F[j],'lines':q(j),'class':C,
                      'merged':min(C,conv),'H':hp,'Aut':aut})
    fibres=Counter(r['class'] for r in trows);merged=Counter(r['merged'] for r in trows)
    for r in trows:need(fibres[r['class']]*r['Aut']==r['H'],'H/Aut')
    need(sorted(fibres.values())==[1,1,1,5],'ordinary fibres')
    need(sorted(merged.values())==[1,2,5],'merged fibres')
    pair_failures=[]
    for i,j in combinations(range(8),2):
        if trows[i]['class']==trows[j]['class'] and trows[F[i]]['class']!=trows[F[j]]['class']:
            pair_failures.append([i,j,F[i],F[j]])
    need(bool(pair_failures),'class descent hostile')
    merged_failures=[]
    for i,j in combinations(range(8),2):
        if trows[i]['merged']==trows[j]['merged'] and trows[F[i]]['merged']!=trows[F[j]]['merged']:
            merged_failures.append([i,j,F[i],F[j]])
    need(bool(merged_failures),'merged descent hostile')
    state_primes=[2,3,5,7,11,13,17,19];clock=23;tag=29
    program=[Fraction(clock*state_primes[F[j]],p) for j,p in enumerate(state_primes)]
    traces=[]
    for j in range(8):
        n=state_primes[j]*tag**j;state=j;tr=[n]
        for t in range(1,10):
            valid=[f for f in program if (n*f.numerator)%f.denominator==0]
            need(len(valid)==1,('one hot guard',j,t))
            n=n*valid[0].numerator//valid[0].denominator
            state=F[state]
            need(n==state_primes[state]*tag**j*clock**t,('clock source',j,t))
            tr.append(n)
        traces.append(tr)
    # Off-promise multi-control inputs are explicitly not semantic states.
    need(sum((2*3*f.numerator)%f.denominator==0 for f in program)==2,'multi-control hostile')
    for m in range(4,2002,2):wheel(m)
    joint_tests=0
    for m in [8,10,12,18,30,290]:
        P=[p for p in prime_factors(m) if p!=2]
        L=m//(2*prod(P))
        for mask in range(1,128):
            D=[i for i in range(7) if mask>>i&1]
            actual=sum(all(gcd(a+2*d,m)==1 for d in D) for a in range(1,m,2))
            predicted=L*prod(p-len({d%p for d in D}) for p in P)
            need(actual==predicted,('joint ports',m,D));joint_tests+=1
    extra={m:wheel(m) for m in [290,5510,2310,30030,1160,1890]}
    P=primes_to(6100)
    first_plus=next(p for p in P if p%290==1)
    first_minus=next(p for p in P if p%290==289)
    pairs=[(p,p+2) for p in P if p%290==289 and p+2 in P]
    need([first_plus,first_minus,pairs[0]]==[1451,2029,(6089,6091)],'prime clocks')
    result={'parameter_range':[-99,100005,8],'generating_function_checks':12502,
       'small_parameters':small,'layer_basis':layer,
       'binary_line_patterns_checked':checked,'hourglass_hostile':hostile,
       'quadratic_cube':F,'quotient':[[list(k),list(v)] for k,v in sorted(Q.items())],
       'image_sizes':[8,len(set(F)),len({F[x] for x in F})],
       'second_image_fibres':sorted(Counter(F[x] for x in F).values()),
       'tournaments':trows,'class_descent_hostile':pair_failures[0],
       'merged_descent_hostile':merged_failures[0],
       'fractran':[str(f) for f in program],'source_clock_traces':traces,
       'joint_port_masks':joint_tests,
       'wheel_universe':'all even m=4..2000 plus recorded extras','wheels':extra,
       'prime_onsets':[first_plus,first_minus,list(pairs[0])],
       'eligible_prime_index':sum(gcd(p,290)==1 for p in P if p<=1451)}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
