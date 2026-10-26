"""Normal-basis marking and exact modular/valuation cycle gates; stdlib only."""
from collections import Counter
from fractions import Fraction
from itertools import product
import json

def need(b,s):
    if not b:raise RuntimeError(s)
def polyrem(a,b):
    while a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
    return a
def irreducible(f,m):
    return all(polyrem(f,g) for d in range(1,m//2+1) for g in range((1<<d)|1,1<<(d+1),2))
def mul(a,b,f):
    x=0
    while b:
        if b&1:x^=a
        a<<=1;b>>=1
    return polyrem(x,f)
def rank(vs):
    piv={}
    for x in vs:
        while x:
            b=x.bit_length()-1
            if b not in piv:piv[b]=x;break
            x^=piv[b]
    return len(piv)
def rot(w,m):return ((w<<1)|(w>>(m-1)))&((1<<m)-1)
def degree(x,m,f):
    y=x
    for d in range(1,m+1):
        y=mul(y,y,f)
        if y==x:return d
    raise RuntimeError('field period')

def field_checks():
    out=[]
    for m in [3,5,7,9]:
        f=next(f for f in range((1<<m)|1,1<<(m+1),2) if irreducible(f,m))
        for alpha in range(1,1<<m):
            basis=[alpha]
            for _ in range(m-1):basis.append(mul(basis[-1],basis[-1],f))
            if rank(basis)==m:break
        coords=[]
        for w in range(1<<m):
            x=0
            for i,b in enumerate(basis):
                if w>>i&1:x^=b
            coords.append(x)
        need(len(set(coords))==1<<m,'normal basis bijection')
        need(coords[-1]==1,'normal trace equals one')
        hist=Counter()
        for w,x in enumerate(coords):
            need(mul(x,x,f)==coords[rot(w,m)],'Frobenius is decoration rotation')
            d=degree(x,m,f);hist[d]+=1
            # A basis-vector port forces the complete Frobenius orbit.
            pairs=[];a=x;b=alpha
            for _ in range(m):pairs.append((a,b));a=mul(a,a,f);b=mul(b,b,f)
            need(len(set(pairs))==m and (a,b)==(x,alpha),'one marked basis port kills cyclic stabilizer')
            need(m%d==0,'degree divides ambient order')
        inv={x:w for w,x in enumerate(coords)}
        def reflect(x):
            w=inv[x]
            v=sum(1<<((-i)%m) for i in range(m) if w>>i&1)
            return coords[v]
        hostile=next((x,y) for x in range(1<<m) for y in range(1<<m)
                     if reflect(mul(x,y,f))!=mul(reflect(x),reflect(y),f))
        out.append({'m':m,'irreducible_polynomial_bits':f,'normal_element':alpha,'basis':basis,
                    'element_degree_counts':dict(sorted(hist.items())),
                    'necklace_classes':sum(count//d for d,count in hist.items()),
                    'reflection_not_multiplicative':list(hostile)})
    return out

def valuation(x):
    need(x!=0,'nonzero valuation')
    return (abs(x)&-abs(x)).bit_length()-1
def syr(n,sign):
    x=3*n+sign
    return x//(1<<valuation(x))
def cycle_inventory(f):
    done=set();cycles=[]
    for x in range(len(f)):
        trail=[];where={}
        while x not in done and x not in where:
            where[x]=len(trail);trail.append(x);x=f[x]
        if x in where:cycles.append(len(trail)-where[x])
        done.update(trail)
    return sorted(cycles)

def modular_checks():
    crt=0;edges=0;samples=[]
    for m in range(3,102,2):
        q=pow(4,-1,m);k=(m-1)//2
        need(k*k%m==q and (k+1)**2%m==q and k*(k+1)%m==(-q)%m,'central four products')
        need((3*q-1)%m==(-q)%m and (-3*q+1)%m==q,'alternating affine centre cycle')
        for a in range(7):
            plus=[(3*x+1)*pow(pow(2,a,m),-1,m)%m for x in range(m)]
            minus=[(3*x-1)*pow(pow(2,a,m),-1,m)%m for x in range(m)]
            need(all(minus[x]==(-plus[-x%m])%m for x in range(m)),'sign conjugacy at fixed clock')
            need(cycle_inventory(plus)==cycle_inventory(minus),'equal modular sign cycle inventory')
            edges+=m
            if a==0 and m in [5,7,11,13]:samples.append({'modulus':m,'both_signs_cycle_lengths':cycle_inventory(plus)})
        for sign in [-1,1]:
            for a in range(1,9):
                dyadic=1<<(a+1);u=((1<<a)-sign)*pow(3,-1,dyadic)%dyadic
                for r in range(m):
                    n=u+dyadic*((r-u)*pow(dyadic,-1,m)%m)
                    for j in [0,1,3]:
                        z=n+j*m*dyadic
                        need(z>0 and z%2==1 and z%m==r and valuation(3*z+sign)==a,'every clock in every odd residue')
                        crt+=1
    return {'moduli':'odd 3..101','fixed_clock_edges':edges,'crt_lifts':crt,
            'central_identity_same_cycle_inventory_examples':samples,
            'no_residue_map_controls':{'plus_mod5':[[n,syr(n,1)] for n in [1,11]],
                                       'minus_mod5':[[n,syr(n,-1)] for n in [1,31]]}}

def compositions(total,length):
    if length==1:yield (total,);return
    for a in range(1,total-length+2):
        for tail in compositions(total-a,length-1):yield (a,)+tail
def rotations(w):return [w[i:]+w[:i] for i in range(len(w))]
def numerator(w):
    c=0;prefix=0
    for j,a in enumerate(w):c+=3**(len(w)-1-j)*(1<<prefix);prefix+=a
    return c
def cycle_checks():
    words=0;necklaces=0;cycles={-1:set(),1:set()}
    for length in range(1,9):
        for total in range(length,19):
            for w in compositions(total,length):
                words+=1;rots=rotations(w)
                if w!=min(rots) or len(set(rots))!=length:continue
                necklaces+=1;c=numerator(w);d=(1<<total)-3**length
                if c%d:continue
                for sign in [-1,1]:
                    n=sign*c//d
                    if n<=0:continue
                    points=[];z=n
                    for a in w:
                        points.append(z)
                        need(valuation(3*z+sign)==a,'clock gate from odd integer divisor')
                        z=syr(z,sign)
                    need(z==n and len(set(points))==length,'exact primitive positive cycle')
                    cycles[sign].add(tuple(min(rotations(tuple(points)))))
    need(cycles[1]=={(1,)} and cycles[-1]=={(1,),(5,7),(17,25,37,55,41,61,91)},'bounded word census only')
    a=(1,1,1,2,1,1,4);b=(1,1,1,1,1,2,4);d=3**7-2**11
    need(sorted(a)==sorted(b) and numerator(a)%d==0 and numerator(b)%d!=0,'ordered carry separates equal clock multisets')
    return {'maximum_odd_period':8,'maximum_total_valuation':18,'ordered_words':words,
            'primitive_rotation_classes':necklaces,'positive_cycles_by_sign':{str(k):sorted(v) for k,v in cycles.items()},
            'same_multiset_hostile':{'words':[a,b],'denominator':d,'numerators':[numerator(a),numerator(b)],
                                      'minus_start_values':[str(Fraction(numerator(a),d)),str(Fraction(numerator(b),d))]}}

if __name__=='__main__':
    print(json.dumps({'fields':field_checks(),'modular':modular_checks(),'cycle_gate':cycle_checks()},sort_keys=True,indent=2))
