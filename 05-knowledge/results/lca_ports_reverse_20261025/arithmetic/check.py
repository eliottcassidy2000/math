"""Exact arithmetic-map controls; no native TM or infinite-orbit claim."""
from fractions import Fraction as Q
from itertools import product, combinations
import hashlib
import json


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def v2(n):
    if not n:
        raise ValueError('valuation of zero is infinite')
    n=abs(n)
    return (n & -n).bit_length()-1


def residue(x, modulus):
    x=Q(x)
    return x.numerator*pow(x.denominator,-1,modulus)%modulus


def parity_step(x, kind):
    if kind=='C':
        return (3*x+1)//2 if x%2 else x//2
    return (3*x-x%2)//2


def odd_step(x, kind):
    value=3*x+(1 if kind=='C' else -1)
    k=v2(value)
    return (value//(1<<k))*(3**(k-1) if kind=='F' else 1),k


def inverse_word(word, tail, kind):
    x=Q(tail)
    for k in reversed(word):
        x=((1<<k)*x-1)/3 if kind=='C' else Q(1,3)+Q(2,3)**k*x
    return x


def macro(a,b):
    if b%2:
        return None if a==0 else (a-1,(3*b+3)//2)
    return a+2,3*b//2+2


def translations():
    checked=0
    for a in range(9):
        for b in range(4,260):
            out=macro(a,b)
            if out is None:
                need(a==0 and b%2,'native halt predicate')
                continue
            ap,bp=out; y=b+4; z=b+3
            need(bp+4==3*y//2 and bp+3==(3*z+1)//2,'integer translations')
            need(ap==a+3*(z%2)-1,'credit cocycle')
            checked+=1
    state=(0,4)
    for _ in range(3):state=macro(*state)
    need(state==(6,23),'blank first odd point')
    need(macro(0,3) is None and odd_step(7,'F')==(15,2),
         'net positive odd-return credit still needs initial guard')
    for x in range(5,1000,2):
        for a in (1,3,6):
            target,k=odd_step(x,'F'); state=(a,x-4)
            for _ in range(k):
                state=macro(*state)
                need(state is not None,'internal odd-return guard')
            need(state==(a+2*k-3,target-4),'accelerated return agrees with macros')
            need(target>x,'strict positive integer growth')
            checked+=1
    return checked


def binary_cylinders():
    checked=0
    for depth in range(1,11):
        for kind in ('C','F'):
            residues=set()
            for bits in product((0,1),repeat=depth):
                x=Q(0)
                for bit in reversed(bits):
                    x=((2*x-1)/3 if bit else 2*x) if kind=='C' else (2*x+bit)/3
                r=residue(x,1<<depth); residues.add(r); y=r
                for bit in bits:
                    need(y%2==bit,'binary cylinder parity')
                    y=parity_step(y,kind)
                checked+=1
            need(len(residues)==1<<depth,'full binary cylinder bijection')
    return checked


def valuation_cylinders_and_lca():
    cylinders=0
    for size in range(1,5):
        for word in product(range(1,5),repeat=size):
            modulus=1<<(sum(word)+1)
            for kind in ('C','F'):
                r=residue(inverse_word(word,1,kind),modulus)
                for lift in (1,3):
                    x=r+lift*modulus
                    for k in word:
                        x,observed=odd_step(x,kind)
                        need(observed==k,'exact valuation cylinder')
                    cylinders+=1
    words=list(product(range(1,5),repeat=4))
    values={kind:[inverse_word(w,-1 if kind=='C' else 1,kind) for w in words]
            for kind in ('C','F')}
    pairs=0
    for i,j in combinations(range(len(words)),2):
        A=0
        for k,l in zip(words[i],words[j]):
            if k!=l:
                expected=A+min(k,l);break
            A+=k
        for kind in ('C','F'):
            difference=values[kind][i]-values[kind][j]
            need(v2(difference.numerator)-v2(difference.denominator)==expected,'weighted LCA')
        pairs+=1
    need(inverse_word((2,),Q(1),'C')==1,'Collatz periodic point')
    need(inverse_word((2,),Q(3,5),'F')==Q(3,5),'noninteger conjugate')
    need(residue(Q(1,3),32)==11 and residue(Q(-1,3),32)==21,'distinct singular branches')
    try:v2(0)
    except ValueError:pass
    else:raise RuntimeError('singular valuation improperly finite')
    return {'exact_lifts':cylinders,'weighted_pairs':pairs}


def fuel_and_extensions():
    fuel=0; family=0; extensions=0
    for x in range(3,3000,2):
        d=v2(x-1);y=x;run=0
        while True:
            nxt,k=odd_step(y,'F')
            if k!=1:break
            run+=1;y=nxt
        need(run==d-1,'marked ancestor fuel')
        fuel+=1
    samples=[]
    for R in range(5,198,6):
        x=((1<<(R+2))+7)//9
        need(9*x==(1<<(R+2))+7 and x%2,'family integrality')
        xp,k=odd_step(x,'F')
        need(v2(x-1)==1 and k==2 and xp==1+(1<<R),'reset family')
        need(v2(xp-1)==R,'unbounded regenerated depth')
        if R in (5,11,17):samples.append({'R':R,'x':x,'next':xp})
        for a in range(1,17):
            if R<a+2:continue
            credit=a+1;y=xp
            for _ in range(a+1):
                need(credit>0,'family prehalt guard')
                y,letter=odd_step(y,'F')
                need(letter==1,'depletion after reset');credit-=1
            need(credit==0,'family hits native barrier');family+=1
    for size in range(5):
        for word in product(range(1,5),repeat=size):
            a=6;valid=True
            for k in word:
                if a<=0:valid=False;break
                a+=2*k-3
            if not valid or a<=0:continue
            extended=word+(1,)*a
            modulus=1<<(sum(extended)+1)
            x=residue(inverse_word(extended,1,'F'),modulus)+modulus
            credit=6
            for k in extended:
                need(credit>0,'extension first halt')
                x,observed=odd_step(x,'F')
                need(observed==k,'integer extension itinerary')
                credit+=2*k-3
            need(credit==0,'integer extension halt');extensions+=1
    return {'fuel_sources':fuel,'family_credit_cases':family,
            'halting_extensions':extensions,'reset_samples':samples}


def quadratic_probability():
    # a+b*rho with rho^2=1-rho, exact rational arithmetic.
    def mul(x,y):
        a,b=x;c,d=y
        return a*c+b*d,a*d+b*c-b*d
    def add(x,y):return x[0]+y[0],x[1]+y[1]
    def scale(t,x):return t*x[0],t*x[1]
    powers=[(Q(1),Q(0))]
    for _ in range(160):powers.append(mul(powers[-1],(Q(0),Q(1))))
    need(powers[6]==(5,-8),'golden probability at credit6')
    for a in range(64):
        need(scale(Q(1,2),add(powers[a],powers[a+3]))==powers[a+1],
             'harmonic equation including ruin boundary')
    checks=0
    for start in range(9):
        bank={start:Q(1)};ruin=Q(0)
        for _ in range(64):
            nxt={}
            for a,p in bank.items():
                for b in (a-1,a+2):
                    if b==-1:ruin+=p/2
                    else:nxt[b]=nxt.get(b,Q(0))+p/2
            bank=nxt; value=(ruin,Q(0))
            for a,p in bank.items():value=add(value,scale(p,powers[a+1]))
            need(value==powers[start+1],'exact stopped martingale balance')
            checks+=1
    return {'harmonic_cases':64,'finite_horizon_balances':checks,'rho_power6':[5,-8]}


def blank_prefix():
    x=27;a=6;steps=0;macro_steps=0;min_credit=a;min_margin=a-v2(x-1);max_letter=0
    for _ in range(4096):
        if a==0:break
        xp,k=odd_step(x,'F'); direct=x;credit=a
        for j in range(k):
            if direct%2:
                need(credit>=1,'direct macro guard');credit-=1
            else:credit+=2
            direct=3*direct//2
        need((direct,credit)==(xp,a+2*k-3),'blank direct return comparison')
        x,a=xp,credit;steps+=1;macro_steps+=k;min_credit=min(min_credit,a)
        min_margin=min(min_margin,a-v2(x-1));max_letter=max(max_letter,k)
    return {'odd_returns':steps,'macro_steps_after_first_odd':macro_steps,
            'halt_seen':a==0,'final_credit':a,'final_scalar_bits':x.bit_length(),
            'minimum_credit':min_credit,'minimum_LCA_margin':min_margin,
            'maximum_letter':max_letter,
            'final_scalar_sha256':hashlib.sha256(x.to_bytes((x.bit_length()+7)//8,'little')).hexdigest()}


def main():
    result={'status':'PASS','translation_checks':translations(),
            'binary_cylinders':binary_cylinders(),
            'valuation_controls':valuation_cylinders_and_lca(),
            'fuel_controls':fuel_and_extensions(),'haar_controls':quadratic_probability(),
            'finite_blank_prefix':blank_prefix(),
            'scope':'abstract arithmetic maps; native macro interface cited; no infinite blank-orbit claim'}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
