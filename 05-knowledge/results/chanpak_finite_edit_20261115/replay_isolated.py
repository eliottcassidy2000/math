from fractions import Fraction
from itertools import product
from random import Random
from replay import edit_pairs, hamilton, singleton_count

def literal(gaps):
    v=[0]
    for t in gaps:
        v.append(v[-1]+t+1)
    return edit_pairs(singleton_count(v[-1]+1),list(zip(v,v[1:])))

def transfer(gaps):
    h=1
    last=-1
    U,V=1,0
    for i,t in enumerate(gaps):
        if not t:
            continue
        a=2**t
        if last==-1:
            U,V=1,a*(2 if i else 1)
        elif i-last==2:
            U,V=U+2*V,2*a*(U+V)
        elif i-last==1:
            h*=U+V
            U,V=1,a
        else:
            raise ValueError('non-isolated zero gaps')
        last=i
    if last>=0:
        h*=U+(2 if last<len(gaps)-1 else 1)*V
    return h

def main():
    k=0
    for q in range(7):
        for gaps in product(range(4),repeat=q):
            if q+sum(gaps)>9 or any(x==y==0 for x,y in zip(gaps,gaps[1:])):
                continue
            if transfer(gaps)!=hamilton(literal(gaps)):
                raise ValueError(('transfer',gaps))
            k+=1
    c=(3239,2289,1547,1017,659,424,271,174)
    for x,v in enumerate(c,1):
        if v**20*(3*2**(x+1)+5)**13 <= 10_000**20 * 3**13:
            raise ValueError(('rational term bound',x))
    geom=Fraction(16,25)
    if geom**20*2**13<=1:
        raise ValueError('geometric bound')
    R=Fraction(sum(c),10_000)+geom**10/(1-geom)
    if R>=Fraction(995,1000):
        raise ValueError('Dirichlet threshold bound')
    print(f'{k} physical corridor controls; rational R < {float(R):.12f} < 0.995')

if __name__=='__main__':
    main()
