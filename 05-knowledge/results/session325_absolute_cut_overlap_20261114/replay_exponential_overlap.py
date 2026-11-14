"""Finite hostile controls for the effective dyadic overlap certificate."""
from itertools import product
from random import Random

def partitions(n):
    if n==0:
        yield ()
        return
    for p in partitions(n-1):
        yield p+((n-1,),)
        for i in range(len(p)):
            yield p[:i]+(p[i]+(n-1,),)+p[i+1:]

def templates(cs):
    out=[]
    for p in partitions(len(cs)):
        rows=[]
        for B in p:
            if len(B)<2:
                break
            G=sum(abs(cs[i]) for i in B).bit_length()
            bound=(len(B)-1)*(G-1)
            valid=[]
            for ds in product(range(bound+1),repeat=len(B)):
                if min(ds)==0 and sum(cs[i]*(1<<d) for i,d in zip(B,ds))==0:
                    valid.append(ds)
            if not valid:
                break
            rows.append((B,set(valid)))
        else:
            out.append(rows)
    return out

def pass_template(ts, es):
    for rows in ts:
        if all(tuple(es[i]-min(es[j] for j in B) for i in B) in opts for B,opts in rows):
            return True
    return False

def check(x):
    if not x:
        raise RuntimeError('bad overlap test')

def main():
    rng=Random(325)
    packets=[(1,-1,3,-3),(1,1,-2),(-2,1,1,1)]
    packets += [tuple(rng.choice((-3,-2,-1,1,2,3)) for _ in range(3)) for _ in range(20)]
    for cs in packets:
        T=templates(cs)
        for es in product(range(7),repeat=len(cs)):
            check(pass_template(T,es)==(sum(c*(1<<e) for c,e in zip(cs,es))==0))
    check({(i,j) for i,j in product(range(8),repeat=2) if 3*(1<<i)-(1<<j)==2}=={(0,0),(1,2)})
    print('dyadic overlap certificate: direct and partition decoders agree')
if __name__=='__main__':
    main()
