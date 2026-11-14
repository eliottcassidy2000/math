"""Compute exact optional-two-slot response maps of the published full-digit bank."""
from replay_corridor import tt, star, components, induce, endtable, need
from pathlib import Path
from json import loads
from fractions import Fraction

def q2(T):
    H=[1]+[sum(a) for a in endtable(T)[1:]]
    return sum(a*b for a,b in zip(H,H[::-1]))

def v2(q):
    def i2(n):
        n=abs(n)
        return (n&-n).bit_length()-1
    return i2(q.numerator)-i2(q.denominator)

def signature(row):
    I,ic,W,wc=tuple(row['I']),row['ic'],tuple(row['W']),row['wc']
    y=[]
    for m in (1,2,3):
        T=star(I,tt(m),W,ic,wc)
        C=[part for part in components(T) if len(part)>1]
        need(len(C)==1 and len(I) in C[0])
        y.append(q2(induce(T,C[0])))
    c=Fraction(y[2]-3*y[1]+2*y[0],6); b=(y[1]-y[0]-6*c)/2
    a=y[0]-2*b-3*c
    r=row['shift']
    return a,b/2**r,c/3**r

def main():
    bank=loads(Path(__file__).with_name('corridor_witnesses.json').read_text())
    sig=[signature(row) for row in bank['rails'].values()]
    scaled={tuple(x/Fraction(2)**v2(next(c for c in v if c)) for x in v) for v in sig}
    need(len(sig)==len(scaled)==99)
    print('99 exact optional-two-slot maps distinct even up to dyadic rescaling')
if __name__=='__main__':
    main()
