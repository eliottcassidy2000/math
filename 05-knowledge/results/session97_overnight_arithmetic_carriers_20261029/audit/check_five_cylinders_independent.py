"""Independent classifier audit, no producer imports."""
from fractions import Fraction as F
import json

def require(x,m):
    if not x:raise ArithmeticError(m)

def run(n):
    q=1;word=[];branch_sources=[]
    while q:
        t=3*n+1;a=0
        while t%2==0:t//=2;a+=1
        if a%3==0:branch_sources.append(n);q+=2
        else:q-=1
        word.append(a);n=t
    return n,word,branch_sources

def main():
    charts=[([1],3,4),([3,1,1,1],125,128),([3,3,1,1,1,1,1],1357,4096),
            ([3,1,3,1,1,1,1],2269,4096),([3,1,1,3,1,1,1],2877,4096)]
    # Direct residue search verifies exact word cylinders; no carry formula used.
    for w,r,m in charts:
        hits=[n for n in range(1,m,2) if run(n)[1]==w]
        require(hits==[r],'unique exact source word in residue period')
    small=[]
    for n in range(13,125,16):
        z,w,bs=run(n);bad=next(i for i,a in enumerate(w) if a not in (1,3))
        q=1
        for a in w[:bad]:
            require(q>0,'forbidden letter not reached before local stop');q+=2 if a==3 else -1
        require(q>0,'forbidden letter consumed while active');small.append([n,w[:bad+1]])
    require(F(3,2)*F(6,7)**3<1,'three branch exclusion')
    require(F(3,2)*F(117,448)<1,'high branch exclusion')
    checked=0
    for base in range(1,4096,2):
        for shift in [0,1,17,10001]:
            n=base+4096*shift;z,w,bs=run(n)
            match=[ww for ww,r,m in charts if n%m==r]
            require((z>n)==bool(match),'iff')
            require((z==n)==(n==1),'equality')
            require(not match or w==match[0],'word identity')
            if all(a in (1,3) for a in w):require(all(v>=125 for v in bs),'low branch exclusion')
            checked+=1
    print(json.dumps({'seven_local_obstructions':small,'hostile_source_count':checked,
                     'expansion_density':str(F(1,2)+F(1,64)+3*F(1,2048))},sort_keys=True))
if __name__=='__main__':main()
