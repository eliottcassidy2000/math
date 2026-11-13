"""Literal finite controls; all-length claims are proved in the note."""
from itertools import product

def ck(test, receipt):
    if not test:
        raise AssertionError(receipt)

def step(n):
    y=3*n+1
    a=(abs(y)&-abs(y)).bit_length()-1
    return y//2**a,a

def walk(n, steps):
    out=[]
    for _ in range(steps):
        n,a=step(n)
        out.append(a)
    return n,tuple(out)

def carry(word):
    b=s=0
    for a in word:
        b=3*b+2**s
        s+=a
    return b,s

def run():
    a_bank=(-5,-7)
    b_bank=(-17,-25,-37,-55,-41,-61,-91)
    tested=0
    for anchors in (a_bank,b_bank):
        bank=[]
        for z in anchors:
            x,w=walk(z,len(anchors))
            ck(x==z, (x,z))
            bank.append(w)
        for k in range(1,5):
            for ids in product(range(len(anchors)),repeat=k):
                w=tuple(a for j in ids for a in bank[j])
                B,s=carry(w)
                gap=2**s-3**len(w)
                ck((B%gap==0)==(len(set(ids))==1),ids)
                if B%gap==0:
                    ck(B//gap==anchors[ids[0]],ids)
                tested+=1
    for k,M in product(range(1,11),range(1,8)):
        for z,word in ((-5,(1,2)),(-7,(2,1))):
            x,w=walk(2*8**k*M+z,2*k)
            ck(x==2*9**k*M+z and w==word*k,(k,M,z))
    print("finite two- and seven-phase concatenations:",tested)
    print("positive native shadow pairs:",70)

if __name__=='__main__':
    run()
