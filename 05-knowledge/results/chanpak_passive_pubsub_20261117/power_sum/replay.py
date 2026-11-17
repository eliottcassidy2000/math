import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"observer"))
from replay import covers,slots,path_table,require
from math import prod

pairs=(
 (4287,((241,(126,72,122,192,170,138,176,7)),
        (257,(94,116,104,34,204,209,8,79)),
        (273,(178,72,163,245,198,210,133,2))),
       (241,273),(257,257),(20581887,131586,512)),
 (4803,((293,(48,113,187,147,96,200,141,19)),
        (301,(112,129,131,183,134,22,62,97)),
        (317,(190,72,242,196,170,74,17,98)),
        (325,(46,192,186,162,107,130,173,17))),
       (293,317,317),(301,301,325),(133559139123,58886554,4096))
)

def strong(a):
    for v in range(len(a)):
        t=1<<v
        while True:
            z=t
            for i in range(len(a)):
                if t>>i&1:z|=a[i]
            if z==t:break
            t=z
        if t!=(1<<len(a))-1:return False
    return True

def cyclic(l,B):
    m=len(B);q=l+m
    return tuple(((1<<q)-(1<<(l)))|((1<<l)-(1<<(v+1))) for v in range(l))+tuple(
           (a<<l)|(1<<q) for a in B)+((1<<l)-1,)

def main():
    for l in range(1,7):
        for B in [(0,), (2,0), (2,4,1)]:
            C=covers(B);h,p=slots(C,1),slots(C,2)
            a=cyclic(l,B)
            require(strong(a) and path_table(a)[-1]==(2**l-1)*(p-h)+2*h,l,B)
    for d,rows,lo,hi,expected in pairs:
        for h,a in rows:
            require(all(i!=j and ((a[i]>>j)^(a[j]>>i))&1 or i==j and not(a[i]>>i&1)
                        for i in range(8) for j in range(8)) and strong(a),a)
            c=covers(a)
            require((slots(c,1),slots(c,2))==(h,h+d),a)
        for j in range(1,len(lo)):
            require(sum(x**j for x in lo)==sum(x**j for x in hi),lo,hi)
        D=lambda a:prod(h+d for h in a)-prod(a)
        require(D(lo)==D(hi),lo,hi)
        got=(D(lo),2*prod(lo),2*(prod(hi)-prod(lo)))
        require(got==expected,got)
        print('verified common high topic, first low topic, gap:',got)

if __name__=='__main__':main()
