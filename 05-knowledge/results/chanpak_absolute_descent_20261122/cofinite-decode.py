"""Exact rooted-unique-start count decoder for the signed wrap grammar."""
from argparse import ArgumentParser

def p(y):
    s=0
    while y:
        y,r=divmod(y,2)
        if not r:s+=y+(y&1)
    return s

def rooted_histories(N):
    """Yield every strong rooted wrap history with absolute height N."""
    if N<1 or N%2==0:return
    if N==1:
        yield ()
        return
    def descend(h,c,e,word):
        if h==N:
            yield word
            return
        if word:
            a,b,f=h+c-e,h+c,h-e
            if a<=N:yield from descend(a,b,f,word+((0,0),))
        # If y >= 2^(length-1), h' >= c(y+1).
        end=N//c
        for length in range(1,end.bit_length()+1):
            B=1<<length
            for y in range(B//2,min(B,end)):
                h1=(1+p(y))*h+(y+1)*c-(not y%2)*e
                if h1>N:continue
                reverse=int(f"{y:0{length}b}"[::-1],2)
                c1=B*c+(B-reverse)*h
                yield from descend(h1,c1,h,word+((length,y),))
    yield from descend(1,1,1,())

def graph(history):
    """Build each literal graph from the selected exponent/contact word."""
    rows=[0];root=0
    for length,y in history:
        n=len(rows);t=n+length
        mask=((1<<length)-1)<<n
        rows=[r|mask|((1<<t) if i!=root else 0) for i,r in enumerate(rows)]
        rows.extend((mask&~((1<<(i+1))-1)) | ((1<<t) if y>>(i-n)&1 else 0)
                    for i in range(n,t))
        r=(1<<root)|sum(1<<(n+i) for i in range(length) if not y>>i&1)
        rows.append(r);root=t
    return tuple(rows),root

if __name__=="__main__":
    args=ArgumentParser(description=__doc__)
    args.add_argument("target",type=int)
    args.add_argument("--all",action="store_true")
    a=args.parse_args()
    work=rooted_histories(a.target)
    for w in work:
        print("history",w,"graph",graph(w))
        if not a.all:break
    else:
        if not a.all:print("no solution in normalized uniquely-started strong grammar")
