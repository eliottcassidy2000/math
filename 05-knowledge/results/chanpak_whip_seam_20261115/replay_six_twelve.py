# Literal endpoint-subset DP, independent of the SCC search program.
def tables(r):
    n=len(r);f=(1<<n)-1
    a=[[0]*n for _ in range(f+1)];h=[1]+[0]*f
    for s in range(1,f+1):
        w=s
        while w:
            b=w&-w;w^=b;i=b.bit_length()-1;z=s^b
            if not z:a[s][i]=1
            else:
                t=z&r[i]
                while t:
                    j=t&-t;t^=j;a[s][i]+=a[z][j.bit_length()-1]
        h[s]=sum(a[s])
    return a,h
def add(x,y):
    a,m=x;b,k=y
    n=len(a);T=((1<<len(b))-1)<<n
    return tuple(v|T for v in a)+tuple(v<<n for v in b),m|(k<<n)
def roots(seq):
    x=(),0
    for y in seq:x=add(x,y)
    return x
def bits(w):return [((0,),int(b)) for b in w]
def state(x):
    r,m=x;n=len(r);a,h=tables(r);f=len(h)-1
    S=d=0
    for s in range(1,f+1):
        z=sum(a[s][i] for i in range(n) if m>>i&1)*h[f^s]
        S+=z;d+=z if s.bit_count()&1 else -z
    out=tuple(v|(0 if m>>i&1 else (1<<n)) for i,v in enumerate(r))+(m,)
    return h[-1],tables(out)[1][-1],S,d
def p(t):
    if not t:return 0
    u=t//2
    return p(u)+(u+(u%2) if t%2==0 else 0)
F=(36,1,2,7,15,30)
B=(8,1,3,6)
C=(4,1,2)
bank=[
 [(F,29)]+bits('100'),
 bits('1')+[(B,6)]+bits('1001')+[(C,7)],
 [(F,16)]+bits('1100'),
 bits('101')+[(F,48)],
 [(C,1),(B,6)]+bits('100'),
 bits('1')+[(C,3),(B,2)]+bits('0'),
]
def main():
    for i,seq in enumerate(bank):
        h,n,s,d=state(roots(seq))
        if (h,n,s)!=(15,979+2*i,9609 if i==1 else 1188):raise AssertionError((i,h,n,s,d))
        print('literal alignment base',i,(h,n,s,d))
        if not seq[0][1]:raise AssertionError('missing entry to first SCC')
    for g in range(4,17):
        A=1203<<g;S=2*(A-15);c={0:6,1:12,2:-6,3:-12}[g%4]
        a=(S-c)//15
        if S-15*a!=c or a%2:raise AssertionError(g)
        x=2*a
        def N(t):return 979+A+S*t+15*p(t)
        if N(x+1)-N(x)!=c:raise AssertionError(g)
        W={N(t)+2*i for t in (x,x+1) for i in range(6)}
        if len(W)!=6+abs(c)//2 or sorted(W)!=list(range(min(W),max(W)+1,2)):
            raise AssertionError(g)
        if g<=7:print('fold height',g,'addresses',x,x+1,'window',min(W),max(W),len(W))
    print('13 exact fold-clock controls; literal base counts obtained independently')
if __name__=='__main__':main()
