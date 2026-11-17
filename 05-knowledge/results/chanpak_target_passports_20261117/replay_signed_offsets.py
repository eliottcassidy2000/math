from math import isqrt

def need(x):
    if not x:
        raise AssertionError('variable-offset replay failed')

TILES=(
    (5,5,(2,4,9,3)),
    (19,8,(40,1,3,6,15,30)),
    (223,12,(2,40,11,49,71,213,143,31)),
    (599,13,(96,21,25,51,225,198,142,15)),
)
A=(4,5,72,19,39,79,155,63)

def table(a):
    n=len(a); f=(1<<n)-1
    need(all(0<=r<=f and not(r>>i&1) for i,r in enumerate(a)))
    need(all((r>>j&1)+(a[j]>>i&1)==1 for i,r in enumerate(a) for j in range(i)))
    d=[[0]*n for _ in range(f+1)];h=[1]+[0]*f
    for i in range(n):
        d[1<<i][i]=1
    for s in range(1,f+1):
        for i,r in enumerate(a):
            v=d[s][i]
            h[s]+=v
            if not v:
                continue
            z=r&~s
            while z:
                b=z&-z;z^=b
                d[s|b][b.bit_length()-1]+=v
    return d,h

def run():
    for u,e,a in TILES:
        _,h=table(a); f=len(h)-1
        need((h[-1],sum(z*h[f^s] for s,z in enumerate(h)))==(u,1<<e))
    for m,a in ((1,A),(2,(260,261,72,275,295,335,411,319,4))):
        d,h=table(a)
        f=len(h)-1
        b=sum(d[s][5]*h[f^s] for s in range(1,len(h)) if s>>5&1)
        need(d[-1][5]>0 and (h[-1],b)==(64*(1<<m)+1,96*(1<<m)-2))
    for u,k,p in ((5,12,4091),(19,15,32749),(599,28,268434857),(95,24,16777121)):
        need(p==(1<<k)-u and all(p%d for d in range(2,isqrt(p)+1)))
        print('prime quotient',p,'u',u,'height',k)
    for u,k,p in ((223,20,1048799),(599,27,134218327),(95,25,33554527),(25,20,1048601)):
        need(p==(1<<k)+u and all(p%d for d in range(2,isqrt(p)+1)))
        print('prime quotient',p,'u',u,'positive height',k)
    print('Exact tile interfaces and variable offset checks: PASS')
if __name__=='__main__':
    run()
