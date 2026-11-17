from itertools import product

def need(ok):
    if not ok:
        raise AssertionError("gap coefficient replay failed")

def F(w):
    return 1+sum((1<<(j-i-1))+(j==i+1)
        for i in range(len(w)) for j in range(i+1,len(w))
        if w[i]=='1' and w[j]=='0')

def L(w): return int(w or '0',2)
def R(w): return sum((b=='0')<<i for i,b in enumerate(w))
def G(w): return L(w)-F(w)+1-(L(w)&1)

def case(P,Q,S,s,t):
    a,b,r,p,l=L(P),R(S),L(Q),R(Q),len(Q)
    if (s,t)==('0','0'):
        return a*(1<<l)*(b+1), -a*L(Q[::-1]), r*(b+1)
    if (s,t)==('0','1'):
        return a*b*(1<<l), a*(p+1), b*(r+1)
    return (a+1)*(b+1)*(1<<l),(a+1)*(p-(1<<l)),(b+1)*(r-(1<<l))

def literal(w):
    n=len(w); a=[sum(1<<j for j in range(i+1,n))+((b=='0')<<n)
                for i,b in enumerate(w)]+[L(w[::-1])]
    n+=1; f=(1<<n)-1
    d=[[0]*n for _ in range(f+1)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,f+1):
        for i,r in enumerate(a):
            z=r&~s
            while z:
                b=z&-z;z^=b
                d[s|b][b.bit_length()-1]+=d[s][i]
    return sum(d[-1])

def run():
    B=['']+[''.join(w) for k in range(1,4) for w in product('01',repeat=k)]
    for w in B:
        need(G(w)==sum((1<<(j-i-1))+(j==i+1)
                for i in range(len(w)) for j in range(i+1,len(w)) if w[i]==w[j]=='1'))
    for P,Q,S in product(B,repeat=3):
        for s,t in [('0','0'),('0','1'),('1','0')]:
            a,b,d=case(P,Q,S,s,t)
            C=F(P+s+Q+t+S)-4*a-2*b-2*d
            for m,n in [(1,2),(2,1),(3,2)]:
                x,y=1<<m,1<<n
                need(F(P+s*m+Q+t*n+S)==a*x*y+b*x+d*y+C)
            if s==t=='0':
                need(C==F(S)-G(P)-G(Q))
    for w in ('10011','10111010','101011100','11000100100','0011000110'):
        need(literal(w)==F(w))
    for r in range(8):
        for t in range(8):
            a,b=1<<r,7*(1<<t)-1
            need(F(f'{b:b}')-G(f'{a:b}')==5)
            a,b=7*(1<<r),(1<<t)-1
            need(F(f'{b:b}')-G(f'{a:b}')==-5)
    print('Exact same-symbol constants, mixed-gap coefficients and literal checks: PASS')
if __name__=='__main__':
    run()
