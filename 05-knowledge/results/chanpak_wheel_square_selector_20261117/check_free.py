from math import isqrt
from itertools import product
from fractions import Fraction
from check import need, apply
from ray.replay import ends, check_graph, strong, prime

D=(54,32,2,119,6,20,55)
mask=26
def TT(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def apex(d,q):
    return tuple(row|(((q>>j&1)^1)<<len(d)) for j,row in enumerate(d))+(q,)
def make(u,w):
    s=bin(u)[2:];t=''.join(str(1-int(b)) for b in reversed(bin(w)[2:]))
    a=TT(len(s)); c=TT(len(t)); n=len(a)+len(D)
    x=((1<<len(D))-1)<<len(a); y=((1<<len(c))-1)<<n
    old=tuple(r|x|y for r in a)+tuple((r<<len(a))|y for r in D)+tuple(r<<n for r in c)
    tags=sum(int(b)<<j for j,b in enumerate(s))|(mask<<len(a))|sum(int(b)<<(n+j) for j,b in enumerate(t))
    return apex(old,tags)
def F(x):
    return apply(bin(x)[2:])[1]
def H(u,w):
    return 91+5*(F(u)+F(w)-2)+142*(u+w)+256*u*w+2*(w%2)
def main():
    check_graph(D)
    h,e=ends(D)
    _,s=ends(tuple(sum((r>>i&1)<<j for j,r in enumerate(D)) for i in range(len(D))))
    f=len(h)-1
    R=[0,0]; L=[0,0]
    for v in range(f+1):
        j=v.bit_count()%2
        R[j]+=h[f^v]*sum(c for i,c in enumerate(e[v]) if not mask>>i&1)
        L[j]+=h[f^v]*sum(c for i,c in enumerate(s[v]) if mask>>i&1)
    need((h[-1],sum(h[i]*h[f^i] for i in range(f+1)),sum(L),sum(R),L[0]-L[1],R[0]-R[1],ends(apex(D,mask))[0][-1])==(5,256,142,142,-2,0,91))
    for u,w in list(product((1,3,5,7),(2,4,6)))+[(81,2)]:
        for a,b in [(u,w),(w,u)]:
            graph=make(a,b);check_graph(graph);need(strong(graph))
            need(ends(graph)[0][-1]==H(a,b))
        need(H(w,u)==H(u,w)+2)
    need(H(81,2)==53549 and prime(53549) and prime(53551))
    results=[]
    for k in range(4,51,2):
        N=10**k+30011
        Y=1<<max(2,(N.bit_length()-1)//10)
        A=256*Y**2;B=142*Y+5*(Y//2-1);C=398*Y+5*(Y-1)
        p,q,E=Fraction(C,A),Fraction(B,A),223-Fraction(B*C,A)
        m=isqrt(4*N//A)
        while A*(m+p+q)**2/4+E<N:m+=1
        while m and A*(m-1+p+q)**2/4+E>=N:m-=1
        if m<10:continue
        G=lambda a: A*a*(m-a)+B*a+C*(m-a)+223
        l,r=0,int((m+q-p)/2)+1
        while r-l>1:
            a=(l+r)//2
            if G(a)>=N:r=a
            else:l=a
        opts=[a for a in (l,r,r-1,r+1) if 1<=a<m]
        a=min(opts,key=lambda a:abs(N-H(Y*a+1,Y*(m-a))))
        b=m-a; lo,hi=H(Y*a+1,Y*b),H(Y*b,Y*a+1)
        need(hi==lo+2 and abs(N-lo)**5 <= 1000**5*N*N)
        results.append(k)
    print(f'Full signed profile, 26 literal DPs including twin prime counts, {len(results)} integer scale controls pass.')
if __name__=='__main__':
    main()
