"""Exact complementary-carry placement; no floating point or external packages."""
from math import isqrt

def F(x):
    if x<0: raise ValueError('nonnegative address required')
    z,p,old=1,0,0
    for b in bin(x)[2:]:
        bit=int(b)
        if not bit:z+=old+p
        old,p=2*old+bit,bit
    return z

def count(a,b):
    if min(a,b)<1:raise ValueError('positive addresses required')
    if (a+b)%2!=1:raise ValueError('opposite address parities required')
    return 256*a*b+142*(a+b)+81+5*(F(a)+F(b))

def ceilsqrt(n):
    x=isqrt(n);return x+(x*x<n)

def scale(N):
    if N<1:raise ValueError('positive target required')
    x=isqrt(N//64)
    return x+(64*x*x<N)

def bracket(N,L=None):
    if N<2**18:raise ValueError('proved universal threshold is 2^18')
    M=scale(N)
    if L is None:
        L=1
        while (2*L)**3<=M:L*=2
    if L<2 or L&(L-1):raise ValueError('dyadic L>=2 required')
    m=((M+1+L-1)//L)*L-1
    delta=m-M
    d=ceilsqrt(2*m*(delta+1))
    right=(m-1)//2;left=right-d
    if left<1:raise ValueError('adaptive window leaves positive cone')
    hl=count(left,m-left);hr=count(right,m-right);evaluations=2
    if not hl<N<=hr:raise ArithmeticError('proved endpoint bracket failed')
    while right-left>1:
        mid=(left+right)//2;hm=count(mid,m-mid);evaluations+=1
        if hm<N:left,hl=mid,hm
        else:right,hr=mid,hm
    pairs=[]
    for a,h in [(left,hl),(right,hr)]:
        pairs.append({'a':a,'b':m-a,'u':a if a%2 else m-a,'w':m-a if a%2 else a,'lower':h,'upper':h+2})
    best=min(pairs,key=lambda p:abs(p['lower']-N))
    return {'target':N,'M':M,'L':L,'m':m,'delta':delta,'window_halfwidth':d,
            'pairs':pairs,'error':abs(best['lower']-N),'best':best,'evaluations':evaluations,
            'bracket_width':hr-hl}

if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser();p.add_argument('target',type=int);p.add_argument('--L',type=int)
    a=p.parse_args();print(json.dumps(bracket(a.target,a.L),indent=2))
