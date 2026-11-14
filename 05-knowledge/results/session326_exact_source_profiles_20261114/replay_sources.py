"""Literal certificates for 99*2^m + [13,15,...,35] and the independent overlaps."""
from itertools import product

def need(q):
    if not q:
        raise RuntimeError('incorrect corridor certificate')

def tt(n): return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))
def ordinal(A,B):return tuple(x|(((1<<len(B))-1)<<len(A)) for x in A)+tuple(x<<len(A) for x in B)
def add(A,mask): return tuple(x|((1<<len(A)) if not mask>>i&1 else 0) for i,x in enumerate(A))+(mask,)
def star(I,U,W,ic,wc):return add(ordinal(ordinal(I,U),W),ic|(wc<<(len(I)+len(U))))
def reach(T,i):
    seen={i}; todo=[i]
    while todo:
        a=todo.pop()
        for j in range(len(T)):
            if T[a]>>j&1 and j not in seen:
                seen.add(j); todo.append(j)
    return seen
def components(T):
    R=[reach(T,i) for i in range(len(T))]
    return {tuple(j for j in sorted(ri) if i in R[j]) for i,ri in enumerate(R)}
def induce(T,vs):return tuple(sum(1<<k for k,j in enumerate(vs) if T[i]>>j&1) for i in vs)
def validate(T):
    need(all(not (r>>i&1) and r>>len(T)==0 for i,r in enumerate(T)))
    need(all((T[i]>>j&1)+(T[j]>>i&1)==1 for i in range(len(T)) for j in range(i+1,len(T))))
def endtable(T):
    n=len(T); M=(1<<n)-1; E=[[0]*n for _ in range(M+1)]
    for i in range(n): E[1<<i][i]=1
    for S in range(1,M+1):
        for i,a in enumerate(E[S]):
            if a:
                r=T[i]&~S
                while r:
                    b=r&-r; r-=b; E[S|b][b.bit_length()-1]+=a
    return E
def profile(T,mask):
    M=(1<<len(T))-1
    e=endtable(T); p=endtable(tuple(M^r^(1<<i) for i,r in enumerate(T)))
    h=[1]+[sum(r) for r in e[1:]]
    W=[sum(k for i,k in enumerate(r) if mask>>i&1) for r in p]
    F=[sum(k for i,k in enumerate(r) if mask>>i&1) for r in e]
    E=[sum(k for i,k in enumerate(r) if not(mask>>i&1)) for r in e]
    return (h[M],
            sum(W[s]*h[M^s] for s in range(1,M+1)),
            sum(F[s]*W[M^s] for s in range(1,M)),
            h[M]+sum(E[s]*h[M^s] for s in range(1,M+1)),
            sum(endtable(add(T,mask))[-1]))
def count(T): return sum(endtable(T)[-1])
def seam(A,v,j):
    n=len(A); B=tt(j+1)
    return (tuple(r|(((1<<(j+1))-1)<<n if i!=v else ((1<<j)-1)<<n)
                  for i,r in enumerate(A))+
            tuple(r<<n|(1<<v if i==j else 0) for i,r in enumerate(B)))
def busch(m):
    P=tt(m+4); bits=(1,0,1,0)+(1,)*m+(0,)
    U=add(ordinal(P,(0,)),sum(x<<i for i,x in enumerate(bits)))
    A=add(U,1<<(len(U)-1))
    return seam(A,len(A)-1,1)


from math import gcd,isqrt
from fractions import Fraction as F
from pathlib import Path
import json
MINUS='10101010100';PLUS='10111011100'
def word(k,sign):return tuple(map(int,'1'*k+(PLUS if sign=='+' else MINUS)+'0'*k))
def rows(w):
    n=len(w)
    return tuple(sum(1<<j for j in range(i+1,n))|((1-w[i])<<n)for i in range(n))+(sum(a<<i for i,a in enumerate(w)),)
def P(w):
    a=[0]*(len(w)-1)
    for i in range(len(w)):
        for j in range(i+1,len(w)):
            if w[i] and not w[j]:a[j-i-1]+=1
    return a

def boundary(w):return [1-w[i]+w[-1-i] for i in range(len(w))]
def evaluate(a,t):return sum(c*t**i for i,c in enumerate(a))
def slot(w,z):
    n=len(w);return z**(n+1)+sum(a*z**(n-d-1)*((z+1)**d+(z-1)**d)for d,a in enumerate(P(w)))
def closed(k,sign):
    h=2048*4**k-1026*2**k+(131 if sign=='+' else 119)
    q=354294*9**k-242340*6**k+50336*4**k-2052*2**k+2
    return h,q

def hq(T):
    H=[1]+[sum(a)for a in endtable(T)[1:]]
    return H[-1],sum(a*b for a,b in zip(H,H[::-1]))
def catalog(cap):
    ans=[];k=0
    while closed(k,'+')[1]<=cap:
        for sign in ('-','+'):
            if sign=='-' and k%10 in (1,4):continue
            h,q=closed(k,sign);ans.append((q,h,k,sign))
        k+=1
    return ans

def inverse(N):
    if N<=0 or N%2==0:return None
    cap=N//3168;letters=catalog(cap)
    def dfs(start,h,q,parts):
        c=(N*pow(h,-1,198))%198
        need(1<=c<=197 and c%2)
        rem=N-c*h;den=99*q
        if rem>0 and rem%den==0:
            z=rem//den
            if z>=32 and z&(z-1)==0:return dict(c=c,t=z.bit_length()-1,h=h,q=q,letters=parts)
        for i in range(start,len(letters)):
            a,b,k,sign=letters[i]
            if q*a>cap:break
            r=dfs(i,h*b,q*a,parts+[(k,sign)])
            if r is not None:return r
        return None
    return dfs(0,1,1,[])

def prime(n):return n>=2 and all(n%p for p in range(2,isqrt(n)+1))

def main_sources():
    small=[]
    for k in range(3):
        for sign in ('-','+'):
            w=word(k,sign);T=rows(w);validate(T)
            need(len(components(T))==1)
            got=hq(T);need(got==closed(k,sign))
            small.append(dict(k=k,sign=sign,order=len(T),H=got[0],Q2=got[1]))
    a=P(word(0,'+'));b=P(word(0,'-'))
    need([x-y for x,y in zip(a,b)]==[-2,2,-1,0,-1,1,0,0,0,0])
    need(boundary(word(0,'+'))==boundary(word(0,'-')))
    B=boundary(word(0,'+'))
    need((evaluate(B,2),2**11*evaluate(B,F(3,2)),2**11*evaluate(B,F(1,2)))==(3070,466248,2048))
    formula=0
    for k in range(41):
        for sign in ('-','+'):
            w=word(k,sign);h,q=closed(k,sign)
            need((slot(w,1),slot(w,2))==(h,q))
            need(gcd(h,99)==1 if sign=='+' else ((gcd(h,99)==1)==(k%10 not in (1,4))))
            need(q>=4**(k+1));formula+=1
        need(slot(word(k,'+'),3)-slot(word(k,'-'),3)==-72900*9**k)
        need(gcd(closed(k,'+')[0],closed(k,'-')[0])==1)
    need(all(closed(k+1,'+')[1]>=4*closed(k,'+')[1]for k in range(41)))
    examples=[(0,17,507659921,None),(0,41,507687593,72526799),(1,19,6120213613,None),(1,127,6120890881,874412983)]
    bank=json.loads(Path(__file__).with_name('selected_templates.json').read_text());out=[]
    for k,c,N,p in examples:
        h,q=closed(k,'+');need(N==99*32*q+c*h)
        need(prime(N) if p is None else (N==7*p and prime(p)))
        ans=inverse(N);need(ans is not None and N==99*2**ans['t']*ans['q']+ans['c']*ans['h'])
        z=bank[str(c)];I=tuple(z['I']);W=tuple(z['W']);ic=z['ic'];wc=z['wc'];r=z['shift']
        hi,alpha,lam,_,_=profile(I,ic);hw,_,_,beta,d=profile(W,wc)
        need((alpha*beta,hi*d-hw*lam)==(99*2**r,c) and 0<=r<=4)
        U=ordinal(rows(word(k,'+')),tt(5-r));T=star(I,U,W,ic,wc)
        validate(T);C=[S for S in components(T)if len(S)>1]
        need(len(C)==1 and all(len(I)+j in C[0]for j in range(len(U))))
        core=induce(T,C[0]);need(len(components(core))==1)
        # Finite phase template itself is also counted directly on a singleton input.
        T0=star(I,tt(1),W,ic,wc);C0=[S for S in components(T0)if len(S)>1]
        need(len(C0)==1 and count(induce(T0,C0[0]))==99*2**(r+1)+c)
        out.append(dict(k=k,c=c,N=N,prime_cofactor=p,shift=r,strong_order=len(core),inverse=ans))
    # In this bounded range only the empty multiset can enter; compare all odd targets.
    negpos=0
    expected={99*2**t+c for t in range(5,8)for c in range(1,198,2)}
    for N in range(1,6501,2):need((inverse(N)is not None)==(N in expected));negpos+=1
    tests=0
    for k,sign in ((0,'-'),(0,'+'),(1,'+')):
        h,q=closed(k,sign)
        for c in range(1,198,2):
            N=99*32*q+c*h;z=inverse(N);need(z is not None);tests+=1
    # Exercise repeated and mixed letters, not only one-letter acceptance.
    for exponents in ((2,0),(1,1),(0,2),(2,1)):
        a,b=exponents;h=1141**a*1153**b;q=160240**(a+b);N=99*32*q+23*h
        z=inverse(N);need(z is not None);tests+=1
    print(json.dumps(dict(status='PASS',literal_profiles=small,formula_checks=formula,finite_target_tests=negpos,
        positive_entry_tests=tests,examples=out),sort_keys=True,indent=2))
if __name__=='__main__':main_sources()
