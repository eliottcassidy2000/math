from itertools import combinations,product
from math import isqrt,gcd
from functools import lru_cache

def read_profile(T):
    n=len(T); F=(1<<n)-1
    d=[[0]*n for _ in range(F+1)]
    for i in range(n):d[1<<i][i]=1
    for m in range(1,F+1):
        for i,c in enumerate(d[m]):
            if not c: continue
            s=T[i]&~m
            while s:
                b=s&-s;s-=b
                d[m|b][b.bit_length()-1]+=c
    h=[1]+[sum(a) for a in d[1:]]
    return h[F],sum(h[m]*h[F^m] for m in range(F+1))

def one_apex(w):
    m=len(w);u=1;v=2**(m+1)
    for i in range(m):
        if w[i]!=1:continue
        for j in range(i+1,m):
            if w[j]!=0:continue
            d=j-i
            u+=1<<max(1,d-1)
            v+=(2**(m-d))*(3**(d-1)+1)
    return u,v

def atom_bank(A):
    ans=set();m=2
    while 2**(m+1)+2*(3**(m-2)+1)<=A:
        for w in product((0,1),repeat=m-2):
            u,v=one_apex((1,)+w+(0,))
            if v<=A:ans.add((u,v))
        m+=1
    return tuple(sorted(ans,key=lambda x:(x[1],x[0])))

def solve_small(N):
    # q>=2h; all singleton powers can fold into s.
    vmax=(2*N)//(2*3024-33)
    bank=atom_bank(vmax)
    reached={(1,1)}; pending=[(1,1)]
    for u,v in pending:
        for a,b in bank:
            if v*b>vmax:break
            x=(a*u,b*v)
            if x not in reached:reached.add(x);pending.append(x)
    ans=[]
    for u,v in reached:
        for I in (-33,-31,-29):
            q,r=divmod(N-I*u,756*v)
            if r==0 and q>=4 and not (q&(q-1)):
                ans.append((u,v,I,q.bit_length()-1))
    return ans,len(reached)

def prime(x):return x>1 and all(x%i for i in range(2,isqrt(x)+1))
def two_triangles(T):
    c=[set(q) for q in combinations(range(len(T)),3)
       if (T[q[0]]>>q[1]&1)==(T[q[1]]>>q[2]&1)==(T[q[2]]>>q[0]&1)]
    return next(((a,b) for a,b in combinations(c,2) if not a&b),None)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1:
        cand={}
        for line in open(sys.argv[1]):
            if not line.strip() or line.startswith('h\t'):continue
            h,n,rows=line.strip().split('\t')
            if int(n)>7:break
            T=tuple(map(int,rows.split(',')));z=two_triangles(T)
            if not z:continue
            u,v=read_profile(T)
            for s in range(2,6):
                for I in (-33,-31,-29):
                    N=756*(1<<s)*v+I*u
                    if prime(N) or (N%7==0 and prime(N//7)):cand[N]=(T,u,v,s,I,z)
        for N,t in sorted(cand.items()):
            x,k=solve_small(N)
            if not x:
                print(N,t,k)
                break

def controls():
    from fractions import Fraction as F
    ceil_coeff=[107,115,115,108,97,84,70,56]
    for m,c in enumerate(ceil_coeff,2):
        L=2**(m+1)+2*(3**(m-2)+1)
        assert (1000*2**(m-2))**10<c**10*L**9
    assert 3**9*25**10>67**10 and 2**9*537**10>1000**10
    rho=F(50,67)
    assert F(sum(ceil_coeff),1000)+F(537,1000)*rho**8/(1-rho)<1
    T=(2,4,25,19,35,15)
    assert two_triangles(T)
    assert read_profile(T)==(23,272)
    N=821861
    assert N==3024*272-29*23 and prime(N)
    assert solve_small(N)==([],9)
    # Check the derived one-apex profile formula against direct subset DP.
    for m in range(2,8):
        for w in product((0,1),repeat=m):
            rows=tuple(((1<<m)-(1<<(j+1)))|((1-w[j])<<m) for j in range(m))
            rows+=(sum(x<<j for j,x in enumerate(w)),)
            assert read_profile(rows)==one_apex(w)
    print('Exact rich-block target:',N,'achieved by',T,'with profile',(23,272))
    print('One-apex component policy rejected it after exhausting 9 profiles.')
if __name__=='__main__':controls()

if __name__=='__main__':
    assert {one_apex((1,)+w+(0,))[1]
            for m in range(2,7) for w in product((0,1),repeat=m-2)
            if one_apex((1,)+w+(0,))[0]==23}=={296}
    print('Every strong one-apex atom with H=23 has P(2)=296.')
