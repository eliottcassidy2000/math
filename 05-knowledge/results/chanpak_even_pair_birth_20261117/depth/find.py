from solve import *
import sys,json
from math import isqrt
def rooted_canon(a):
    from itertools import permutations,groupby
    groups=[list(g) for _,g in groupby(sorted(range(len(a)),key=lambda i:a[i].bit_count()),key=lambda i:a[i].bit_count())]
    return min(tuple(sum((a[i]>>j&1)<<t for t,j in enumerate(o)) for i in o)
      for parts in product(*(list(permutations(g)) for g in groups)) for o in (sum(parts,()),))

@cache
def tournaments(n):
    if n==1:return {(0,)}
    return {rooted_canon(apex(a,q)) for a in tournaments(n-1) for q in range(1<<(n-1))}

@cache
def banks(max_size):
    bank={}
    for n in range(3,max_size+1):
        for a in sorted(tournaments(n)):
            Ain=htable(a);Aout=htable(converse(a));M=(1<<n)-1;small=[1]+[sum(r) for r in Ain[1:]]
            # Strong deletion protects the intended geometry; singleton/external is omitted.
            reach=a[0]|1
            while True:
                r=reach
                for j in range(n):
                    if r>>j&1:reach|=a[j]
                if reach==r:break
            if reach!=M:continue
            reach=1
            while True:
                r=reach
                for j in range(n):
                    if a[j]&r:reach|=1<<j
                if reach==r:break
            if reach!=M:continue
            P=sum(small[s]*small[M^s] for s in range(M+1))
            p=[[[sum(A[s][i]*small[M^s] for s in range(M+1) if s.bit_count()%2==j) for i in range(n)] for j in range(2)] for A in (Aout,Ain)]
            for q in range(1,M):
                l,r=[[sum(v for i,v in enumerate(t) if q>>i&1==j) for t in B] for B,j in zip(p,(1,0))]
                h=small[-1];k=H(apex(a,q))
                S=(h,P,sum(l),sum(r),l[0]-l[1],r[0]-r[1],k)
                bank.setdefault(S,(a,q))
    return bank

def search(N,maxcore,budget,report=False):
    bank=banks(maxcore);T=0;C=0;G=0;Rpr=0
    by_size={i:[(S,a,q) for S,(a,q) in bank.items() if len(a)==i] for i in range(3,maxcore+1)}
    for order in range(4,budget+1):
        for i in range(3,min(maxcore,order-1)+1):
            t=order-i-1
            for S,a,q in by_size[i]:
                h,P,L,R,ls,rs,k=S
                B=N-k+max(ls,0)+max(rs,0)
                if B<0:continue
                for m in range(t+1):
                    n=t-m
                    # A graph actually at this order has full-length words.
                    u0=(1<<(m-1)) if m else 0;w0=(1<<(n-1)) if n else 0
                    U=(1<<m)-1;W=(1<<n)-1
                    if R*u0+L*w0+P*u0*w0>B:continue
                    if N>k+(h+R)*U+(h+L)*W+P*U*W+max(-ls,0)+max(-rs,0):continue
                    z,stats=dfs(S,N,m,n,full=True)
                    T+=stats[0];C+=stats[1];G+=stats[2];Rpr+=stats[3]
                    if z is None:continue
                    X=literal(a,q,*z);assert len(X)==order
                    if order<=16:assert H(X)==N
                    answer={'target':N,'order':order,'D':a,'q':q,'u':z[0],'w':z[1],'profile':S,'rows':X,'nodes':T,'cut_eval':C,'cut_prune':G,'residue_prune':Rpr}
                    if report:print(json.dumps(answer),flush=True)
                    return answer
    if report:print(json.dumps({'target':N,'miss_budget':budget,'nodes':T,'cut_eval':C,'cut_prune':G,'residue_prune':Rpr}),flush=True)

if __name__=='__main__':
    maxcore,budget=map(int,sys.argv[1:3])
    for N in map(int,sys.argv[3:]):search(N,maxcore,budget,True)
