import itertools,json,math
from collections import defaultdict

def step(s,c,p,r):
    k=len(p)
    if r and s[r]==0 and not(c>>p[r-1]&1):return None
    t=list(s);t[r]=(2 if c>>p[r]&1 else 1) if r<k else 1
    return tuple(t)
def accept(s,p,A):
    k=len(p)
    return all(s[r]!=2 and (r==0 or s[r]!=0 or A[p[r-1]]>>p[r]&1) for r in range(k))
def transfer(A,w,block=False):
    total=0;k=len(A)
    for p in itertools.permutations(range(k)):
        d={(0,)*(k+1):1}
        for item in w:
            c,m=item if block else (item,1)
            e=defaultdict(int)
            for s,v in d.items():
                for mask in range(1<<(k+1)):
                    j=mask.bit_count();q=sum((-1)**(j-i)*math.comb(j,i)*i**m for i in range(j+1))
                    if not q:continue
                    t=s
                    for r in range(k+1):
                        if mask>>r&1:
                            t=step(t,c,p,r)
                            if t is None:break
                    if t is not None:e[t]+=v*q
            d=e
        total+=sum(v for s,v in d.items() if accept(s,p,A))
    return total

def rows(A,w):
    k=len(A);n=k+len(w);R=[0]*n
    for a in range(k):R[a]=A[a]
    for i,c in enumerate(w):
        x=k+i
        for a in range(k):
            if c>>a&1:R[a]|=1<<x
            else:R[x]|=1<<a
        for y in range(x+1,n):R[x]|=1<<y
    return R

def dp(R):
    n=len(R)
    if not n:return 1
    d={(1<<a,a):1 for a in range(n)}
    for size in range(1,n):
        e=defaultdict(int)
        for (mask,a),v in d.items():
            choices=R[a]&~mask
            while choices:
                b=choices&-choices;choices-=b;e[mask|b,b.bit_length()-1]+=v
        d=e
    return sum(d.values())
def check(x,y,msg):
    if x!=y:raise RuntimeError((msg,x,y))
def run():
    cases=0
    for k,maxm in [(0,5),(1,6),(2,5),(3,3)]:
        pairs=list(itertools.combinations(range(k),2))
        for bits in range(1<<len(pairs)):
            A=[0]*k
            for j,(a,b) in enumerate(pairs):
                if bits>>j&1:A[a]|=1<<b
                else:A[b]|=1<<a
            for m in range(maxm+1):
                for w in itertools.product(range(1<<k),repeat=m):
                    check(transfer(A,w),dp(rows(A,w)),('DP',A,w));cases+=1
    A=[2,4,1]
    vals=[transfer(A,[(7,3),(0,m)],True) for m in range(11)]
    checks=0
    for bits in range(8):
        A=[0]*3
        for j,(a,b) in enumerate(itertools.combinations(range(3),2)):
            if bits>>j&1:A[a]|=1<<b
            else:A[b]|=1<<a
        for c in range(8):
            for m in range(6):
                check(transfer(A,[(7,2),(c,m),(0,2)],True),transfer(A,[7]*2+[c]*m+[0]*2),'block');checks+=1
    return {'status':'PASS','exhaustive_subset_dp_cases':cases,'block_expansion_checks':checks,'cyclic_apex_rows':[2,4,1],'word':'7^3 0^m','values_m_0_through_10':vals}


# Two-block C3 specialization, independently evaluated as exponential polynomials.
def mul(P,Q):
    R=defaultdict(int)
    for (a,b),v in P.items():
        for (c,d),w in Q.items():R[a+c,b+d]+=v*w
    return dict(R)
def egf(A,r,s):
    k=len(A);Z={(1,1):1,(1,0):-1,(0,1):-1,(0,0):1}
    F={(0,0):1,(1,1):1,(1,0):-1};G={(0,0):1,(1,1):1,(0,1):-1}
    coeff=[0]*k
    for p in itertools.permutations(range(k)):
        d=sum(not(A[p[i]]>>p[i+1]&1) for i in range(k-1))
        for t in range(d,k):coeff[t]+=math.comb(k-1-d,t-d)
    P=mul(F,G);ans=0
    for t,v in enumerate(coeff):
        ans+=v*sum(c*a**r*b**s for (a,b),c in P.items());P=mul(P,Z)
    return ans

def extra():
    algebra=0
    for k in range(5):
        p=tuple(range(k))
        for s0 in itertools.product(range(3),repeat=k):
            for z in range(2):
                state=s0+(z,)
                for c in range(1<<k):
                    def E(s,r):return None if s is None else step(s,c,p,r)
                    for i in range(k+1):
                        check(E(E(state,i),i),E(state,i),'idempotent');algebra+=1
                        for j in range(i):check(E(E(state,i),j),E(E(state,j),i),'commute');algebra+=1
    ecases=0;leading_checks=0
    for k in range(1,5):
        pairs=list(itertools.combinations(range(k),2))
        for bits in range(1<<len(pairs)):
            A=[0]*k
            for j,(a,b) in enumerate(pairs):
                if bits>>j&1:A[a]|=1<<b
                else:A[b]|=1<<a
            # Apply (shift-1)...(shift-k) to positive-clock values at m=1.
            v=[egf(A,k,m) for m in range(1,k+2)]
            for j in range(1,k+1):v=[v[i+1]-j*v[i] for i in range(len(v)-1)]
            check(v[0],math.factorial(k)**3*(k+1),'highest clock coefficient');leading_checks+=1
            for r in range(4):
                for s in range(4):
                    check(egf(A,r,s),transfer(A,[((1<<k)-1,r),(0,s)],True),'EGF');ecases+=1
    A=[2,4,1]
    for m in range(1,21):check(egf(A,3,m),87-189*2**m+72*3**m+36*4**m,'cyclic formula')
    w=[7,1,2,4,0];R=rows(A,w);h=dp(R);check(h,transfer(A,w),'mixed mask')
    return {'highest_clock_coefficient_checks':leading_checks,'idempotence_commutation_checks':algebra,'EGF_vs_block_checks':ecases,'closed_formula_checks':20,'mixed_mask_witness':{'apex_rows':A,'word':w,'rows':R,'H':h}}

if __name__=='__main__':
    result=run();result.update(extra());print(json.dumps(result,indent=2,sort_keys=True))
