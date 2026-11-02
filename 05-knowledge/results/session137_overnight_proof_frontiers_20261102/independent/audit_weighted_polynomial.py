"""Independent exact controls for averaged polynomial weighted stability."""
from fractions import Fraction as F
from itertools import combinations,product
from hashlib import sha256
import json,random
from audit_weighted_stability import matrix,uv,require

def powers_two(e):
    return F(2**e) if e>=0 else F(1,2**(-e))

def c(n):
    return (24-powers_two(3-n)*(n*n+2*n+3))/F(n*(n-1))

def beta(n):
    return (24-powers_two(4-n)*(n*n-n+3))/F(n*(n-1))

def all_h(a):
    n=len(a);N=1<<n;dp=[[F(0)]*n for _ in range(N)];h=[F(0)]*N;h[0]=F(1)
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,N):
        h[s]=sum(dp[s])
        for i,value in enumerate(dp[s]):
            if value:
                for j in range(n):
                    if not(s>>j&1):dp[s|1<<j][j]+=value*a[i][j]
    return h

def run():
    cases=weak_roots=strong_roots=exact_three=0
    digest=sha256();ratios={}
    require(c(2)==1 and beta(3)==1 and beta(4)==F(3,4),'coefficient bases')
    for n in range(3,201):
        t=powers_two(2-n)*(n-1);d=powers_two(3-n)
        require(c(n)==((n-2)*c(n-1)+2*t)/n,'weak recurrence closed form')
        require(beta(n)==((n-2)*c(n-1)+2*d)/n,'strong average closed form')
        require(beta(n)>=F(6,n*(n-1)),'simplified universal polynomial constant')
        if n>=4:
            require(n*(n-1)*beta(n)-(n-1)*(n-2)*beta(n-1)==powers_two(4-n)*(n*n-5*n+7),
                    'numerator successive difference')
    def audit(a):
        nonlocal cases,weak_roots,strong_roots,exact_three
        n=len(a);full=(1<<n)-1;h=all_h(a);u,v=uv(a);delta=h[full]-u;gap=delta-abs(v)
        V=sum(a[i][j]*a[j][i] for i,j in combinations(range(n),2))
        if n>=2:
            require(delta>=c(n)*V,'polynomial weak stability')
            for x in range(n):
                vs=[i for i in range(n) if i!=x];b=[[a[i][j] for j in vs] for i in vs]
                up,vp=uv(b);parent_delta=h[full^(1<<x)]-up
                incident=sum(a[x][i]*a[i][x] for i in range(n) if i!=x)
                require(delta>=parent_delta+powers_two(2-n)*(n-1)*incident,'root weak recurrence bound')
                weak_roots+=1
                if n>=3:
                    require(gap>=parent_delta+powers_two(3-n)*incident,'root strong recurrence bound')
                    strong_roots+=1
        if n>=3:
            require(gap>=beta(n)*V,'polynomial strong stability')
            require(gap>=F(6,n*(n-1))*V,'simplified polynomial stability')
            if V:ratios[n]=min(ratios.get(n,gap/V),gap/V)
        if n==3:
            require(gap==V,'exact order-three identity');exact_three+=1
        digest.update(str((n,delta,gap,V)).encode());cases+=1
    for n in range(1,5):
        for values in product((F(0),F(1,2),F(1)),repeat=n*(n-1)//2):audit(matrix(n,values))
    rng=random.Random(137245)
    for n in range(5,9):
        for _ in range(24):
            values=[]
            for _ in combinations(range(n),2):
                d=rng.choice((2,3,4,5,7));values.append(F(rng.randrange(d+1),d))
            audit(matrix(n,values))
    # Sharp n4: ordinal sum of two symmetric weighted pairs.
    a=matrix(4,(F(1,2),F(1),F(1),F(1),F(1),F(1,2)))
    h=all_h(a)[-1];u,v=uv(a)
    require((h,u,v)==(F(1),F(9,16),F(1,16)),'sharp order-four coordinates')
    require(h-u==c(4)*F(1,2) and h-u-abs(v)==beta(4)*F(1,2),'sharp n4 both constants')
    # Strong formula must not be extended down to order two.
    a=matrix(2,(F(1,2),));u,v=uv(a)
    require(all_h(a)[-1]-u-abs(v)==0 and beta(2)*F(1,4)>0,'n2 strong hostile')
    print(json.dumps({'scope':'all ternary weighted matrices n1..4;24 fixed-seed rational samples each n5..8;coefficient identities n3..200',
                      'matrices':cases,'weak_root_checks':weak_roots,'strong_root_checks':strong_roots,
                      'exact_order_three_checks':exact_three,
                      'sample_min_gap_over_variance':{n:str(x) for n,x in ratios.items()},
                      'sharp_n4':{'H':'1','u':'9/16','v':'1/16','V':'1/2','delta':'7/16','gap':'3/8'},
                      'semantic_sha256':digest.hexdigest(),
                      'verdict':'PASS polynomial coefficients, root recurrences, averaged bounds, n3 identity, n4 sharpness and n2 hostile'},sort_keys=True,indent=2))

if __name__=='__main__':run()
