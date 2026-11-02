"""Exact rational weighted-tournament audit, including beneficial 2-cycles.

All n<=4 with unordered-edge probabilities in {0,1/2,1}; sixteen fixed-seed
rational matrices at each n=5,6,7. Full Fraction Gaussian determinants,
weighted subset path DP, and literal insertion paths are separate routes.
"""
from fractions import Fraction as F
from itertools import combinations,product
from hashlib import sha256
import random,json
from audit_feedback2 import determinant

def require(c,m):
    if not c:raise RuntimeError(m)

def run():
    count=roots=paths=two_cycles=empty_terms=subsets=0;digest=sha256()
    def audit(a):
        nonlocal count,roots,paths,two_cycles,empty_terms,subsets
        n=len(a);N=1<<n;full=N-1;h=[F(0)]*N;u=[F(0)]*N;v=[F(0)]*N
        h[0]=u[0]=v[0]=F(1);dp=[[F(0)]*n for _ in range(N)]
        for i in range(n):dp[1<<i][i]=1
        for s in range(1,N):
            h[s]=sum(dp[s])
            for j in range(n):
                if dp[s][j]:
                    for k in range(n):
                        if not(s>>k&1):dp[s|1<<k][k]+=dp[s][j]*a[j][k]
            vs=[i for i in range(n) if s>>i&1]
            b=[[a[i][j] for j in vs] for i in vs]
            v[s]=(-1)**len(vs)*determinant(b)
            u[s]=determinant([[b[i][j]+(i==j) for j in range(len(vs))] for i in range(len(vs))])
            require(u[s]+v[s]>0 and u[s]-v[s]>0,'weighted algebraic positivity')
            require(h[s]>=u[s]+abs(v[s]),'weighted inequality');subsets+=1
        for x in range(n):
            parent=full^(1<<x);eu=ou=sv=eh=gp=gm=F(0)
            def visit(last,used,weight):
                nonlocal eu,ou,sv,eh,gp,gm,paths,two_cycles,empty_terms
                if a[last][x]:
                    length=used.bit_count();rest=parent^used;w=weight*a[last][x]
                    sv+=w*v[rest];paths+=1;two_cycles+=length==1;empty_terms+=rest==0
                    if length%2==0:
                        eu+=w*u[rest];eh+=w*h[rest]
                        gp+=w*(2*h[rest]-u[rest]+v[rest]);gm+=w*(2*h[rest]-u[rest]-v[rest])
                    else:
                        ou+=w*u[rest];gp+=w*(u[rest]+v[rest]);gm+=w*(u[rest]-v[rest])
                for j in range(n):
                    if j!=x and not(used>>j&1) and a[last][j]:visit(j,used|1<<j,weight*a[last][j])
            for j in range(n):
                if a[x][j]:visit(j,1<<j,a[x][j])
            require(u[full]==u[parent]+eu-ou,'weighted u insertion')
            require(v[full]==-sv,'weighted v insertion')
            require(h[full]==h[parent]+2*eh,'weighted H insertion')
            require(h[full]-u[full]-v[full]==h[parent]-u[parent]+gp,'weighted plus gap')
            require(h[full]-u[full]+v[full]==h[parent]-u[parent]+gm,'weighted minus gap')
            roots+=1
        fractional=any(a[i][j] not in (0,1) for i,j in combinations(range(n),2))
        triangle=any(a[i][j]==a[j][k]==a[k][i]==1
                     for i in range(n) for j in range(n) for k in range(n))
        expected_P_equality=(n<=2 or (not fractional and (not triangle or n==3)))
        require((h[full]==u[full]+abs(v[full]))==expected_P_equality,'weighted full equality classification')
        require((h[full]==u[full]+v[full])==(not fractional and not triangle),'weighted disc equality classification')
        digest.update(json.dumps([n,str(h[full]),str(u[full]),str(v[full])]).encode());count+=1
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for weights in product((F(0),F(1,2),F(1)),repeat=len(pairs)):
            a=[[F(0)]*n for _ in range(n)]
            for (i,j),p in zip(pairs,weights):a[i][j]=p;a[j][i]=1-p
            audit(a)
    rng=random.Random(1371950)
    for n in (5,6,7):
        for _ in range(16):
            a=[[F(0)]*n for _ in range(n)]
            for i,j in combinations(range(n),2):
                denominator=rng.choice([2,3,4,5]);p=F(rng.randrange(denominator+1),denominator)
                a[i][j]=p;a[j][i]=1-p
            audit(a)
    print(json.dumps({'scope':'FINITE-EXACT rational weighted tournaments: all ternary edge weights n<=4;16 fixed-seed samples each n5..7',
                      'matrices':count,'root_insertions':roots,'positive_weight_old_paths':paths,
                      'two_cycle_terms':two_cycles,'empty_complement_terms':empty_terms,
                      'induced_determinants':subsets,'semantic_sha256':digest.hexdigest(),
                      'verdict':'PASS weighted insertion identities, H>=det(I+A)+|detA|, and complete weighted equality classification'},sort_keys=True,indent=2))
if __name__=='__main__':run()
