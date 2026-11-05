#!/usr/bin/env python3
from verify_unit_offset_rigidity import *
from fractions import Fraction
import random

def identity(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def matmul(a,b):return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def matvec(a,x):return [sum(v*w for v,w in zip(row,x)) for row in a]

def smith(matrix):
    a=[row[:] for row in matrix];m=len(a);n=len(a[0]);left=identity(m);right=identity(n)
    def rowadd(i,j,c):
        a[i]=[x+c*y for x,y in zip(a[i],a[j])];left[i]=[x+c*y for x,y in zip(left[i],left[j])]
    def coladd(i,j,c):
        for row in a:row[i]+=c*row[j]
        for row in right:row[i]+=c*row[j]
    def rowswap(i,j):a[i],a[j]=a[j],a[i];left[i],left[j]=left[j],left[i]
    def colswap(i,j):
        for row in a:row[i],row[j]=row[j],row[i]
        for row in right:row[i],row[j]=row[j],row[i]
    rank=0
    for k in range(min(m,n)):
        found=next(((i,j) for i in range(k,m) for j in range(k,n) if a[i][j]),None)
        if found is None:break
        rowswap(k,found[0]);colswap(k,found[1])
        while True:
            i=next((i for i in range(k+1,m) if a[i][k]),None)
            if i is not None:
                rowadd(i,k,-(a[i][k]//a[k][k]))
                if a[i][k]:rowswap(i,k)
                continue
            j=next((j for j in range(k+1,n) if a[k][j]),None)
            if j is not None:
                coladd(j,k,-(a[k][j]//a[k][k]))
                if a[k][j]:colswap(j,k)
                continue
            bad=next(((i,j) for i in range(k+1,m) for j in range(k+1,n) if a[i][j]%a[k][k]),None)
            if bad is None:break
            rowadd(k,bad[0],1)
        if a[k][k]<0:
            a[k]=[-x for x in a[k]];left[k]=[-x for x in left[k]]
        rank+=1
    return a,left,right,rank

def certify_smith(M,D,E,F,rank):
    require(matmul(matmul(E,M),F)==D,'Smith matrix identity')
    require(abs(det(E))==abs(det(F))==1,'Smith unimodularity')
    require(all(D[i][j]==0 for i in range(len(D)) for j in range(len(D[0])) if i!=j),'Smith diagonal')
    require(all(D[i][i]>0 for i in range(rank)),'positive Smith factors')
    require(all(D[i][i]%D[i-1][i-1]==0 for i in range(1,rank)),'Smith divisibility chain')
    require(all(not any(row) for row in D[rank:]),'zero Smith rows')

def solve_integer(M,t):
    D,E,F,rank=smith(M);certify_smith(M,D,E,F,rank);b=matvec(E,t)
    for i in range(len(M)):
        modulus=D[i][i] if i<rank else 0
        if (modulus and b[i]%modulus) or (not modulus and b[i]):
            row=E[i]
            vals=[sum(row[k]*M[k][j] for k in range(len(M))) for j in range(len(M[0]))]
            require(all(v%modulus==0 for v in vals) if modulus else not any(vals),'obstruction annihilates lattice')
            require(sum(x*y for x,y in zip(row,t))==b[i],'obstruction target value')
            return {'solvable':False,'matrix':M,'target':t,'D':D,'E':E,'F':F,'functional':row,'modulus':modulus,'target_value':b[i]}
    y=[0]*len(M[0])
    for i in range(rank):y[i]=b[i]//D[i][i]
    z=matvec(F,y);require(matvec(M,z)==t,'integer solution certificate')
    return {'solvable':True,'matrix':M,'target':t,'D':D,'E':E,'F':F,'z':z}

def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def valuations(q):
    require(q>0,'positive rational only');out=factors(q.numerator)
    for p,c in factors(q.denominator).items():out[p]=out.get(p,0)-c
    return {p:c for p,c in out.items() if c}

def solve_ratios(bank,kappa):
    require(bank and all(bank),'nonempty physical bank')
    ratios=[Fraction(profile(a)[1],profile(a)[0]) for a in bank]
    vs=list(map(valuations,ratios));target=valuations(kappa);primes=sorted(set(target).union(*(set(v) for v in vs)))
    M=[[v.get(p,0) for v in vs] for p in primes];t=[target.get(p,0) for p in primes]
    result=solve_integer(M,t);result.update(primes=primes,ratios=list(map(str,ratios)),kappa=str(kappa))
    return result

def construct_modules(bank,z):
    u=[max(x,0) for x in z];v=[max(-x,0) for x in z]
    if not sum(u) or not sum(v):u[0]+=1;v[0]+=1
    U=();V=()
    for a,cu,cv in zip(bank,u,v):
        for _ in range(cu):U=join(U,a)
        for _ in range(cv):V=join(V,a)
    require(U and V,'nonempty constructed modules');validate(U);validate(V)
    return U,V,u,v

def leaf_case(P,R):
    hp,bp,_=profile(P);hr,br,_=profile(R);dp,dr=bp-hp,br-hr;ep,er=4*hp-bp,4*hr-br
    require(dp>0 and dr>0,'positive leaf denominators')
    if ep==er==0:return 'all',None
    if ep*er>=0:return 'none',None
    return 'ratio',Fraction(-dr*ep,dp*er)

def run():
    K=(0,);W=(8,1,3,6);bank=[K,W];P=K;R=transitive(4)
    kind,kappa=leaf_case(P,R);require((kind,kappa)==('ratio',Fraction(5,2)),'physical positive target')
    yes=solve_ratios(bank,kappa);require(yes['solvable'] and yes['z']==[4,-1],'constructive lattice solution')
    U,V,u,v=construct_modules(bank,yes['z']);require(coeff(P,U,R,V)==(15360,0,-65),'physical coefficient recipe')
    graphs=[]
    for m in (1,2):
        a=graph(transitive(m),P,U,R,V);h=h_subsets(a)[-1];require(h==15360*3**m-65,'physical inverse witness')
        graphs.append({'m':m,'H':h,'rows':a})
    kind,kappa2=leaf_case(K,transitive(3));require(kappa2==Fraction(7,2),'physical missing-prime target')
    no_support=solve_ratios(bank,kappa2);require(not no_support['solvable'],'missing prime obstruction')
    no_index=solve_ratios([transitive(3),W],kappa);require(not no_index['solvable'],'Smith index obstruction')
    require((-1+5)%3!=0 and 3%3==0 and (5-5)%3==0,'explicit congruence functional')
    zero=solve_ratios(bank,Fraction(1));require(zero['solvable'],'zero valuation target')
    zU,zV,zu,zv=construct_modules(bank,zero['z']);require(zU and zV and Fraction(profile(zU)[1],profile(zU)[0])==Fraction(profile(zV)[1],profile(zV)[0]),'common padding')
    require(leaf_case(transitive(2),transitive(2))[0]=='all','both zero leaf case')
    require(leaf_case(K,K)[0]=='none','same-sign leaf case')
    require(leaf_case(transitive(2),K)[0]=='none','one-zero leaf case')
    rng=random.Random(197601);matrices=[[[2,0],[0,3]],[[6,0,0],[0,10,0],[0,0,15]],[[0,2],[4,6]],[[0,0],[0,0]],[[1,2],[2,4]]]
    for _ in range(200):
        m=rng.randrange(1,5);n=rng.randrange(1,5);matrices.append([[rng.randrange(-5,6) for _ in range(n)] for _ in range(m)])
    generated=arbitrary=0
    for M in matrices:
        z=[rng.randrange(-3,4) for _ in range(len(M[0]))];t=matvec(M,z)
        ans=solve_integer(M,t);require(ans['solvable'],'known integer target must solve');generated+=1
        ans=solve_integer(M,[rng.randrange(-7,8) for _ in range(len(M))]);arbitrary+=1
    print(json.dumps({'positive_case':yes,'module_copy_counts':[u,v],'physical_witnesses':graphs,'missing_prime_case':no_support,'index_obstruction_case':no_index,'common_padding_counts':[zu,zv],'known_integer_targets':generated,'arbitrary_certified_targets':arbitrary,'status':'PASS'},sort_keys=True,indent=2))
if __name__=='__main__':run()
