#!/usr/bin/env python3
"""Finite-exact fixed-ray obstruction for <=4 operations; standard library only."""
import json
from pathlib import Path
from fractions import Fraction
from math import isqrt
from collections import Counter
from migration import compiler, tournament, multiply, extend
ROOT=Path(__file__).resolve().parent
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
OUTPUT=(2,0,1,1)

def require(condition,message):
    if not condition:raise RuntimeError(message)

def matmul(A,B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)) for i in range(4))

def determinant(A,lam=0):
    a=[[A[i][j]-(lam if i==j else 0) for j in range(4)]for i in range(4)]
    denominator=1;sign=1
    for k in range(3):
        if not a[k][k]:
            for j in range(k+1,4):
                if a[j][k]:a[k],a[j]=a[j],a[k];sign=-sign;break
            else:return 0
        pivot=a[k][k]
        for i in range(k+1,4):
            for j in range(k+1,4):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                require(numerator%denominator==0,'nonexact Bareiss division')
                a[i][j]=numerator//denominator
        denominator=pivot
    return sign*a[-1][-1]

def characteristic(A):
    B=I;answer=[1]
    for k in range(1,5):
        B=matmul(A,B);numerator=-sum(B[i][i]for i in range(4))
        require(numerator%k==0,'nonintegral characteristic coefficient')
        ck=numerator//k;answer.append(ck)
        B=tuple(tuple(B[i][j]+(ck if i==j else 0)for j in range(4))for i in range(4))
    require(B==((0,0,0,0),)*4,'Cayley-Hamilton failed')
    return answer

def evaluate(p,x):
    total=0
    for c in p:total=total*x+c
    return total

def positive_integer_roots(p):
    # Monic rational-root theorem after removing zero roots: every integer
    # root divides the remaining nonzero constant coefficient.
    p=list(p)
    while len(p)>1 and not p[-1]:p.pop()
    if len(p)==1:return []
    constant=abs(p[-1]);divisors=set()
    for d in range(1,isqrt(constant)+1):
        if constant%d==0:divisors.add(d);divisors.add(constant//d)
    return sorted(d for d in divisors if d>=2 and evaluate(p,d)==0)

def nullspace(A,lam):
    rows=[[Fraction(A[i][j]-(lam if i==j else 0))for j in range(4)]for i in range(4)]
    pivot_columns=[];r=0
    for c in range(4):
        pivot=next((i for i in range(r,4)if rows[i][c]),None)
        if pivot is None:continue
        rows[r],rows[pivot]=rows[pivot],rows[r]
        q=rows[r][c];rows[r]=[x/q for x in rows[r]]
        for i in range(4):
            if i!=r:
                q=rows[i][c];rows[i]=[x-q*y for x,y in zip(rows[i],rows[r])]
        pivot_columns.append(c);r+=1
    answer=[]
    for c in range(4):
        if c in pivot_columns:continue
        v=[Fraction(0)]*4;v[c]=Fraction(1)
        for i,p in enumerate(pivot_columns):v[p]=-rows[i][c]
        answer.append(tuple(v))
    return answer

def positive_matrix(A):
    P=((0,0,1,1),(2,0,1,1),(2,2,1,0),(0,0,1,0))
    cols=[]
    for j in range(4):
        h,F,q,g=multiply(A,tuple(P[i][j]for i in range(4)))
        require((F-h)%2==0 and (q-F+h-g)%2==0,'nonintegral conjugate')
        cols.append(((F-h)//2,(q-F+h-g)//2,g,h-g))
    return tuple(tuple(cols[j][i]for j in range(4))for i in range(4))

def regenerate_bank():
    bank={((2,1,1,0),(0,1,0,0),(0,0,0,0),(0,0,1,1)):'S0',
          ((1,0,0,0),(1,2,1,0),(0,0,1,1),(0,0,0,0)):'S1'}
    for n in range(2,5):
        for code in range(1<<(n*(n-1)//2)):
            local=tournament(n,code)
            if sorted(sum(local[i]>>j&1 for j in range(1,n))for i in range(1,n))!=list(range(n-1)):continue
            M=positive_matrix(compiler(local,n-1))
            require(all(x>=0 for row in M for x in row),'negative matrix')
            bank.setdefault(M,[n,code,n-1])
    frozen={tuple(map(tuple,x['matrix'])):x['recipe']for x in json.loads((ROOT/'bank.json').read_text())}
    require(bank==frozen,'frozen/compiler bank mismatch')
    require(len(bank)==32,'bank size changed')
    return bank

def actual_counts(A,cut):
    # Independent endpoint DP, not transfer formulas.
    n=len(A);size=1<<n;end=[[0]*n for _ in range(size)];start=[[0]*n for _ in range(size)]
    H=[0]*size;H[0]=1
    for mask in range(1,size):
        for v in range(n):
            if not mask>>v&1:continue
            rest=mask^(1<<v)
            if not rest:end[mask][v]=start[mask][v]=1
            else:
                end[mask][v]=sum(end[rest][u]for u in range(n)if rest>>u&1 and A[u]>>v&1)
                start[mask][v]=sum(start[rest][u]for u in range(n)if rest>>u&1 and A[v]>>u&1)
        H[mask]=sum(end[mask])
    q=sum(H[(size-1)^Q]*sum(start[Q][v]for v in range(n)if cut>>v&1)for Q in range(1,size))
    g=sum(end[-1][v]for v in range(n)if cut>>v&1)
    full=list(A)+[cut]
    for i in range(n):
        if not cut>>i&1:full[i]|=1<<n
    m=n+1;D=[[0]*m for _ in range(1<<m)]
    for mask in range(1,1<<m):
        for v in range(m):
            if mask>>v&1:
                rest=mask^(1<<v)
                D[mask][v]=1 if not rest else sum(D[rest][u]for u in range(m)if rest>>u&1 and full[u]>>v&1)
    return H[-1],sum(D[-1]),q,g

def positive_survivor(bank):
    word=[[3,0,2],[4,29,3],[3,5,2]];M=I;A=(0,);cut=0;v=(0,0,0,1);stages=[]
    for n,code,z in word:
        local=tournament(n,code);B=positive_matrix(compiler(local,z));M=matmul(B,M)
        A,cut=extend(A,cut,local,z);v=multiply(B,v)
        actual=actual_counts(A,cut);C,L,g,t=v
        require(actual==(g+t,2*C+g+t,2*C+2*L+g,g),'literal survivor count mismatch')
        stages.append({'recipe':[n,code,z],'positive_state':v,'native_state':actual,'parent_order':len(A)})
    require(M==((64,44,46,18),(0,0,0,0),(0,0,0,0),(56,38,41,17)),'survivor matrix changed')
    require(multiply(M,(1,0,0,0))==(64,0,0,56) and multiply(M,(0,0,0,1))==(18,0,0,17),'survivor face changed')
    v=(0,0,0,1);values=[]
    for k in range(12):
        C,L,g,t=v;F=2*C+g+t
        require(L==g==0 and 8*C-9*t==-9,'lost invariant')
        require(79*F==52*80**k+27,'family formula mismatch')
        require(79*C==18*(80**k-1) and 79*t==16*80**k+63,'state formula mismatch')
        values.append(F);v=multiply(M,v)
    return {'word':word,'matrix':M,'full_order_first_iterate':len(A)+1,'parent_adjacency_masks':A,'cut':cut,'stages':stages,'first_twelve_outputs':values}

def main():
    bank=regenerate_bank();cur={I:[]};eigen_one={};layers=[]
    for depth in range(1,5):
        nxt={}
        for M,word in cur.items():
            for B,recipe in bank.items():nxt.setdefault(matmul(B,M),word+[recipe])
        selected={M:w for M,w in nxt.items()if determinant(M,1)==0}
        for M,w in selected.items():eigen_one.setdefault(M,w)
        layers.append({'depth':depth,'distinct_matrices':len(nxt),'eigenvalue_one':len(selected)})
        cur=nxt
    require(layers==[{'depth':1,'distinct_matrices':32,'eigenvalue_one':7},{'depth':2,'distinct_matrices':690,'eigenvalue_one':15},{'depth':3,'distinct_matrices':14396,'eigenvalue_one':49},{'depth':4,'distinct_matrices':297239,'eigenvalue_one':143}],'finite universe changed')
    candidates=[];by_key={};root_summary=Counter();integer_pairs=0;zero_output=0;sign_changing=0
    for M,word in eigen_one.items():
        polynomial=characteristic(M)
        require(evaluate(polynomial,1)==0,'det/charpoly mismatch')
        for B in positive_integer_roots(polynomial):
            integer_pairs+=1
            require(determinant(M,B)==0,'integer root/det mismatch')
            basis=nullspace(M,B)
            require(len(basis)==1,'unhandled multidimensional eigenspace')
            a=basis[0];out=sum(c*x for c,x in zip(OUTPUT,a))
            if not out:zero_output+=1;continue
            a=tuple(x/out for x in a)
            if min(a)<0:sign_changing+=1;continue
            require(multiply(M,a)==tuple(B*x for x in a),'eigenray identity mismatch')
            key=(B,a)
            require(key not in by_key or by_key[key]==M,'shared same-radix positive eigenray')
            by_key[key]=M;root_summary[B]+=1
            candidates.append({'radix':B,'normalized_ray':[str(x)for x in a],'word':word,'matrix':M,'characteristic_coefficients':polynomial})
    require(integer_pairs==88 and zero_output==4 and sign_changing==0,'integer eigenspace census changed')
    require(len(candidates)==84,'candidate count changed')
    require(sum(x['radix']==8 for x in candidates)==7,'radix8 output-positive count changed')
    require(all(all(isinstance(r,str)for r in x['word'])for x in candidates if x['radix']==8),'new radix8 mixed macro')
    # Positive and hostile controls for helpers and scope.
    require(determinant(I,1)==0 and characteristic(I)==[1,-4,6,-4,1],'identity control')
    require(positive_integer_roots([1,-82,161,-80,0])==[80],'root-factor control')
    require(len(nullspace(I,1))==4,'nullspace dimension control')
    result={'status':'PASS','scope':'Words of length at most four in this fixed authenticated 32-matrix grammar. No all-length fixed-ray or unrestricted spectrum claim.',
      'bank_size':len(bank),'layers':layers,'distinct_eigenvalue_one_union':len(eigen_one),
      'integer_radix_eigenspaces':integer_pairs,'all_integer_radix_eigenspaces_one_dimensional':True,
      'zero_output_integer_eigenrays':zero_output,'sign_changing_normalized_integer_eigenrays':sign_changing,
      'positive_output_integer_radix_candidates':len(candidates),'radix_histogram':dict(sorted(root_summary.items())),
      'same_radix_normalized_ray_collisions':0,'candidates':sorted(candidates,key=lambda x:(x['radix'],x['normalized_ray'])),
      'positive_survivor':positive_survivor(bank)}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
