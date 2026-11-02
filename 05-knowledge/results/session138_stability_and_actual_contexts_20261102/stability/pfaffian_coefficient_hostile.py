from fractions import Fraction as F
import json

def pf(K,vs):
    if not vs:return F(1)
    a=vs[0]
    return sum((-1)**(j+1)*K[a][vs[j]]*pf(K,vs[1:j]+vs[j+1:]) for j in range(1,len(vs)))
def det(a):
    a=[list(map(F,row)) for row in a];d=F(1);n=len(a)
    for i in range(n):
        j=next((j for j in range(i,n) if a[j][i]),None)
        if j is None:return F(0)
        if j!=i:a[i],a[j]=a[j],a[i];d=-d
        x=a[i][i];d*=x
        for j in range(i+1,n):
            t=a[j][i]/x
            for k in range(i,n):a[j][k]-=t*a[i][k]
    return d

def matrix(x):
    K=[[F(0)]*6 for _ in range(6)]
    vals={(0,1):x,(0,2):-1,(0,3):-1,(0,4):-1,(0,5):-1,
          (1,2):1,(1,3):-1,(1,4):1,(1,5):1,
          (2,3):1,(2,4):1,(2,5):1,(3,4):1,(3,5):-1,(4,5):1}
    for (i,j),v in vals.items():K[i][j]=F(v);K[j][i]=-F(v)
    return K
x=F(2,3);K=matrix(x);A=[[(F(i!=j)+K[i][j])/2 for j in range(6)] for i in range(6)]
u=det([[A[i][j]+(i==j) for j in range(6)]for i in range(6)]);v=det(A)
require=lambda c,m: None if c else (_ for _ in ()).throw(RuntimeError(m))
require(pf(K,[2,3,4,5])==3,'coefficient3')
require(pf(matrix(0),list(range(6)))==-2,'interceptminus2')
require(pf(matrix(1),list(range(6)))==1,'slope3')
require(pf(K,list(range(6)))==0,'termwisehostile')
require(u==F(32,9) and v==F(-5,3),'actualmargins')
print(json.dumps({'scope':'FINITE-EXACT single weighted6-tournament, one uncertain edge',
 'skew_upper':[[str(K[i][j]) for j in range(i+1,6)]for i in range(6)],
 'full_pfaffian':'3x-2','x':str(x),'full_pfaffian_at_x':'0','complement_pfaffian':'3',
 'u':str(u),'v':str(v),'margin':str(u-abs(v)),
 'proposed_total_margin_bound':str(1-F(5,36)),
 'verdict':'REFUTES termwise Pfaffian square>=x^2; does NOT refute total-margin refinement'},sort_keys=True,indent=2))
