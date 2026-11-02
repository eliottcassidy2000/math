"""Independent Fraction elimination of the two-feedback-vertex determinant bound."""
from fractions import Fraction as F
from tournament_kernel import *
import json

def determinant(a):
    n=len(a);a=[list(map(F,row)) for row in a];answer=F(1)
    for c in range(n):
        p=next((r for r in range(c,n) if a[r][c]),None)
        if p is None:return F(0)
        if p!=c:a[p],a[c]=a[c],a[p];answer=-answer
        pivot=a[c][c];answer*=pivot
        a[c]=[x/pivot for x in a[c]]
        for r in range(c+1,n):
            factor=a[r][c]
            a[r]=[x-factor*y for x,y in zip(a[r],a[c])]
    return answer

def tournament(m,word):
    n=m+2;a=[[0]*n for _ in range(n)]
    for i,j in combinations(range(m),2):a[i][j]=1
    a[m][m+1]=1
    for i in range(m):
        for k in range(2):
            b=(word>>(2*i+k))&1;a[i][m+k]=b;a[m+k][i]=1-b
    return tuple(map(tuple,a))

def run():
    rows=[];cases=0
    for m in range(0,6):
        strong_cases=0;min_slack=None;maxu=0;maxv=0
        for w in range(4**m):
            a=tournament(m,w);n=m+2
            u=determinant([[a[i][j]+int(i==j) for j in range(n)] for i in range(n)])
            v=(-1)**n*determinant(a);P=u+abs(v)
            b=[[a[i][m+k] for k in range(2)] for i in range(m)]
            S=[[int(i==j)+a[m+i][m+j]-sum((1-b[l][i])*(b[l][j]-(b[l+1][j] if l+1<m else 0)) for l in range(m)) for j in range(2)] for i in range(2)]
            require(u==S[0][0]*S[1][1]-S[0][1]*S[1][0],'Schur determinant')
            for k in range(2):require(S[k][k]==1+sum(b[l][k]==0 and b[l+1][k]==1 for l in range(m-1)),'01 transition diagonal')
            require(all(abs(S[i][j])<=m+1 for i,j in [(0,1),(1,0)]),'offdiagonal bound')
            require(u<=(1+m//2)**2+(m+1)**2,'u bound')
            require(abs(v)<=4*(m+1),'v bound')
            h=hp(a);require(h>=P,'H>=P finite feedback2')
            if strong(a):
                strong_cases+=1;slack=h-P;min_slack=slack if min_slack is None else min(min_slack,slack)
            maxu=max(maxu,u);maxv=max(maxv,abs(v));cases+=1
        rows.append({'n':m+2,'presentations':4**m,'strong_presentations':strong_cases,'min_strong_H_minus_P':None if min_slack is None else str(min_slack),'max_u':str(maxu),'max_abs_v':str(maxv)})
    print(json.dumps({'scope':'FINITE-EXACT every oriented two-apex attachment word with fixed apex arc, transitive spine length0..5',
                      'presentations':cases,'summary':rows},sort_keys=True,indent=2))
if __name__=='__main__':run()
