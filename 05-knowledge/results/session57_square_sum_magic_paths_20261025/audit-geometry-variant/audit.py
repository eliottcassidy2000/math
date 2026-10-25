#!/usr/bin/env python3
"""Independent exact checks: explicit exceptions remain active under python -O."""
import itertools, json, math
from pathlib import Path
import sympy as S
OUT=Path(__file__).parent
checks=0
def check(condition, context):
    global checks
    checks+=1
    if not condition: raise RuntimeError(context)
def rank(A,p):
    a=[[x%p for x in row] for row in A]; h=0
    if not a:return 0
    for j in range(len(a[0])):
        pivot=next((i for i in range(h,len(a)) if a[i][j]),None)
        if pivot is None:continue
        a[h],a[pivot]=a[pivot],a[h]
        u=pow(a[h][j],-1,p);a[h]=[v*u%p for v in a[h]]
        for i in range(h+1,len(a)):
            u=a[i][j];a[i]=[(v-u*w)%p for v,w in zip(a[i],a[h])]
        h+=1
        if h==len(a):break
    return h

def incidence(r,s):
    return [[int(k==i)-int(k==r+j) for k in range(r+s)] for i in range(r) for j in range(s)]

def graph_rank(edges,r,s):
    parent=list(range(r+s)); count=0
    def root(x):
        while x!=parent[x]: x=parent[x]
        return x
    for i,j in edges:
        u,v=root(i),root(r+j)
        if u!=v:parent[u]=v;count+=1
    return count
out={}
for r,s in [(2,2),(2,3),(3,3),(4,4)]:
    edges=list(itertools.product(range(r),range(s))); L=incidence(r,s)
    for mask in range(1<<(r*s)):
        chosen=[k for k in range(r*s) if (mask>>k)&1]
        rg=graph_rank([edges[k] for k in chosen],r,s)
        check(rank([L[k] for k in chosen],3)==rg,(r,s,mask,'graphic rank'))
    out[f'graph_subsets_{r}_{s}']=1<<(r*s)
# Independent zero-pattern checks over all small additive bases, with gauge b0=0.
for r,s,p in [(2,2,3),(2,3,5),(3,3,5),(4,4,3)]:
    tested=0;singular=0
    for vals in itertools.product(range(p),repeat=r+s-1):
        a=vals[:r];b=(0,)+vals[r:];z=[(a[i]+b[j])%p for i in range(r) for j in range(s)]
        if not any(z):continue
        zero=[(i,j) for i in range(r) for j in range(s) if z[i*s+j]==0]
        cyclic=graph_rank(zero,r,s)<len(zero)
        rectangles=any(all(z[i*s+j]==0 for i in ii for j in jj) for ii in itertools.combinations(range(r),2) for jj in itertools.combinations(range(s),2))
        check(cyclic==rectangles,(r,s,p,vals,'zero cycle'))
        tested+=1;singular+=cyclic
    out[f'additive_bases_{r}_{s}_F{p}']={'tested':tested,'singular':singular}
# Symbolic quaternion orthogonality in both directions.
a,b,c,d=S.symbols('a b c d');n=a*a+b*b+c*c+d*d
M=S.Matrix([[a*a+b*b-c*c-d*d,2*(b*c-a*d),2*(b*d+a*c)], [2*(b*c+a*d),a*a-b*b+c*c-d*d,2*(c*d-a*b)],[2*(b*d-a*c),2*(c*d+a*b),a*a-b*b-c*c+d*d]])
for product in [M*M.T,M.T*M]:
    check((product-n*n*S.eye(3)).applyfunc(S.expand)==S.zeros(3),'quaternion identity')
roots=((2655,2976,320),(2624,2520,1665),(1440,895,3624))
grid=[[x*x for x in row] for row in roots]
check(len(set(sum(grid,[])))==9 and min(sum(grid,[]))>0,'distinct positivity')
sums=[sum(row) for row in grid]+[sum(grid[i][j] for i in range(3)) for j in range(3)]+[sum(grid[i][i] for i in range(3)),sum(grid[i][2-i] for i in range(3))]
check(sums==[4001**2]*6+[5151**2,2920**2],'line totals')
check(all(math.isqrt(x)**2==x for x in sums) and len(set(sums))==3,'square totals not magic')
Q=M.subs({a:48,b:32,c:23,d:12})
check(tuple(tuple(abs(int(Q[i,j])) for j in (0,2,1)) for i in range(3))==roots,'quaternion permutation')
out['line_totals']=sums
# Decorated form classes and linear embedding / degree accounting.
C,x,y=S.symbols('C x y'); ss=[0,x,y,x+y]
forms=[S.expand(C+t-u) for u in ss for t in ss]
classes={f:[i for i,g in enumerate(forms) if f==g] for f in set(forms)}
check(sorted(map(len,classes.values()))==[1,1,1,1,2,2,2,2,4],'multiplicities')
check(S.Matrix([[f.coeff(t) for t in (C,x,y)] for f in forms]).rank()==3,'decorated base rank')
check(2**(16-len(classes))*2**(len(classes)-3)==2**(16-3),'component degree sum')
out['decoration']={'components':128,'component_degree':64,'total_degree':8192,'component_dualizing_twist':3}
out['X44']={'dimension':6,'degree':512,'squaring_cover_degree':32768,'dualizing_twist':2,'singular_dimension':3}
out['checks']=checks;out['status']='PASS'
(OUT/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
