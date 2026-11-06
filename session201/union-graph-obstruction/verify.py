"""Portable exact controls for strict central increase and the union minimum."""
from pathlib import Path
from functools import lru_cache
from itertools import product
import hashlib
import json

P=((1,0,2,1),(0,0,1,0),(1,1,2,2),(0,0,1,1))
Q=((0,0,0,1),(0,1,0,-1),(0,0,1,-1),(1,-1,-1,1))
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
T=(1,0,2,1)

def need(c,m):
    if not c:raise RuntimeError(m)
def tr(A):return tuple(zip(*A))
def mm(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)) for i in range(4))
def mv(A,x):return tuple(sum(a*b for a,b in zip(row,x)) for row in A)
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def sub(A,B):return tuple(tuple(A[i][j]-B[i][j] for j in range(4)) for i in range(4))
def step(s,c):
    D,E,U,d=s
    if c=='O':return D,0,2*U,d
    if c=='Z':return D-E+U,d,2*U-d,d
    if c=='J':return D-E+U,D-E,D+U,D
    raise ValueError(c)
def positive(s):
    D,E,U,d=s
    return D-E-d,U-D+E-d,E,d-E
GEN={c:tr(tuple(step(x,c) for x in I)) for c in 'OZJ'}
def matrix(w):
    M=I
    for c in w:M=mm(GEN[c],M)
    return M
def state(w,s=T):
    for c in w:s=step(s,c)
    return s
def K(s,t):return dot(s,mv(Q,t))
@lru_cache(None)
def g(n):
    if n==1:return 1
    if n%2==0:return g(n//2)
    h=(n+1)//2;return g(h)+h-h%2

def tau(n):return (g(n),n%2,n,1)
def main():
    root=Path(__file__).resolve().parent;path=root/'INPUT-MENU.json'
    menu=json.loads(path.read_text());residuals={};cone={}
    for c,M in GEN.items():
        C=tr(tuple(positive(mv(M,mv(P,x))) for x in I))
        need(mm(P,C)==mm(M,P),'exact positive-coordinate transform')
        need(all(z>=0 for row in C for z in row),'cone preserved')
        # New d=e+p coefficientwise dominates old d.
        need(all(C[2][j]+C[3][j]>=int(j>=2) for j in range(4)),'root size does not decrease')
        R=mm(mm(tr(P),Q),mm(sub(M,I),P))
        need(all(z>=0 for row in R for z in row),'nonnegative increment matrix')
        need(all(R[i][j]>=2 for i in (2,3) for j in (2,3)),'strict lower-right block')
        residuals[c]=R;cone[c]=C
    expected_J=((1,0,2,1),(0,0,1,1),(2,1,3,3),(1,1,3,2))
    need(residuals['J']==expected_J,'inherited J residual checked literally')
    for entry in menu:
        need([list(row) for row in mm(Q,matrix(entry['word']))]==entry['matrix'],'local menu is the physical menu')
    need(len(menu)==12 and sum(e['word']=='' for e in menu)==1,'twelve forms include empty exactly once')
    physical=[]
    for length in range(6):
        for letters in product('OZJ',repeat=length):
            word=''.join(letters);s=state(word);y=positive(s)
            need(min(y)>=0 and y[2]+y[3]>=1 and mv(P,y)==s,'actual physical cone')
            physical.append(s)
    increment_checks=0
    for s in physical:
        for t in physical[::13]:
            for c in 'OZJ':
                need(K(s,step(t,c))-K(s,t)>=2*s[3]*t[3],'actual strict increase')
                increment_checks+=1
    minchecks=0
    for u in range(2,34):
        for v in range(2,34):
            vals=[K(tau(u),state(e['word'],tau(v))) for e in menu]
            empty=K(tau(u),tau(v))
            need(sum(n==empty for n in vals)==1 and min(vals)==empty,'empty is unique minimum of union fiber')
            for e,n in zip(menu,vals):need(n-empty>=2*len(e['word']),'word-length strict increase')
            if u==2:need(empty==g(v)+v-v%2,'nonsynchronized minimum readout')
            minchecks+=1
    return {'status':'PASS','positive_coordinate_generators':cone,'increment_matrices':residuals,
            'physical_states_checked':len(physical),'strict_increment_checks':increment_checks,
            'actual_union_minimum_checks':minchecks,
            'frozen_menu_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'scope':'all-order union witness graph obstruction; projected value support remains OPEN'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
