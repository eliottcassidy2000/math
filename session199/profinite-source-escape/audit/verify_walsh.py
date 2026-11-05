"""Independent exact Walsh-source audit. No external packages or repo imports."""
from itertools import product
import json

O=((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1))
Z=((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0))
J=((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))
GEN={'O':O,'Z':Z,'J':J}
b=(0,0,0,1)
PN=((6,4,10,9),(4,3,7,6),(10,7,15,14),(9,6,14,13))
PM=((11,8,20,18),(8,6,14,13),(20,14,34,31),(18,13,31,29))

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def add(*vectors): return tuple(map(sum,zip(*vectors)))
def mul(A,v): return tuple(sum(x*y for x,y in zip(row,v)) for row in A)
def dot(u,v): return sum(x*y for x,y in zip(u,v))
def mm(A,B): return tuple(tuple(sum(A[i][h]*B[h][j] for h in range(4)) for j in range(4)) for i in range(4))
def matrix(*terms): return tuple(tuple(sum(s*A[i][j] for s,A in terms) for j in range(4)) for i in range(4))
def old(v):
    a,f,e,p=v
    return (a+2*e+p,e,a+f+2*e+2*p,e+p)
def state(word):
    y=b
    for c in word:y=mul(GEN[c],y)
    return y

def build(word):
    # Fixed vertices 0->1, root0. Every newly added vertex loses to all
    # existing vertices except the current root for Z/J.
    rows=[2,0]; root=0
    for c in word:
        n=len(rows); rows.append(0)
        for i in range(n):
            if c!='O' and i==root: rows[n] |= 1<<i
            else: rows[i] |= 1<<n
        if c=='J':root=n
    return rows,root

def direct_profile(rows,root):
    n=len(rows); size=1<<n
    end=[[0]*n for _ in range(size)]
    startroot=[[0]*n for _ in range(size)]
    for i in range(n):end[1<<i][i]=1
    startroot[1<<root][root]=1
    for mask in range(1,size):
        for j in range(n):
            if not (mask>>j&1):continue
            prev=mask^(1<<j)
            if prev:
                end[mask][j]=sum(end[prev][i] for i in range(n) if prev>>i&1 and rows[i]>>j&1)
                startroot[mask][j]=sum(startroot[prev][i] for i in range(n) if prev>>i&1 and rows[i]>>j&1)
    H=[1]+[sum(row) for row in end[1:]]
    full=size-1
    U=sum(H[mask]*sum(startroot[full^mask]) for mask in range(size) if not(mask>>root&1))
    return (H[full],end[full][root],U,H[full^(1<<root)]),sum(startroot[full])

def formula_z(k):return (0,2*k,0,1)
def formula_j(k):return (0,3**k-1,0,1)
def formula_zj(k):
    return (-2*k,0,2*k,1) if k%2==0 else (2*k-4,2,-2*k,1)

KZ=matrix((1,O),(-1,Z),(1,J))
KJ=matrix((1,O),(1,Z),(-1,J))
KZJ=matrix((1,O),(-1,Z),(-1,J))
checks={}
for name,K,fn in [('Z',KZ,formula_z),('J',KJ,formula_j),('ZJ',KZJ,formula_zj)]:
    require(fn(0)==b,'base '+name)
    for k in range(200):require(mul(K,fn(k))==fn(k+1),'closed form '+name)
checks['formula_recurrences_checked_through']=200
for k in range(9):
    totals={name:(0,0,0,0) for name in ['Z','J','ZJ']}
    for word in product('OZJ',repeat=k):
        y=state(word)
        for name in totals:
            sign=(-1)**sum(word.count(c) for c in name)
            totals[name]=add(totals[name],tuple(sign*x for x in y))
    for name,fn in [('Z',formula_z),('J',formula_j),('ZJ',formula_zj)]:
        require(totals[name]==fn(k),'word aggregate')
checks['literal_word_layers_checked_through']=8
count=0
for k in range(6):
    for word in product('OZJ',repeat=k):
        prof,starts=direct_profile(*build(word))
        require(prof==old(state(word)),'independent direct graph profile')
        require(starts==1,'physical root start')
        count+=1
checks['direct_graphs_checked']=count
for r in range(20):
    for s in range(20):
        require(dot(formula_z(r),mul(PN,formula_z(s)))==12*(r+1)*(s+1)+1,'N double Walsh')
        require(dot(formula_z(r),mul(PM,formula_z(s)))==24*r*s+26*r+26*s+29,'M double Walsh')
checks['bilinear_grid_checked']=[20,20]
for m in range(10):
    require(formula_z(2**m-1)==state('O'*m),'physical dyadic compiler')
require(mm(J,O)==mm(Z,J),'word matrix overlap')
A,ra=build('OJ');B,rb=build('JZ');perm=(0,1,3,2)
require(perm[ra]==rb,'root isomorphism')
require(all(((A[i]>>j)&1)==((B[perm[i]]>>perm[j])&1) for i in range(4) for j in range(4)),'graph isomorphism')
checks['opposite_sign_rooted_isomorphism']={'word1':'OJ','word2':'JZ','map':perm,'profile':old(state('OJ'))}
checks['nonphysical_Z_aggregate_k2']=old(formula_z(2))
checks['balanced_signed_readout_k4']=dot((3,2,6,5),formula_z(4))
checks['status']='PASS: indexed-word identities; no support or coverage inference'
print(json.dumps(checks,indent=2,sort_keys=True))
