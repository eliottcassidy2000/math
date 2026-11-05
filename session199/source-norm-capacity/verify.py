"""Exact, standalone checks for source minima, norm chamber, and support scope."""
from itertools import product
from pathlib import Path
from math import log, sqrt
import hashlib
import json

HERE=Path(__file__).resolve().parent
O=((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1))
Z=((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0))
J=((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))
G={'O':O,'Z':Z,'J':J}; UROW=(1,1,2,2); SEED=(0,0,0,1)
P={
'N':((6,4,10,9),(4,3,7,6),(10,7,15,14),(9,6,14,13)),
'M':((11,8,20,18),(8,6,14,13),(20,14,34,31),(18,13,31,29)),
'A':((14,10,21,20),(10,7,16,15),(21,16,33,31),(20,15,31,29)),
'B':((14,10,26,23),(11,8,20,18),(20,15,35,32),(19,14,34,31)),
'C':((23,18,44,40),(18,14,34,31),(44,34,80,73),(40,31,73,67)),
'J0':((12,12,24,24),(12,12,21,21),(15,12,30,27),(15,12,30,27)),
'J1':((12,12,21,21),(12,12,18,18),(18,15,32,30),(18,15,31,29)),
}
P['T']=tuple(zip(*P['B']));P['J2']=tuple(zip(*P['J1']));P['J3']=tuple(zip(*P['J0']))
EXPECTED=[
{UROW},
{(2,2,3,3),(2,1,4,3)},
{(3,2,6,5),(4,4,5,5)},
{(5,4,10,9),(6,5,9,9),(6,3,11,9)},
{(9,6,18,15),(10,9,15,15),(10,5,19,15)},
{(15,10,30,25),(18,15,27,27),(18,9,33,27),(20,19,25,25)},
{(25,20,50,45),(27,18,54,45),(30,25,45,45),(30,15,55,45)},
{(45,30,90,75),(50,45,75,75),(50,25,95,75)},
]
CHECKS=0

def check(condition,msg):
    global CHECKS
    CHECKS+=1
    if not condition:raise RuntimeError(msg)

def rowmul(v,A):return tuple(sum(v[i]*A[i][j]for i in range(4))for j in range(4))
def colmul(A,v):return tuple(sum(a*b for a,b in zip(row,v))for row in A)
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def bilinear(v,A,w):return dot(v,colmul(A,w))
def old(v):
    a,f,e,p=v
    return a+2*e+p,e,a+f+2*e+2*p,e+p

def successors(v):return [colmul(A,v)for A in G.values()]

def minimal(rows):
    out=set()
    for r in sorted(set(rows),key=lambda a:(sum(a),a)):
        if not any(all(x<=y for x,y in zip(s,r))for s in out):out.add(r)
    return out

def step(front):return minimal(rowmul(r,A)for r in front for A in G.values())

def source(word):
    v=SEED
    for c in word:v=colmul(G[c],v)
    return v

def normcoords(s):
    D,E,U,d=s
    check((U+E)%2==0,'integral norm chart')
    return D-E+U,(U+E)//2,U-d

def rawN(v,w):
    D,E,U,d=v;g,e,L,f=w
    return (D-E+U)*g-(D-E+d)*e+(D+3*U-d)*L-(U+E)*f

def norm(v):
    x,y,z=v
    return x*x+4*y*z

def append(T,c,label):
    n=len(T); out=list(T)+[0]
    for v in range(n):
        if label!='O' and v==c:out[n]|=1<<v
        else:out[v]|=1<<n
    return out,n if label=='J' else c

def graph(word):
    T,c=[2,0],0
    for label in word:T,c=append(T,c,label)
    return T,c

def inventories(T,c):
    n=len(T);full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    rp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    rp[1<<c][c]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                todo=T[v] & (full^mask)
                while todo:
                    bit=todo&-todo;todo-=bit;w=bit.bit_length()-1
                    dp[mask|bit][w]+=dp[mask][v]
                    rp[mask|bit][w]+=rp[mask][v]
    H=[sum(row)for row in dp];H[0]=1
    rootH=[sum(row)for row in rp]
    U=sum(H[full^mask]*rootH[mask]for mask in range(1,1<<n))
    return (H[full],dp[full][c],U,H[full^(1<<c)]),rootH[full]

def HC(T):
    n=len(T);full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)];fixed=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    fixed[1][0]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                todo=T[v]&(full^mask)
                while todo:
                    bit=todo&-todo;todo-=bit;w=bit.bit_length()-1
                    dp[mask|bit][w]+=dp[mask][v]
                    fixed[mask|bit][w]+=fixed[mask][v]
    return sum(dp[full]),sum(fixed[full][v]for v in range(n)if T[v]&1)

def diagonal_N_graph(T,c):
    left,a=append(T,c,'J');m=len(left);n=len(T);r=m+c
    right=[sum((not(T[i]>>j&1))<<j for j in range(n)if j!=i)for i in range(n)]
    out=list(left)+[row<<m for row in right]
    for i in range(m):
        for j in range(n):
            if j==c and i in(a,c):out[m+j]|=1<<i
            else:out[i]|=1<<(m+j)
    return out

def main():
    front={UROW};fronts=[front]
    for k in range(1,61):
        front=step(front);fronts.append(front)
        if k<len(EXPECTED):check(front==EXPECTED[k],'complete printed Pareto state')
        if k>=7:check(front=={tuple(5*x for x in r)for r in fronts[k-3]},'exact scaled state recurrence')
    # Separate exhaustive universe: all 3^k rows, with no pruning of predecessor rows.
    exhaustive_counts=[]
    for k in range(11):
        rows=[]
        for word in product(G,repeat=k):
            r=UROW
            for c in word:r=rowmul(r,G[c])
            rows.append(r)
        check(minimal(rows)==fronts[k],'full unpruned word-row comparison')
        exhaustive_counts.append(len(rows))
    minimum=[]
    for k in range(61):
        m=k+2;mi=min(r[3]for r in fronts[k]);minimum.append(mi)
        if m==2:want=2
        else:
            j,r=divmod(m-3,3);want=(3,5,9)[r]*5**j
        check(mi==want,'sharp minimum by order')
        check(mi**3>=5**(m-1),'uniform cubic cover bound')
    for base in ('J','JZ','JZZ'):
        v=source(base)
        for repeat in range(20):
            h,e,u,d=old(v)
            check(h==u and e==d,'balanced extremal witness')
            check(u==minimum[len(base)+3*repeat],'extremal source recipe')
            w=v
            for c in 'JZZ':w=colmul(G[c],w)
            check(old(w)==(5*h,h,5*h,h),'all displayed clock witnesses')
            v=w
    for name,A in P.items():
        check(all(A[i][j]>=3*UROW[i]*UROW[j]for i in range(4)for j in range(4)),'all ten product bounds: '+name)
    # Raw four-coordinate basis confirms bilinear identity without a sampled-profile premise.
    for i,j in product(range(4),repeat=2):
        v=tuple(int(a==i)for a in range(4));w=tuple(int(a==j)for a in range(4))
        D,E,U,d=v;g,e,L,f=w
        X=D-E+U;Y=U+E;Z0=U-d;XX=g-e+L;YY=L+e;ZZ=L-f
        check(rawN(v,w)==X*XX+Y*ZZ+Z0*YY,'norm polarization exact coefficients')
    bank={SEED};queue=[SEED]
    for v in queue:
        D,E,U,d=old(v);x,y,z=normcoords((D,E,U,d))
        check(2*y<x<4*y and x%2==1 and y>=1 and z>=1,'physical chamber')
        check(rawN((D,E,U,d),(D,E,U,d))==norm((x,y,z)),'diagonal norm')
        out=(x-2*y,x-y+z,y)
        check(min(out)>0 and norm(out)==norm((x,y,z)) and out[0]<2*out[1],'outer involution leaves physical cone')
        X,Y,Z0=out
        check((X+2*Z0,Z0,Y-X-Z0)==(x,y,z),'inverse branch')
        for w in successors(v):
            check(old(w)[2]>U,'strict physical U growth')
            if old(w)[2]<=512 and w not in bank:bank.add(w);queue.append(w)
    # Literal path/cover counts independently authenticate all words through length four.
    literal=0
    for k in range(5):
        for word in product(G,repeat=k):
            T,c=graph(word);v,s=inventories(T,c)
            check(v==old(source(word))and s==1,'literal source inventories');literal+=1
    controls=[]
    for word in ('','J','JZ','JZZ','O','OJ','JJ','ZZ','JOZ'):
        T,c=graph(word);v,s=inventories(T,c);h,cy=HC(diagonal_N_graph(T,c))
        check(h==rawN(v,v)and cy==2,'literal diagonal C2')
        controls.append({'source_word':word,'source_order':len(T),'profile':v,'N':h,'C':cy})
    # Independent exact integer comparison for the strict entropy exponent.
    check(4*4*5<(25-9)**2,'9+4sqrt(5)<25')
    receipt={
        'status':'PASS','checks':CHECKS,'exhaustive_word_counts_by_depth':exhaustive_counts,
        'pareto_states_0_to_7':[sorted(s)for s in EXPECTED],
        'minimum_U_orders_2_to_22':minimum[:21],
        'periodic_frontier_depth_checked':60,'complete_source_bank_U_cap':512,
        'source_profiles':len(bank),'literal_source_words':literal,'literal_C2_controls':controls,
        'rho_display':3*log((3+sqrt(5))/2)/(2*log(5)),
        'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in [HERE/'PROOF.md',HERE/'verify.py']},
        'dependency_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((HERE/'dependencies').glob('*'))},
    }
    (HERE/'VALIDATION.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps(receipt,sort_keys=True,indent=2))

if __name__=='__main__':main()
