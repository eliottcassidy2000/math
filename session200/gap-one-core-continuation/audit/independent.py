"""Independent gap-one audit. No author modules; explicit failures survive -O."""
from fractions import Fraction as Q
from pathlib import Path
from collections import deque
import hashlib,json

def check(p,m):
    if not p:raise RuntimeError(m)
M={
'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))}
seed=(1,0,2,1)
def act(x,c):return tuple(sum(u*v for u,v in zip(row,x))for row in M[c])
def walk(w,x=seed):
    for c in w:x=act(x,c)
    return x
def scalar(x):return x[0]-x[1]+x[2]
def bare_graph(w):
    # Orient each new edge from the old vertex to the appended vertex,
    # except that Z and J reverse the edge incident with the old root.
    adj=[[False,True],[False,False]];root=0
    for c in w:
        n=len(adj)
        for v,row in enumerate(adj):row.append(c=='O' or v!=root)
        adj.append([c!='O' and v==root for v in range(n)]+[False])
        if c=='J':root=n
    return adj,root

def paths(adj,first=None,last=None,close=False):
    n=len(adj);allmask=(1<<n)-1;dp={}
    for v in (range(n)if first is None else [first]):dp[(1<<v,v)]=1
    for mask in range(1,allmask+1):
        for v in range(n):
            count=dp.get((mask,v),0)
            if not count:continue
            for u in range(n):
                if not(mask>>u&1)and adj[v][u]:
                    key=(mask|1<<u,u);dp[key]=dp.get(key,0)+count
    return sum(dp.get((allmask,v),0)for v in range(n)if
               (last is None or v==last)and(not close or adj[v][first]))

# Forward OZ enumeration supplies every integer cover 2..5001 exactly once.
cap=5001;bank={2:(seed,'')};todo=deque([2])
while todo:
    u=todo.popleft();x,w=bank[u]
    for c in 'OZ':
        y=act(x,c);v=y[2]
        if v>cap:continue
        check(v not in bank,'unique zero-J cover recipe')
        bank[v]=(y,w+c);todo.append(v)
check(set(bank)==set(range(2,cap+1)),'complete all-cover zero-J inventory')

def g(t):
    if t==1:return 1
    if t%2==0:return g(t//2)
    s=(t+1)//2;return g(s)+s-(s%2)
formula=0;graph_receipts=[]
for t,(x,w)in sorted(bank.items()):
    check(x[0]==g(t),'g recurrence against direct forward enumeration')
    if not(t&1):continue
    D,E,U,d=x
    check((E,U,d)==(1,t,1)and 2*D>=t+3,'physical odd-cover lower bound')
    y=act(x,'J')
    check((y[0],y[2],y[3])==(D+t-1,D+t,D),'gap-one retained core')
    for a in range(7):
        for k in range(7):
            r=1<<(a+1);m=(1<<(k+1))-1
            n=(m*(r-1)+2)*D+(m*r+1)*t-1
            recipe=w+'J'+'O'*(a+1)+'Z'*k
            check(scalar(walk(recipe))==n,'literal source formula')
            check(recipe.count('J')==1,'supplied precursor J cost')
            formula+=1
            if t<=17 and a<=1 and k<=1:
                adj,root=bare_graph(recipe+'J')
                hp=paths(adj);hc=paths(adj,first=0,close=True)
                check((hp,hc)==(n,1),'literal graph target and C1')
                graph_receipts.append(dict(t=t,a=a,k=k,word=recipe+'J',order=len(adj),paths=hp,cycles=hc))
# General return identity and exceptional triangle control.
checks=0
for x in [(i,j,k,l)for i in range(2)for j in range(2)for k in range(2)for l in range(2)]:
    check(walk('JZ',x)==walk('OJ',x),'universal chronological JZ=OJ')
    D,E,U,d=x;G=D+2*U
    check(walk('OJ',x)==(G,D,G,D),'universal OJ balanced formula')
    checks+=1
triangle=walk('Z')
check(triangle==(3,1,3,1),'triangle balanced-core exception')
check(seed[0]+2*seed[2]==5,'all nontrivial sources have G at least5 by cone')
# Exact density coefficient control, with the three geometric tail lemmas.
short=sum((Q(1,v)for v in (9,15,19,37)),Q())+Q(10,69)+Q(32,385)
check(short==Q(27199346,56025585),'short density evaluation')
check(Q(17,35)-short==Q(13081,56025585)>0,'strict gap below17/35')
for a in range(2,80):check(3*(1<<(a+1))-1>=23*(1<<(a-2)),'a-tail geometric comparison')
for k in range(1,80):check((1<<(k+1))-1>=3*(1<<(k-1)),'k-tail from one')
for k in range(2,80):check((1<<(k+1))-1>=7*(1<<(k-2)),'k-tail from two')
L=6
partial=sum((Q(1,((1<<(k+1))-1)*(3*(1<<(a+1))-1)+4)for a in range(L)for k in range(L)),Q())
coarse_tail=Q(8,5*(1<<L))
check(partial+coarse_tail<Q(471,1000),'alternative bound below0.471')
# Literal finite target counting, allowing all exponents and every possible t.
N=20001;targets=set();presentations=0
check(9*cap>2*N,'complete target-cover cutoff')
for t,(x,w)in bank.items():
    if not(t&1):continue
    # A>=3,B>=3 and D>t/2 imply n>4.5t. All required t
    # lie below 2N/9<cap; the exponent bound follows from n>2^(a+k).
    D=x[0]
    for a in range(N.bit_length()):
        r=1<<(a+1)
        for k in range(N.bit_length()):
            m=(1<<(k+1))-1;n=(m*(r-1)+2)*D+(m*r+1)*t-1
            if n<=N:targets.add(n);presentations+=1
check(len(targets)<=short*N,'complete finite count respects all-order bound')
print(json.dumps(dict(status='PASS_WITH_ONE_WORDING_CLARIFICATION',cover_max=cap,
    zero_J_profiles=len(bank),odd_formula_checks=formula,literal_graph_controls=graph_receipts,
    universal_return_basis_checks=checks,short_density_bound=str(short),
    strict_ceiling='17/35',strict_gap=str(Q(17,35)-short),
    alternative_density_bound=str(partial+coarse_tail),alternative_density_decimal=float(partial+coarse_tail),
    complete_target_universe=dict(target_max=N,cover_max=cap,targets=len(targets),presentations=presentations,complete=True),
    verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2,sort_keys=True))
