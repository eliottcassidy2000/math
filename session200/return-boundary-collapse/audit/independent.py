"""Independent matrix and exhaustive-permutation audit; no author imports."""
from fractions import Fraction
from itertools import permutations, product
from collections import Counter
import json

T=(1,0,2,1); B=(3,1,3,1); Z0=(1,0,2,0)
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
GEN={
 'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
 'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
 'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))}


def need(condition,message):
    if not condition:raise RuntimeError(message)


def mul(a,b):
    return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)) for i in range(4))


def apply(a,x):
    return tuple(sum(row[j]*x[j] for j in range(4)) for row in a)


def matrix(word):
    a=I
    for c in word:a=mul(GEN[c],a)
    return a


def scale(k,v):return tuple(k*x for x in v)
def plane(x):return x[0]==x[2] and x[1]==x[3]
def ray(x):return x==scale(x[1],B)


def rank(a):
    a=[[Fraction(x) for x in row] for row in a]
    r=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        q=a[r][col];a[r]=[x/q for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                q=a[i][col];a[i]=[x-q*y for x,y in zip(a[i],a[r])]
        r+=1
    return r


# Chronological JZ=OJ is an identity on the entire raw vector space.
need(mul(GEN['Z'],GEN['J'])==mul(GEN['J'],GEN['O']),'JZ=OJ')
returns=[];count=0;rank_counts=Counter();basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
plane_basis=[(1,0,1,0),(0,1,0,1)]
stack=[('',I)]
while stack:
    word,a=stack.pop()
    if word:
        count+=1;out=apply(a,B)
        if ray(out):
            lam=out[1];need(lam>=3 and lam%2==1,'return multiplier')
            # Reproduce the proof's global suffix normalization independently.
            zeros=len(word)-len(word.rstrip('Z'))
            if zeros:
                p=word[:-zeros]
                need(p and p[-1]=='J','last non-Z letter must be J')
                normalized=p[:-1]+'O'*zeros+'J'
            else:normalized=word
            need(normalized.endswith('OJ'),'return has normalized OJ suffix')
            need(matrix(normalized)==a,'global matrix normalization')
            prefix=matrix(normalized[:-2])
            need(all(plane(apply(prefix,x)) for x in plane_basis),'prefix preserves plane')
            need(all(plane(apply(a,x)) for x in basis),'whole-space image lies in plane')
            need(all(ray(apply(a,x)) for x in plane_basis),'plane image lies on triangle ray')
            r=rank(a);need(r in (1,2),'return matrix rank')
            rank_counts[r]+=1;returns.append((word,a,lam))
    if len(word)<10:
        for c in 'OZJ':stack.append((word+c,mul(GEN[c],a)))
need(count==88572 and len(returns)==5476,'complete declared word universe')
short=[item for item in returns if len(item[0])<=5]
mu_values=set()
for _,r1,l1 in short:
    for _,r2,l2 in short:
        composition=mul(r2,r1)
        need(rank(composition)==1,'two-return full matrix rank one')
        q=apply(composition,Z0);need(ray(q),'two-return auxiliary collapse')
        mu=q[1];need(l1*l2<2*mu<2*l1*l2,'physical joint ratio bounds')
        mu_values.add(mu)
        for k in range(1,4):
            bracket=l1*l2+2*((1<<k)-1)*mu
            need(bracket>=9,'nontrivial target factor')
        for _,r3,l3 in short:
            need(apply(r3,q)==scale(mu*l3,B),'later return preserves normalized slope')

# Exact two-letter menu: not an exponential family after two blocks.
A=matrix('OJOJ');Bmap=matrix('ZZOJ');states={Z0};state_counts={}
for depth in range(1,13):
    states={apply(a,q) for a in (A,Bmap) for q in states}
    if depth>=2:need(states=={scale(c*9**(depth-2),B) for c in (45,69,71)},'AB exact image')
    state_counts[depth]=len(states)

# Positive control uses a separate dense adjacency matrix and exhaustive permutations.
def literal_graph(word):
    adjacent=[[False,True],[False,False]];root=0
    for c in word:
        n=len(adjacent)
        for row in adjacent:row.append(False)
        adjacent.append([False]*(n+1))
        for i in range(n):
            incoming=(c!='O' and i==root)
            adjacent[n][i]=incoming;adjacent[i][n]=not incoming
        if c=='J':root=n
    return adjacent


def path(adj,p):return all(adj[a][b] for a,b in zip(p,p[1:]))
C=matrix('ZOJ');need(apply(C,B)==scale(5,B),'single return physical boundary')
need(apply(C,Z0)==(11,3,11,3) and not ray(apply(C,Z0)),'one return need not collapse')
need(rank(C)==2,'one-return threshold sharp in rank')
q=apply(C,Z0);h=q[0]-q[1]+q[2];m=25
need((m,h,m+6*h)==(25,19,139),'single-return target arithmetic')
adj=literal_graph('JZOOZZJ');n=len(adj)
for i in range(n):
    need(not adj[i][i],'no loops')
    for j in range(i):need(adj[i][j]!=adj[j][i],'tournament orientation')
paths=sum(path(adj,p) for p in permutations(range(n)))
cycles=sum(path(adj,(0,)+p) and adj[p[-1]][0] for p in permutations(range(1,n)))
need((n,paths,cycles)==(9,139,1),'independent literal Hamilton counts')
need(all(139%d for d in range(2,12)),'prime control')
result={
 'status':'PASS',
 'word_length_cap':10,'complete_words':count,'return_words':len(returns),
 'return_matrix_rank_counts':dict(rank_counts),
 'global_normalization_checks':len(returns),
 'short_return_menu':len(short),'two_return_rank_one_checks':len(short)**2,
 'three_return_scale_checks':len(short)**3,'AB_state_counts':state_counts,
 'positive_control':{'word':'JZOOZZJ','vertices':n,'H':paths,'C':cycles,'source_joint_readout':[m,h],'single_return_rank':rank(C),'method':'exhaustive vertex permutations'},
 'scope':'All-order proof audit plus exact finite matrix and graph controls; prime exclusion only for at least two return blocks; no universal target coverage.'
}
print(json.dumps(result,sort_keys=True,indent=2))
