#!/usr/bin/env python3
"""Independent exact audit: no producer/canonical imports and no Git operations."""
import itertools, math, json
from collections import defaultdict, Counter
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
OUT=Path(__file__).resolve().parent

def need(ok,label):
    if not ok: raise RuntimeError(label)

def mat(n,code):
    z=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(itertools.combinations(range(n),2)):
        z[i][j]=1 if (code>>b)&1 else -1;z[j][i]=-z[i][j]
    return z

def codeof(z):
    return sum(int(z[i][j]==1)<<b for b,(i,j) in enumerate(itertools.combinations(range(len(z)),2)))

def curl(z):
    r=len(z)-1
    return tuple(z[i][j]+z[j][r]+z[r][i] for i,j in itertools.combinations(range(r),2))

def H(z):
    n=len(z)
    @lru_cache(None)
    def f(left,v):
        if not left:return 1
        return sum(f(left^(1<<w),w) for w in range(n) if left>>w&1 and z[v][w]==1)
    return sum(f(((1<<n)-1)^(1<<v),v) for v in range(n))

def components(z):
    # Tarjan, deliberately not the producers' transitive-closure algorithm.
    n=len(z); ix={}; low={};stack=[];on=set();found=[]
    def visit(v):
        ix[v]=low[v]=len(ix);stack.append(v);on.add(v)
        for w in range(n):
            if z[v][w]!=1:continue
            if w not in ix:visit(w);low[v]=min(low[v],low[w])
            elif w in on:low[v]=min(low[v],ix[w])
        if low[v]==ix[v]:
            block=[]
            while True:
                w=stack.pop();on.remove(w);block.append(w)
                if w==v:break
            found.append(tuple(sorted(block)))
    for v in range(n):
        if v not in ix:visit(v)
    return list(reversed(found))

def rotate(z,blocks,s):
    blocks=blocks[s:]+blocks[:s];rank={v:i for i,b in enumerate(blocks) for v in b}
    return [[z[v][w] if rank[v]==rank[w] else (1 if rank[v]<rank[w] else -1) for w in range(len(z))] for v in range(len(z))]

def assignment_masks(n):
    r=n-1; pairs=list(itertools.combinations(range(r),2));domains=[]
    for i,j in pairs:
        domains.append({v:sum(1<<a for a in range(1<<r) if abs(v+(1 if a>>i&1 else -1)-(1 if a>>j&1 else -1))==1) for v in (-3,-1,1,3)})
    return domains

def graph_connected(n,b):
    if n==1:return True
    adj=[0]*n
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        if b>>k&1:adj[i]|=1<<j;adj[j]|=1<<i
    seen=todo=1
    while todo:
        v=(todo&-todo).bit_length()-1;todo&=todo-1
        fresh=adj[v]&~seen;seen|=fresh;todo|=fresh
    return seen==(1<<n)-1

def curl_checks():
    rows=[]; total_H=[1]; class_H=[0]
    for n in range(1,7):
        groups=defaultdict(list);counts={}; comps={}
        for b in range(1<<(n*(n-1)//2)):
            z=mat(n,b);k=curl(z);groups[k].append(b)
            if k not in counts:counts[k]=H(z);comps[k]=len(components(z))
            else:need(H(z)==counts[k],('H fibre',n,b))
        for k,values in groups.items():
            z=mat(n,values[0]);blocks=components(z)
            need(sorted(codeof(rotate(z,blocks,s)) for s in range(len(blocks)))==values,('rotation fibre',n,k))
            need(len(values)==comps[k],('cardinality',n,k))
        domains=assignment_masks(n);valid=0;total=0;distribution=Counter()
        for k in itertools.product((-3,-1,1,3),repeat=len(domains)):
            feasible=(1<<(1<<(n-1)))-1
            for v,d in zip(k,domains):feasible &= d[v]
            size=feasible.bit_count()
            need(size==len(groups.get(k,())),('input assignment enumeration',n,k,size))
            if size:valid+=1;total+=size;distribution[size]+=1
        connected=sum(graph_connected(n,b) for b in range(1<<(n*(n-1)//2)))
        need(connected==valid,('connected graph identity',n))
        total_H.append(math.factorial(n)*2**((n-1)*(n-2)//2))
        class_H.append(total_H[n]-sum(math.comb(n-1,j-1)*class_H[j]*total_H[n-j] for j in range(1,n)))
        need(sum(counts.values())==class_H[n],('Hamiltonian weighted logarithm',n))
        row=dict(n=n,tournaments=total,curl_classes=valid,connected_graphs=connected,sum_H_over_curl_classes=class_H[n],all_alphabet_arrays=4**len(domains),fibre_distribution=dict(sorted(distribution.items())))
        rows.append(row);print('curl',row,flush=True)
    return rows

def Fprofile(z):
    n=len(z); f=[0]*(n+1)
    for p in itertools.permutations(range(n)):
        bad=sum(z[a][b]!=1 for a,b in zip(p,p[1:]))
        for k in range(bad+1,n+1):f[k]+=math.comb(n-1-bad,k-1-bad)
    need(all(f[k]%math.factorial(k)==0 for k in range(1,n+1)),('cover factorial divisibility',codeof(z)))
    return f

def context_weights(outside,sides,maxd):
    q=len(outside); weights=[0]*(maxd+1)
    for d in range(1,maxd+1):
        for pos in itertools.combinations(range(q+d),d):
            if any(y==x+1 for x,y in zip(pos,pos[1:])):continue
            pos=set(pos)
            for p in itertools.permutations(range(q)):
                it=iter(p);word=[-1 if j in pos else next(it) for j in range(q+d)]
                def edge(a,b):
                    if a==-1:return sides[b]
                    if b==-1:return not sides[a]
                    return outside[a][b]==1
                if all(edge(a,b) for a,b in zip(word,word[1:])):weights[d]+=1
    return weights

def extend(z,outside,sides):
    n=len(z);q=len(outside);full=[[0]*(n+q) for _ in range(n+q)]
    for i in range(n):
        for j in range(n):full[i][j]=z[i][j]
    for i in range(q):
        for j in range(q):full[n+i][n+j]=outside[i][j]
        for v in range(n):full[v][n+i]=1 if sides[i] else -1;full[n+i][v]=-full[v][n+i]
    return full

def cyclic(z,m):
    if m==0:return z
    out=mat(2*m,(1<<(m*(2*m-1)))-1)
    # All first block -> second block; both internally transitive.
    return extend(z,out,[True]*m+[False]*m)

@lru_cache(None)
def S(n,k):
    if not n:return int(not k)
    if k<1 or k>n:return 0
    return k*S(n-1,k)+S(n-1,k-1)

def cyclic_weights(m,n):
    if m==0:return [0,1]+[0]*(n-1)
    def g(k):return math.factorial(k)*S(m,k) if k>=0 else 0
    return [0]+[3*g(k)**2+2*g(k)*g(k+1)+g(k+1)**2+g(k-1)**2+2*g(k-1)*g(k) for k in range(1,n+1)]

def polyQ(f):
    n=len(f)-1;out=[0]*(n+1);fall=[1]
    for d in range(1,n+1):
        nxt=[0]*(len(fall)+1)
        for j,x in enumerate(fall):nxt[j+1]+=x;nxt[j]-=(d-1)*x
        fall=nxt
        for j,x in enumerate(fall):out[j]+=Fraction(f[d],math.factorial(d))*x
    return out

def det(a):
    z=[list(map(Fraction,row)) for row in a];out=1
    for j in range(len(z)):
        p=next((i for i in range(j,len(z)) if z[i][j]),None)
        if p is None:return 0
        if p!=j:z[p],z[j]=z[j],z[p];out=-out
        pivot=z[j][j];out*=pivot
        for i in range(j+1,len(z)):
            f=z[i][j]/pivot
            for k in range(j+1,len(z)):z[i][k]-=f*z[j][k]
    return out

def context_checks():
    generic=0;cyc=0;reconstruction=0;rank=[]
    for n in range(1,6):
        for b in range(1<<(n*(n-1)//2)):
            z=mat(n,b);f=Fprofile(z);q=polyQ(f);r=(n+1)//2
            need(q[n]==1 and all(not q[k] for k in range(n+1) if k%2!=n%2) and q[0]==0,('parity monicity',n,b))
            for m in range(max(0,r-1)):
                w=cyclic_weights(m,n);actual=H(cyclic(z,m));need(actual==sum(x*y for x,y in zip(w,f)),('cyclic formula',n,b,m));cyc+=1
            if n<=4:
                for k in range(4):
                    for ob in range(1<<(k*(k-1)//2)):
                        outside=mat(k,ob)
                        for side in itertools.product((False,True),repeat=k):
                            w=context_weights(outside,side,n)
                            actual=H(extend(z,outside,side));expected=sum(x*y for x,y in zip(w,f))
                            need(actual==expected,('general module context',n,b,k,ob,side));generic+=1
        print('context order',n,'generic',generic,'cyclic',cyc,flush=True)
    for n in range(1,9):
        r=(n+1)//2;rows=[]
        for j in range(r):
            z=mat(n,(1<<(n*(n-1)//2))-1)
            if j:z[0][2*j]=-1;z[2*j][0]=1
            f=Fprofile(z);q=polyQ(f);expected=[0]*(n+1);expected[n]=1
            for ell in range(1,j+1):expected[n-2*ell]=2*math.comb(2*j-1,2*ell-1)
            need(q==expected,('near transitive Q',n,j))
            rows.append([sum(x*y for x,y in zip(cyclic_weights(m,n),f)) for m in range(r)])
        d=det(rows);need(d!=0,('connection rank',n));rank.append(dict(n=n,rank=r,determinant=str(d)))
    return dict(arbitrary_module_contexts=generic,cyclic_contexts=cyc,near_transitive_rank_witnesses=rank)

if __name__=='__main__':
    report={'curl':curl_checks(),'contexts':context_checks(),'status':'PASS'}
    (OUT/'independent_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
