"""Independent exact audit. Standard library only; no producer imports/data.
Arc encoding: lexicographic pair bit 1 means lower label beats higher.
"""
import itertools as it
import json

def require(condition):
    if not condition:
        raise RuntimeError("Independent mathematical audit check failed")

def tournament(n,code):
    edges=set()
    for bit,(u,v) in enumerate(it.combinations(range(n),2)):
        edges.add((u,v) if code & (1<<bit) else (v,u))
    return edges

def paths(n,E):
    # Enumerate every directed simple path recursively, not subset counting DP.
    out=[]
    def visit(p,mask):
        out.append((p,mask))
        for v in range(n):
            if not mask>>v&1 and (p[-1],v) in E: visit(p+(v,),mask|1<<v)
    for v in range(n): visit((v,),1<<v)
    return out

def data(n,E):
    ps=paths(n,E); h=[0]*(1<<n); h[0]=1
    for p,m in ps: h[m]+=1
    K=[[0]*n for _ in range(n)]; full=(1<<n)-1
    cycles=[]
    for p,m in ps:
        if len(p)%2==0: K[p[0]][p[-1]]+=h[full^m]
        elif len(p)>=3 and p[0]==min(p) and (p[-1],p[0]) in E: cycles.append(m)
    # Enumerate families of directed cycles. Equal supports can encode distinct cycles.
    atoms=[0]*(1<<n)
    def pack(i,used,weight):
        atoms[used]+=weight
        for j in range(i,len(cycles)):
            if not used & cycles[j]: pack(j+1,used|cycles[j],2*weight)
    pack(0,0,1)
    for S in range(1<<n):
        require(h[S]==sum(w for U,w in enumerate(atoms) if U & S==U))
    return h,K,atoms

def cut(n,E,S): return sum(1 for u,v in E if S>>u&1 and not S>>v&1)
def semidegree(n,E): return min(min(sum(u==x for u,v in E),sum(v==x for u,v in E)) for x in range(n))
def direct_h(n,E): return sum(all((p[i],p[i+1]) in E for i in range(n-1)) for p in it.permutations(range(n)))

report={'bases':0,'attachments':0,'mixed_differences':0,'curvature_faces':0,'mincut_birth_certificates':0,'connectivity_tournaments':0}
Sprev={1}
for n in range(1,6):
    records=[]; spectrum=set()
    for code in range(1<<(n*(n-1)//2)):
        E=tournament(n,code); h,K,atoms=data(n,E); V=(1<<n)-1
        spectrum.add(h[V]); values=[]
        for A in range(1<<n):
            F=E|{(n,u) if A>>u&1 else (u,n) for u in range(n)}
            val=direct_h(n+1,F); values.append(val)
            c=sum(K[u][v] for u in range(n) for v in range(n) if A>>u&1 and not A>>v&1)
            require(val==h[V]+2*c and val not in (7,21))
            require((c==0)==(cut(n,E,A)==0))
            report['attachments']+=1
        for u,v in it.combinations(range(n),2):
            require(K[u][v]+K[v][u]>=h[V^(1<<u)^(1<<v)]>=1)
            for A in range(1<<n):
                if A & ((1<<u)|(1<<v)): continue
                d=values[A|1<<u|1<<v]-values[A|1<<u]-values[A|1<<v]+values[A]
                require(d==-2*(K[u][v]+K[v][u]))
                report['curvature_faces']+=1
        for B in range(1<<n):
            for R in range(1,1<<n):
                if B&R: continue
                lhs=sum((-1)**(R.bit_count()-Q.bit_count())*h[B|Q] for Q in range(1<<n) if Q&R==Q)
                rhs=sum(w for U,w in enumerate(atoms) if U&R==R and U&(B|R)==U)
                require(lhs==rhs>=0)
                report['mixed_differences']+=1
        if n>=2 and all(cut(n,E,A)>0 for A in range(1,V)):
            records.append((h[V],min(values[1:V]),values[1:V]))
        report['bases']+=1
    for base,smallest,values in records:
        require(smallest>base)
        if smallest>max(spectrum):
            require(all(v not in spectrum for v in values))
            report['mincut_birth_certificates']+=1
    require(Sprev<=spectrum)
    Sprev=spectrum
report['residual_words']={}; total=0
for m in range(2,7):
    vals=set()
    for inner in it.product((0,1),repeat=m-2):
        w=(1,)+inner+(0,)
        val=sum((1 if j-i==1 else 2**(j-i-2)) for i in range(m) for j in range(i+1,m) if w[i] and not w[j])
        require(val not in (3,10)); vals.add(val); total+=1
    report['residual_words'][m]=sorted(vals)
require(total==31)
report['residual_word_count']=total
for n in range(2,7):
    for code in range(1<<(n*(n-1)//2)):
        E=tournament(n,code); d=semidegree(n,E)
        lam=min(cut(n,E,S) for S in range(1,(1<<n)-1))
        for k in range(1,n):
            if n<=4*k+1: require((d>=k)==(lam>=k))
        report['connectivity_tournaments']+=1
report['sharp_examples']=[]
for k in range(1,5):
    q=2*k+1; n=2*q; E=set()
    for offset in (0,q):
        E|={(offset+i,offset+(i+d)%q) for i in range(q) for d in range(1,k+1)}
    E|={(i,j) for i in range(q) for j in range(q,n)}
    if k>=2: E.remove((0,q)); E.add((q,0))
    require(semidegree(n,E)==k)
    # Two strong blocks, one backward bridge: proof gives lambda=1 (k>=2).
    r=cut(n,E,((1<<q)-1)<<q)
    require(r==(1 if k>=2 else 0))
    report['sharp_examples'].append({'k':k,'n':n,'semidegree':k,'backward_cut':r})
print(json.dumps(report,indent=2,sort_keys=True))
