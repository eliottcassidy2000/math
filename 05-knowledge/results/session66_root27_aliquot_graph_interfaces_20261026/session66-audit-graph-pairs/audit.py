#!/usr/bin/env python3
"""Independent exact audit. Standard library; no producer routines imported.
Generic graph isomorphism first quotients identical adjacency rows (false twins),
then enumerates weighted quotient isomorphisms. Strict decoders reconstruct every
edge and reject malformed images. Gates remain live under python -O.
"""
from itertools import combinations, permutations, product
from collections import Counter, defaultdict
from functools import lru_cache
from math import factorial, comb
import hashlib, json, os, pathlib, random, subprocess

HERE=pathlib.Path(__file__).resolve().parent
PROD=pathlib.Path(os.environ.get('S66_PRODUCER',str(HERE.parent/'session66-graph-pairs')))
REPO=pathlib.Path(os.environ['S66_REPO']) if os.environ.get('S66_REPO') else None

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def orient(n,bits):
    a=[[0]*n for _ in range(n)]
    for k,(u,v) in enumerate(combinations(range(n),2)):
        a[u][v]=int(bool(bits&(1<<k)));a[v][u]=1-a[u][v]
    return a

def valid_tournament(a):
    n=len(a)
    return all(len(r)==n for r in a) and all(a[i][i]==0 for i in range(n)) and all(a[i][j] in (0,1) and a[j][i]==1-a[i][j] for i,j in combinations(range(n),2))

def transport(a,p):
    z=[[0]*len(a) for _ in a]
    for u,r in enumerate(a):
        for v,x in enumerate(r):z[p[u]][p[v]]=x
    return z

def signature(a):
    n=len(a)
    return min(tuple(a[p[i]][p[j]] for i,j in combinations(range(n),2)) for p in permutations(range(n)))

def autos(a):
    return [p for p in permutations(range(len(a))) if all(a[u][v]==a[p[u]][p[v]] for u in range(len(a)) for v in range(len(a)))]

def make(a,kind):
    n=len(a);width=7 if kind=='joint' else 9
    edges={}
    def edge(u,v,k=1):
        u,v=sorted((u,v));edges[u,v]=edges.get((u,v),0)|k
    for v in range(n):
        base=width*v
        if kind=='joint':
            for path in ([0,1],[0,2,3],[0,4,5,6]):
                for u,w in zip(path,path[1:]):edge(base+u,base+w)
        else:
            edge(base,base+1);edge(base,base+2)
            for core,leafbank in ((0,[3]),(1,[4,5]),(2,[6,7,8])):
                for leaf in leafbank:edge(base+core,base+leaf)
    for u in range(n):
        for v in range(n):
            if a[u][v]:edge(width*u+1,width*v+(3 if kind=='joint' else 2),2 if kind=='joint' else 1)
    g=[[0]*(n*width) for _ in range(n*width)]
    for (u,v),c in edges.items():g[u][v]=g[v][u]=c
    return g

def matrix_ok(g,maxcolor):
    N=len(g)
    return all(len(r)==N and all(isinstance(x,int) for x in r) for r in g) and all(g[i][i]==0 for i in range(N)) and all(g[i][j]==g[j][i] and g[i][j] in range(maxcolor+1) for i,j in combinations(range(N),2))

def components(adj):
    todo=set(range(len(adj)));cc=[]
    while todo:
        s={min(todo)};front=list(s)
        while front:
            u=front.pop()
            for v in adj[u]-s:s.add(v);front.append(v)
        todo-=s;cc.append(s)
    return cc

def strict_decode(g,kind):
    """Return recovered tournament, or ValueError exactly off the codec image."""
    def reject(cond):
        if not cond:raise ValueError('outside codec image')
    reject(matrix_ok(g,3 if kind=='joint' else 1))
    N=len(g);mapping=[];tails=[];heads=[]
    if kind=='joint':
        adj=[{v for v in range(N) if g[u][v]&1} for u in range(N)]
        for cc in components(adj):
            reject(len(cc)==7)
            candidates=[u for u in cc if len(adj[u])==3];reject(len(candidates)==1)
            c=candidates[0];arms={};seen={c}
            for v in sorted(adj[c]):
                path=[];prev=c
                while True:
                    reject(v not in seen);seen.add(v);path.append(v)
                    nxt=adj[v]-{prev};reject(len(nxt)<=1)
                    if not nxt:break
                    prev,v=v,next(iter(nxt))
                reject(len(path) not in arms);arms[len(path)]=path
            reject(seen==cc and set(arms)=={1,2,3})
            loc=[c]+arms[1]+arms[2]+arms[3];mapping+=loc
            tails.append(arms[1][0]);heads.append(arms[2][1])
        a=[[int(bool(g[t][h]&2)) for h in heads] for t in tails]
    else:
        adj=[{v for v in range(N) if g[u][v]} for u in range(N)]
        leaves={u for u in range(N) if len(adj[u])==1};cores=set(range(N))-leaves
        reject(all(len(adj[u])>=3 for u in cores))
        typ={u:len(adj[u]&leaves) for u in cores};reject(all(t in (1,2,3) for t in typ.values()))
        for c in sorted(u for u in cores if typ[u]==1):
            ns=adj[c]&cores;reject(len(ns)==2)
            ts=[u for u in ns if typ[u]==2];hs=[u for u in ns if typ[u]==3]
            reject(len(ts)==len(hs)==1);t,h=ts[0],hs[0]
            tails.append(t);heads.append(h)
            mapping += [c,t,h]+sorted(adj[c]&leaves)+sorted(adj[t]&leaves)+sorted(adj[h]&leaves)
        a=[[int(bool(g[t][h])) for h in heads] for t in tails]
    reject(len(mapping)==N and set(mapping)==set(range(N)))
    reject(valid_tournament(a))
    expected=make(a,kind)
    reject(all(expected[i][j]==g[mapping[i]][mapping[j]] for i in range(N) for j in range(N)))
    return a

def twin_quotient(g):
    banks=defaultdict(list)
    for u,r in enumerate(g):banks[tuple(r)].append(u)
    groups=list(banks.values());w=[len(v) for v in groups]
    q=[[g[a[0]][b[0]] for b in groups] for a in groups]
    return q,w

def graph_isomorphisms(g,h,count=False):
    """Generic, untyped, edge-colored isomorphism; independent of both decoders."""
    if len(g)!=len(h):return 0
    q,w=twin_quotient(g);r,z=twin_quotient(h);n=len(q)
    if len(r)!=n:return 0
    def invariants(a,weights):
        return [(weights[u],tuple(sorted((a[u][v],weights[v]) for v in range(n) if a[u][v]))) for u in range(n)]
    x,y=invariants(q,w),invariants(r,z)
    if Counter(x)!=Counter(y):return 0
    domains=[{v for v in range(n) if x[u]==y[v]} for u in range(n)]
    def search(dom,assignment):
        if len(assignment)==n:return 1
        u=min((u for u in range(n) if u not in assignment),key=lambda u:(len(dom[u]),-sum(bool(c) for c in q[u])))
        total=0
        for v in sorted(dom[u]):
            new={k:s.copy() for k,s in dom.items()};new[u]={v};good=True
            for k in range(n):
                if k==u or k in assignment:continue
                new[k]={j for j in dom[k] if j!=v and q[u][k]==r[v][j]}
                if not new[k]:good=False;break
            if good:
                got=search(new,{**assignment,u:v});total+=got
                if got and not count:return 1
        return total
    value=search(dict(enumerate(domains)),{})
    return value*__import__('functools').reduce(lambda a,b:a*factorial(b),w,1) if count else int(bool(value))

def literal_cycles(a):
    ans=[];n=len(a)
    # DFS rooted at the smallest cycle vertex is independent from producer's permutation loop.
    def extend(path):
        start=path[0];last=path[-1]
        if len(path)>=3 and len(path)%2 and a[last][start]:ans.append(tuple(path))
        for v in range(start+1,n):
            if v not in path and a[last][v]:extend(path+[v])
    for u in range(n):extend([u])
    return ans

def conflict(cs):
    return [[int(i!=j and bool(set(c)&set(d))) for j,d in enumerate(cs)] for i,c in enumerate(cs)]

def ind2(g):
    @lru_cache(None)
    def rec(vertices):
        if not vertices:return 1
        u,*others=vertices
        return rec(tuple(others))+2*rec(tuple(v for v in others if not g[u][v]))
    return rec(tuple(range(len(g))))

def hp(a):
    n=len(a)
    @lru_cache(None)
    def suffix(v,unused):
        if not unused:return 1
        return sum(suffix(w,unused^(1<<w)) for w in range(n) if unused&(1<<w) and a[v][w])
    return sum(suffix(v,((1<<n)-1)^(1<<v)) for v in range(n)) if n else 1

def clone_audit(a):
    n=len(a)
    j=[[int(u!=v and all(a[u][x]==a[v][x] for x in range(n) if x not in (u,v))) for v in range(n)] for u in range(n)]
    cc=components([{v for v in range(n) if j[u][v]} for u in range(n)])
    # Separate exhaustive subset definition of every transitive module.
    modules=[]
    for bits in range(1,1<<n):
        s={u for u in range(n) if bits&(1<<u)}
        uniform=all(len({a[u][v] for u in s})==1 for v in set(range(n))-s)
        trianglefree=all(not(a[u][v] and a[v][w] and a[w][u]) for u,v,w in permutations(s,3))
        if uniform and trianglefree:modules.append(s)
    maximal=[s for s in modules if not any(s<t for t in modules)]
    require({frozenset(s) for s in cc}=={frozenset(s) for s in maximal},'maximal transitive modules')
    order=[sorted(s,key=lambda u:-sum(a[u][v] for v in s)) for s in cc]
    q=[[a[s[0]][t[0]] if i!=j else 0 for j,t in enumerate(order)] for i,s in enumerate(order)]
    block=[i for i,s in enumerate(order) for _ in s];flat=[u for s in order for u in s]
    reconstructed=[[int(x<y) if block[x]==block[y] else q[block[x]][block[y]] for y in range(n)] for x in range(n)]
    require(reconstructed==[[a[u][v] for v in flat] for u in flat],'clone reconstruction')
    return q,[len(s) for s in cc],cc

def even_cycles(p):
    remain=set(range(len(p)));total=0
    while remain:
        u=min(remain);v=u;k=0
        while v in remain:remain.remove(v);k+=1;v=p[v]
        total+=not(k%2)
    return total

def run():
    out={'scope':'Independent finite exact tests plus separately recorded all-order proof audit; no producer code imported.'}
    pins=json.loads((PROD/'SOURCE_PINS.json').read_text());verified=[]
    for f in pins['files']:
        data=(PROD/f['local_path']).read_bytes()
        require(hashlib.sha256(data).hexdigest()==f['sha256'],'sha pin '+f['local_path'])
        blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        require(blob==f['git_blob'],'git blob '+f['local_path'])
        if REPO is not None:
            raw=subprocess.check_output(['git','-C',str(REPO),'show',pins['source_commit']+':'+f['source_path']])
            require(data==raw,'source pin bytes '+f['source_path'])
        verified.append(f['source_path'])
    out['source_pins_verified']=len(verified)
    rng=random.Random(66072101);reps={};roundtrips=0;clonecases=0;ocfcases=0;multi=None
    for n in range(6):
        reps[n]={}
        for bits in range(1<<comb(n,2)):
            a=orient(n,bits);s=signature(a);reps[n].setdefault(s,a)
            for kind in ('joint','gamma'):
                g=make(a,kind);require(strict_decode(g,kind)==a,'strict unlabelled roundtrip')
                p=list(range(len(g)));rng.shuffle(p);b=strict_decode(transport(g,p),kind)
                require(signature(b)==s,'arbitrarily relabelled strict roundtrip');roundtrips+=2
            clone_audit(a);clonecases+=1
            cs=literal_cycles(a);require(hp(a)==ind2(conflict(cs)),'literal OCF');ocfcases+=1
            bysupport=Counter(frozenset(c) for c in cs)
            if multi is None and max(bysupport.values(),default=0)>1:
                unique=list(bysupport);bad=ind2(conflict(unique))
                require(bad!=hp(a),'support-only hostile must differ')
                multi={'n':n,'bits':bits,'cycles':cs,'support_multiplicities':sorted(bysupport.values()),'H':hp(a),'support_only_I2':bad}
    require([len(reps[n]) for n in reps]==[1,1,1,2,4,12],'class counts')
    out['all_labelled_0_to_5']={'tournaments':clonecases,'strict_roundtrip_checks':roundtrips,'clones_against_all_transitive_modules':clonecases,'OCF_DFS_vs_path_recursion':ocfcases,'class_counts':[len(reps[n]) for n in reps]}
    pairtests=0;scrambled=0;automorphism_cases=[]
    for n in reps:
        aa=list(reps[n].values())
        for kind in ('joint','gamma'):
            for i,a in enumerate(aa):
                g=make(a,kind)
                scramble=list(range(len(g)));rng.shuffle(scramble)
                require(graph_isomorphisms(g,transport(g,scramble)),'generic scrambled positive');scrambled+=1
                for j in range(i,len(aa)):
                    require(bool(graph_isomorphisms(g,make(aa[j],kind)))==(i==j),'generic codec isomorphism iff');pairtests+=1
                observed=graph_isomorphisms(g,g,True);expected=len(autos(a))*(12**n if kind=='gamma' else 1)
                require(observed==expected,'generic full automorphism count')
                automorphism_cases.append({'n':n,'kind':kind,'tournament_key':list(signature(a)),'count':observed})
    out['generic_graph_iso']={'class_pairs':pairtests,'scrambled_positive_cases':scrambled,'automorphism_cases':len(automorphism_cases),'counts':automorphism_cases}
    out['literal_cycle_support_hostile']=multi
    malformed=Counter();valid_reversals=0
    for n in (1,2):
        a=orient(n,0)
        for kind in ('joint','gamma'):
            base=make(a,kind);N=len(base)
            for layer in ((1,2) if kind=='joint' else (1,)):
                for u,v in combinations(range(N),2):
                    g=[row[:] for row in base];g[u][v]^=layer;g[v][u]^=layer
                    try:strict_decode(g,kind)
                    except ValueError:malformed[kind]+=1
                    else:raise RuntimeError('single-edge malformed image accepted')
            g=[r+[0] for r in base]+[[0]*(N+1)]
            try:strict_decode(g,kind)
            except ValueError:malformed[kind]+=1
            else:raise RuntimeError('extra isolate accepted')
            if n==2:
                b=orient(2,1);require(strict_decode(make(b,kind),kind)==b,'valid arc reversal rejected');valid_reversals+=1
    out['malformed_rejections']={**dict(malformed),'valid_two_edge_arc_reversal_controls':valid_reversals}
    defects=0;minima=0
    for n in range(5):
        for p in permutations(range(n)):
            best=comb(n,2)+1
            for bits in range(1<<comb(n,2)):
                a=orient(n,bits);d=sum(a[u][v]!=a[p[u]][p[v]] for u,v in combinations(range(n),2));best=min(best,d)
                for kind,width in (('joint',7),('gamma',9)):
                    g=make(a,kind);lift=[width*p[u//width]+u%width for u in range(len(g))]
                    delta=sum(g[u][v]!=g[lift[u]][lift[v]] for u,v in combinations(range(len(g)),2))
                    require(delta==2*d,'edge discrepancy');defects+=1
            require(best==even_cycles(p),'minimum equals even vertex cycles');minima+=1
    out['defects']={'lifted_metric_cases':defects,'full_permutation_minimum_cases':minima}
    # Known natural-pair and order-forgetting hostiles.
    t=[[int(i<j) for j in range(3)] for i in range(3)]
    x=[r[:] for r in t];x[0][1]=0;x[1][0]=1
    y=[r[:] for r in t];y[0][2]=0;y[2][0]=1
    def fb(a):return [[0 if u==v else (1 if a[min(u,v)][max(u,v)] else 2) for v in range(len(a))] for u in range(len(a))]
    require(hp(x)==1 and hp(y)==3 and graph_isomorphisms(fb(x),fb(y)),'lost order hostile')
    a=[[0,1,1,1],[0,0,1,0],[0,0,0,1],[0,1,0,0]];b=[list(r) for r in zip(*a)]
    ga=make(a,'joint');gb=make(b,'joint')
    require(signature(a)!=signature(b),'source/converse distinct')
    for layer in (1,2):require(graph_isomorphisms([[int(bool(x&layer)) for x in r] for r in ga],[[int(bool(x&layer)) for x in r] for r in gb]),'independent layers coincide')
    require(not graph_isomorphisms(ga,gb) and not graph_isomorphisms(make(a,'gamma'),make(b,'gamma')),'faithful codes must separate converse hostile')
    invariant=[]
    for mask in range(8):
        g=[[0]*3 for _ in range(3)]
        for k,(u,v) in enumerate(combinations(range(3),2)):g[u][v]=g[v][u]=(mask>>k)&1
        if transport(g,[1,2,0])==g:invariant.append(mask)
    require(invariant==[0,7],'C3 no-go')
    require(ind2([[int(i!=j) for j in range(3)] for i in range(3)])==7,'K3=7')
    require(ind2([[int(abs(i-j)==1) for j in range(4)] for i in range(4)])==21,'P4=21')
    out['hostiles']={'lost_order_n':3,'lost_alignment_n':4,'C3_invariant_graph_masks':invariant,'I_K3_2':7,'I_P4_2':21}
    # Canonical module automorphism quotient, through all small classes.
    cloneaut=0
    for dd in reps.values():
        for a in dd.values():
            q,weights,cc=clone_audit(a)
            qa=[p for p in autos(q) if all(weights[i]==weights[p[i]] for i in range(len(q)))]
            require(len(qa)==len(autos(a)),'weighted quotient automorphism group')
            cloneaut+=1
    out['weighted_quotient_automorphism_cases']=cloneaut
    # C3[TT2,TT2,TT2] gives genuine gap action and unsafe labelled toggle.
    base=[[0,1,0],[0,0,1],[1,0,0]]
    inflated=[[int(u%2<v%2) if u//2==v//2 else base[u//2][v//2] for v in range(6)] for u in range(6)]
    q,w,cc=clone_audit(inflated);aut=autos(inflated)
    require(w==[2,2,2] and len(aut)==3,'genuine three-gap action')
    gap_actions=[tuple(p[2*i]//2 for i in range(3)) for p in aut]
    def moved(s,p):return sum(1<<p[i] for i in range(3) if s&(1<<i))
    canon=lambda s:min(moved(s,p) for p in gap_actions)
    subsets={canon(s) for s in range(8)};poly=[sum(s.bit_count()==i for s in subsets) for i in range(4)]
    legal=[f for f in range(8) if all(canon(s)!=canon(t) or canon(s^f)==canon(t^f) for s in range(8) for t in range(8))]
    require(poly==[1,1,1,1] and legal==[0,7],'labelled gap does not descend')
    out['genuine_clone_gap_control']={'vertices':6,'automorphisms':3,'module_sizes':w,'subset_orbit_polynomial':poly,'legal_fixed_toggle_masks':legal}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__':run()
