"""Exact finite square-sum atlas tests; standard library only."""
from collections import Counter
from itertools import permutations, combinations
from math import isqrt
from pathlib import Path
import json

def check(p, msg):
    if not p: raise RuntimeError(msg)
def edge(a,b): return tuple(sorted((a,b)))
def edges(p): return {edge(a,b) for a,b in zip(p,p[1:])}
def square(a,b): return a != b and isqrt(a+b)**2 == a+b
N=25
adj={v:[w for w in range(1,N+1) if square(v,w)] for v in range(1,N+1)}
# Independent vertex DFS, complete because 18 is a leaf.
check(adj[18]==[7], 'root completeness')
paths=[]
def dfs(p,used):
    if len(p)==N: paths.append(tuple(p)); return
    for w in adj[p[-1]]:
        if w not in used: dfs(p+[w],used|{w})
dfs([18],{18})
E=list(map(edges,paths)); common=set.intersection(*E)
fa={v:[] for v in adj}
for a,b in common: fa[a].append(b); fa[b].append(a)
blocks=[];seen=set()
for v in range(1,N+1):
    if v in seen or len(fa[v])>1: continue
    b=[];prev=None;w=v
    while w is not None:
        b.append(w);seen.add(w);nxt=[x for x in fa[w] if x!=prev]
        prev,w=w,(nxt[0] if nxt else None)
    blocks.append(tuple(b))
check(seen==set(adj),'forest extraction')
owner={v:i for i,b in enumerate(blocks) for v in b}
def encode(p):
    out=[];i=0
    while i<len(p):
        j=owner[p[i]];b=blocks[j];seg=p[i:i+len(b)]
        check(tuple(seg) in (b,b[::-1]),'contiguous block')
        out.append((j,1 if tuple(seg)==b else -1));i+=len(b)
    return tuple(out)
codes=list(map(encode,paths))
# Separate contracted search, only endpoint seam arithmetic. No full paths supplied.
orient={i:[b] if len(b)==1 else [b,b[::-1]] for i,b in enumerate(blocks)}
contracted=[]
def glue(seq,used,p):
    if len(used)==len(blocks): contracted.append(tuple(p));return
    for i in range(len(blocks)):
        if i in used:continue
        for b in orient[i]:
            if square(p[-1],b[0]):glue(seq+[i],used|{i},p+list(b))
r=owner[18];root=next(b for b in orient[r] if b[0]==18)
glue([r],{r},list(root))
check(set(contracted)==set(paths),'contracted gluing mismatch')
# Exact inherited owner color definition, independent of square labels.
c1=[(23,13),(23,2),(25,11),(25,24),(15,1),(15,10),(17,8),(16,9)]
c2=[(13,3),(13,12),(14,2),(14,11),(19,6),(19,17),(20,5),(20,16),(21,4),(21,15),(22,3),(22,14),(24,1),(24,12)]
color={edge(*e):c for c,ee in [(1,c1),(2,c2)] for e in ee}
check(all(square(*e) for e in color),'invalid colored edge')
profiles=[tuple(sum(color.get(e,3)==c for e in ee) for c in (1,2,3)) for ee in E]
poly=Counter(profiles)
check(poly==Counter({(7,11,6):7,(7,10,7):2,(6,12,6):1}),'color polynomial')
# Every ambient signed adjacent-block swap / single block flip.
codeindex={s:i for i,s in enumerate(codes)}
swaps=[];flips=[]
for i,s in enumerate(codes):
    for j in range(len(s)-1):
        t=s[:j]+(s[j+1],s[j])+s[j+2:]
        if t in codeindex and i<codeindex[t]:swaps.append([i,codeindex[t],j])
    for j,(b,sign) in enumerate(s):
        if len(blocks[b])==1:continue
        t=s[:j]+((b,-sign),)+s[j+1:]
        if t in codeindex and i<codeindex[t]:flips.append([i,codeindex[t],j])
check(not swaps and not flips,'frozen atlas')
def components(ee):
    rem=set(range(len(paths)));cc=[]
    while rem:
        done={rem.pop()};todo=list(done)
        while todo:
            v=todo.pop()
            for a,b in ee:
                w=b if a==v else a if b==v else None
                if w in rem: rem.remove(w);done.add(w);todo.append(w)
        cc.append(sorted(done))
    return cc
exchange={}
for d in (1,2,3):
    ee=[(i,j) for i,j in combinations(range(len(paths)),2) if len(E[i]-E[j])<=d]
    exchange[d]={'edges':ee,'components':components(ee)}
check(sorted(map(len,exchange[1]['components']))==[4,6],'E1 census')
check(len(exchange[2]['components'])==1,'E2 connected')
# Reflection labels independently reconstruct all vertices in every path.
for p in paths:
    q=[a+b for a,b in zip(p,p[1:])];v=p[0]
    for j,t in enumerate(q):v=t-v;check(v==p[j+1],'reflection reconstruction')
# Weighted zonotope support law, all orders for 3-4-5 and small fine refinement.
def y(order,weights):
    out=[0]*len(weights);prior=0
    for i in order:out[i]=weights[i]*prior;prior+=weights[i]
    return out
w=(3,4,5);verts=[y(p,w) for p in permutations(range(3))]
for c in permutations((-2,0,3)):
    support=sum(w[i]*w[j]*max(c[i],c[j]) for i,j in combinations(range(3),2))
    check(max(sum(a*b for a,b in zip(c,v)) for v in verts)==support,'support identity')
# Refine the blocks into weights 1+2, 4, 2+3; contiguous fine orders must project.
f=(1,2,4,2,3);groups=((0,1),(2,),(3,4));groupof={a:i for i,g in enumerate(groups) for a in g}
const=[sum(f[a]*f[b] for a,b in combinations(g,2)) for g in groups]
refined_checked=0
for p in permutations(range(5)):
    runs=[]
    for a in p:
        i=groupof[a]
        if not runs or runs[-1]!=i:runs.append(i)
    if len(runs)!=3:continue
    fine=y(p,f);proj=[sum(fine[a] for a in g)-const[i] for i,g in enumerate(groups)]
    check(proj==y(runs,w),'refinement identity');refined_checked+=1
# Small hostile: same order geometry, opposite orientations yield different feasibility.
# Q15, split its known path after the first two vertices.
p15=(8,1,15,10,6,3,13,12,4,5,11,14,2,7,9)
check(all(square(a,b) for a,b in zip(p15,p15[1:])),'Q15 positive')
check(not square(8,15),'orientation-loss hostile')
# Exact color-respecting strong lumpability, for lazy uniform candidate surgery.
# Choose one of the other nine enumerated paths uniformly; reject if distance>d.
initial_cells=[]
for profile in dict.fromkeys(profiles):
    initial_cells.append([i for i,x in enumerate(profiles) if x==profile])
initial_cells.sort(key=min)
kernels={}
for d in (1,2):
    A=[[0]*len(paths) for _ in paths]
    for i,j in exchange[d]['edges']: A[i][j]=A[j][i]=1
    for i in range(len(paths)): A[i][i]=len(paths)-1-sum(A[i])
    kernels[d]=A

def refine(matrices):
    cells=[c[:] for c in initial_cells];history=[cells]
    while True:
        new=[]
        for C in cells:
            groups={}
            for i in C:
                signature=tuple(sum(A[i][j] for j in D) for A in matrices for D in cells)
                groups.setdefault(signature,[]).append(i)
            new.extend(groups.values())
        new.sort(key=min)
        if new==cells: return cells,history
        cells=new;history.append(cells)

markov={}
for d,A in kernels.items():
    cells,history=refine([A])
    Q=[[sum(A[C[0]][j] for j in D) for D in cells] for C in cells]
    for C in cells:
        check(all([sum(A[i][j] for j in D) for D in cells]==Q[cells.index(C)] for i in C),'strong lumpability')
    check(all(sum(row)==9 for row in Q),'stochastic quotient')
    # Symmetry gives stationary masses |C|/10 and exact detailed balance.
    check(all(len(C)*Q[a][b]==len(D)*Q[b][a] for a,C in enumerate(cells) for b,D in enumerate(cells)),'quotient reversibility')
    check(len(cells)==9,'coarsest state count')
    expected=[3,4] if d==1 else [3,6]
    check([C for C in cells if len(C)>1]==[expected],'merged fiber')
    markov[d]={'denominator':9,'initial_cells':initial_cells,'refinement_history':history,
      'coarsest_cells':cells,'quotient_numerator':Q,'stationary_mass_numerators':[len(C) for C in cells],
      'initial_color_transition_numerators':[[sum(A[i][j] for j in D) for D in initial_cells] for i in range(10)]}
common_cells,common_history=refine([kernels[1],kernels[2]])
check(len(common_cells)==10,'joint driver must retain all paths')
markov['joint']={'coarsest_cells':common_cells,'refinement_history':common_history}

receipt={'status':'FINITE-EXACT','n':N,'path_count':len(paths),'paths':paths,'common_edges':sorted(common),
 'blocks':blocks,'signed_codes':codes,'birth_color_profiles':profiles,
 'birth_color_polynomial':[{'exponents':k,'coefficient':v} for k,v in sorted(poly.items())],
 'markov_color_refinement':markov,'adjacent_block_swaps':swaps,'single_block_flips':flips,'exchange':exchange,
 'refinement_orders_checked':refined_checked,'positive_control':'Q15 concatenation; same projected order with first block reversed fails',
 'universes':'All Q25 paths by leaf-rooted DFS; all endpoint-compatible oriented orders of 11 common-forest blocks; all 120 fine orders for weights(1,2,4,2,3).'}
Path(__file__).with_name('certificate.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS: 10 Q25 paths; independent contracted search agrees; 11-block weighted permutohedral atlas frozen under adjacent swaps and individual flips; E1 components 6+4, E2 connected; color polynomial 7*x^7*y^11*z^6+2*x^7*y^10*z^7+x^6*y^12*z^6; refinement controls',refined_checked, '; Markov minimal color-respecting quotients: 9 states each, joint10')
