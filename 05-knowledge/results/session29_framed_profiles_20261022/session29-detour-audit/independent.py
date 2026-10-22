import itertools as it, json, collections, pathlib, hashlib
D=pathlib.Path('/workspace/shared/session29-bridge-extension')
R=json.loads((D/'results.json').read_text()); Z=json.loads((D/'detours-results.json').read_text())
def check(x,s):
 if not x: raise ValueError(s)
def edges(n): return list(it.combinations(range(n),2))
def arc(m,n,u,v):
 if u==v:return False
 k=edges(n).index(tuple(sorted((u,v))))
 return bool((m>>k)&1) == (u<v)
def encode(n,fun):return sum(int(fun(u,v))<<i for i,(u,v) in enumerate(edges(n)))
def induced(m,n,vs):return encode(len(vs),lambda u,v:arc(m,n,vs[u],vs[v]))
def canon(m,n,root=None):
 vs=tuple(i for i in range(n) if i!=root)
 return min(induced(m,n,p+(() if root is None else (root,))) for p in it.permutations(vs))
def tau(m):return encode(5,lambda u,v:arc(m,5,4-v,4-u))
def square(base):
 return [(('S',tuple(sorted((a,c)))),('T',tuple(sorted((b,d))))) for a,b in base for c,d in base]
def bridges(es):
 adj=collections.defaultdict(list)
 for k,(u,v) in enumerate(es):adj[u].append((v,k));adj[v].append((u,k))
 depth={};low={};out=[]
 def visit(u,pe):
  depth[u]=low[u]=len(depth)
  for v,k in adj[u]:
   if k==pe:continue
   if v not in depth:
    visit(v,k);low[u]=min(low[u],low[v])
    if low[v]>depth[u]:out.append(k)
   else:low[u]=min(low[u],depth[v])
 for u in list(adj):
  if u not in depth:visit(u,-1)
 return sorted(out)
tiles=[m for m in range(1024) if m&657==657]
rows={m:canon(m,5) for m in tiles};cols={m:canon(m^366,5,4) for m in tiles}
check([{'tile':m,'source_class':rows[m],'root_class':cols[m]} for m in tiles]==R['base_incidence'],'base graph records')
for m in tiles:
 check(tau(tau(m))==m and tau(m) in tiles,'tau involution on spine tiles')
 check(canon(m^1023,5)==rows[tau(m)],'second source converse normalization')
 check(canon((m^366)^1023,5,0)==cols[tau(m)],'second rooted converse normalization')
base=[(rows[m],cols[m]) for m in tiles];E=square(base);bs=bridges(E)
rd=collections.Counter(rows.values());cd=collections.Counter(cols.values())
pred=[i*65 for i,m in enumerate(tiles) if rd[rows[m]]==1 or cd[cols[m]]==1]
check(bs==pred,'full bridge theorem')
check(len(bs)==19,'19 boundary bridges')
check(len(set(u for u,v in E))==78 and len(set(v for u,v in E))==561,'boundary vertex counts')
sourceleaf=[m for m in tiles if rd[rows[m]]==1]
check(sourceleaf==[767,991,1021,1023],'source leaves')
controls=0
for mult in it.product(range(3),repeat=4):
 b=[(i,j) for (i,j),k in zip(it.product(range(2),repeat=2),mult) for z in range(k)]
 r=collections.Counter(x for x,y in b);c=collections.Counter(y for x,y in b);n=len(b)
 p=[i*(n+1) for i,(x,y) in enumerate(b) if r[x]==1 or c[y]==1]
 check(bridges(square(b))==p,'small multigraph criterion');controls+=1
def address(w):
 e,f=w;f=tau(f)
 return ('S',tuple(sorted((rows[e],rows[f])))),('T',tuple(sorted((cols[e],cols[f]))))
sel=[m for m in tiles if rows[m] in (40,41)]
check(len(sel)==18,'selected alphabet')
words=list(it.product(sel,repeat=2));se=[address(w) for w in words];sb=bridges(se)
check(len(sb)==8,'eight selected boundary bridges')
repaired=set();local=[]
for rec in R['original_eight']:
 w=tuple(rec['actual_original_word']);check(words.index(w) in sb,'old unresolved edge')
 path=rec['alternate_actual_words']
 if path is None:continue
 start,end=address(w);cur=start;seen={cur}
 for q in path:
  check(tuple(q)!=w,'do not reuse removed edge')
  u,v=address(q);check(cur in (u,v),'path incidence');cur=v if cur==u else u
  check(cur not in seen,'simple path');seen.add(cur)
 check(cur==end,'path endpoint');repaired.add(w);local.append({'original':w,'alternate':path})
def target(w):
 a,b=w
 return encode(10,lambda u,v:arc(a^366,5,u,v) if v<5 else arc(b^366,5,u-5,v-5) if u>=5 else (u,v)==(4,5))
def flip(m,p):
 check(sorted(p)==list(range(10)),'Hamiltonian permutation')
 check(all(arc(m,10,u,v) for u,v in zip(p,p[1:])),'Hamiltonian arcs')
 keep={tuple(sorted(q)) for q in zip(p,p[1:])}
 return m ^ sum(1<<i for i,q in enumerate(edges(10)) if q not in keep)
def scc(m,n):
 reach=[[u==v or arc(m,n,u,v) for v in range(n)] for u in range(n)]
 for k in range(n):
  for u in range(n):
   for v in range(n):reach[u][v]|=reach[u][k] and reach[k][v]
 groups=[];unused=set(range(n))
 while unused:
  u=min(unused);g=sorted(v for v in unused if reach[u][v] and reach[v][u]);unused.difference_update(g);groups.append(g)
 snapshot=list(groups);groups.sort(key=lambda g:sum(arc(m,n,g[0],h[0]) for h in snapshot if g!=h),reverse=True)
 return groups
def hp(m):
 dp=[[0]*10 for _ in range(1024)]
 for v in range(10):dp[1<<v][v]=1
 for s in range(1,1024):
  for v in range(10):
   if not dp[s][v]:continue
   for w in range(10):
    if not s>>w&1 and arc(m,10,v,w):dp[s|1<<w][w]+=dp[s][v]
 return sum(dp[-1])
ct=encode(10,lambda u,v:arc(335,5,u,v) if v<5 else arc(842,5,v-5,u-5) if u>=5 else (u,v)==(4,9))
detail=[]
for z in Z:
 rec=next(x for x in R['original_eight'] if x['tile']==z['tile']);w=tuple(rec['actual_original_word'])
 T=target(w);check(T==z['target_mask'],'target literal')
 U=flip(T,z['target_to_U_path']);check(U==z['U_mask'],'first flip')
 V=flip(U,z['U_to_V_path']);check(V==z['V_mask'],'second flip')
 pp=z['V_to_canonical_permutation'];check(sorted(pp)==list(range(10)),'relabel permutation')
 inv=[pp.index(i) for i in range(10)];mm=induced(V,10,inv)
 if z['V_converse']:mm^=(1<<45)-1
 check(mm==ct,'mixed canonical target literal')
 gs=scc(U,10);key=[(len(g),canon(induced(U,10,g),len(g))) for g in gs]
 check(key==[(4,4),(6,81)],'U SCC canonical classes')
 check(len(scc(T,10))==len(scc(V,10))==1,'T and V strong')
 check(address(w)[0]==('S',(41,41)),'original BB source')
 check(address((689,661))==( ('S',(41,41)),('T',(335,842)) ),'return boundary edge')
 check(address(w)[1]!=('T',(335,842)),'distinct strong nodes')
 h=hp(T);check(h=={689:585,721:221}[z['tile']],'H count')
 repaired.add(w);detail.append({'word':w,'T':T,'U':U,'V':V,'H_T':h,'U_SCC_types':key,'old_to_new':pp,'converse':z['V_converse']})
check(repaired=={words[i] for i in sb},'all eight repairs')
result={'status':'PASS','implementation':'independent integer-bit adjacency, Tarjan, subset DP, explicit permutation checks','base_tiles':len(tiles),'base_rows':len(rd),'base_columns':len(cd),'boundary_edges':len(E),'boundary_bridges':len(bs),'source_leaf_tiles':sourceleaf,'remaining_unclassified_target_leaf_tiles':[tiles[i//65] for i in bs if tiles[i//65] not in sourceleaf+[689,721]],'small_multigraph_controls':controls,'selected_edges':len(se),'selected_intrinsic_nonbridges':len(se)-len(sb),'local_repairs':local,'literal_detours':detail,'global_nonbridges_certified':len(se)}
print(json.dumps(result,indent=2,sort_keys=True))
