"""Independent class-signature incidence forest certificates for mirror constraints."""
import collections,csv,itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
D=[]
for row in sorted(csv.DictReader((P/'classes_input_n8.tsv').open(),delimiter='\t'),key=lambda r:int(r['canon'])):
 code=int(row['canon']);out=[0]*8
 for k,(u,v) in enumerate(itertools.combinations(range(8),2)):
  if code>>(27-k)&1:out[u]|=1<<v
  else:out[v]|=1<<u
 inc=[255^(1<<u)^out[u] for u in range(8)];e=[[0]*8 for _ in range(8)]
 for u in range(8):
  for v in range(8):
   if out[u]>>v&1:e[u][v]=e[v][u]=(out[v]&inc[u]).bit_count()
 c=[sum(e[u][v] for v in range(8) if out[u]>>v&1)for u in range(8)]
 D.append((code,out,inc,e,c))
def signature(d,p,q,rev=False):
 code,out,inc,e,c=d;arr=inc if rev else out
 return tuple(sorted(tuple(sorted(((c[v]%p if p else c[v]),(e[u][v]%q if q else e[u][v]))for v in range(8)if arr[u]>>v&1))for u in range(8)))
O={}
for p in range(1,9):
 for q in range(1,8):
  pairs=[(signature(d,p,q),signature(d,p,q,True))for d in D];left=collections.Counter(a for a,b in pairs);right=collections.Counter(b for a,b in pairs)
  par={};cycles=0
  def find(x):
   if x not in par:par[x]=x
   if par[x]!=x:par[x]=find(par[x])
   return par[x]
  for a,b in pairs:
   x=find((0,a));y=find((1,b))
   if x==y:cycles+=1
   else:par[x]=y
  bad=[D[i][0]for i,(a,b)in enumerate(pairs)if left[a]>1 and right[b]>1]
  O[f'{p},{q}']={'alphabet_size':p*q,'direct_signatures':len(left),'joint_signatures':len(set(pairs)),'singleton_or_converse':not bad,'both_nonsingleton_classes':len(bad),'class_incidence_cycle_rank':cycles,'forest_sufficient':cycles==0,'first_bad_canons':bad[:12]}
R={'status':'FINITE-EXACT independent sufficient-condition scan','menus':O,'minimal_forest_alphabet':min(v['alphabet_size']for v in O.values()if v['forest_sufficient'])}
(P/'converse_menu_scan_receipt.json').write_text(json.dumps(R,indent=2)+'\n')
print(json.dumps({k:v for k,v in O.items()if v['alphabet_size']<=12},indent=2))
