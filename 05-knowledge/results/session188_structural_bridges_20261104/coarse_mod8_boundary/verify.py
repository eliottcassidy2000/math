import importlib.util,itertools,math,json,hashlib
from pathlib import Path
s=importlib.util.spec_from_file_location('v',Path(__file__).resolve().parent.parent/'covered_response'/'verify.py');v=importlib.util.module_from_spec(s);s.loader.exec_module(v)
def num(G,a,c):
 A=[i for i in range(3) if a>>i&1];O=[i for i in range(3) if not a>>i&1];C=[i for i in range(3) if c>>i&1];P=[i for i in range(3) if not c>>i&1]
 def left(I):return v.add(v.R,*[G[u][w] for u in range(3) for w in I])
 def right(I):return v.add(v.R,*[G[u][w] for u in I for w in range(3)])
 def mid(I,J):return v.add(*[G[u][w] for u in I for w in J])
 return v.add(v.prod(left(A),v.add(v.R,mid(O,C)),right(P)),v.prod(left(C),mid(P,A),right(O)))
def functional(p,qs):
 z=0
 for k in range(4):
  for e,c in p.items():
   value=c*math.comb(k+2,2)
   for t,q in zip(e,qs):value*=q[k+t] if k+t<len(q) else 0
   z+=value
 return z
Gt=[list(x) for x in zip(*v.G)];D={(a,c):v.add(num(v.G,a,c),v.scale(num(Gt,a,c),-1)) for a,c in itertools.product(range(8),repeat=2)}

ROOT=Path(__file__).resolve().parent
certificates=[]
for (a,c),p in D.items():
 values=[]
 for mask in range(8):
  qs=[[0,1+2*((mask>>i)&1),-2+4*((mask>>i)&1),4] for i in range(3)]
  values.append(functional(p,qs))
 coefficients=[]
 for mask in range(8):
  total=0;sub=mask
  while True:
   total+=(-1)**((mask^sub).bit_count())*values[sub]
   if not sub:break
   sub=(sub-1)&mask
  v.need(total%8==0,'universal multilinear coefficient divisible by eight');coefficients.append(total)
 certificates.append({'block_a_mask':a,'block_c_mask':c,'parameter_values_at_boolean_masks':values,'multilinear_coefficients_by_subset_mask':coefficients})
def build(blocks,reverse=False):
 sizes=list(map(len,blocks));offsets=[sum(sizes[:i]) for i in range(3)];rows=[0]*sum(sizes)
 for i,A in enumerate(blocks):
  j=(i+(2 if reverse else 1))%3
  for u,r in enumerate(A):rows[offsets[i]+u]=(r<<offsets[i])|(((1<<sizes[j])-1)<<offsets[j])
 return rows,offsets
literal=[]
for blocks in ([v.tt(4)]*3,[v.tt(4),[8,1,3,6],v.tt(4)],[[8,1,3,6]]*3):
 T,offs=build(blocks);U,_=build(blocks,True);profiles=[v.profile(A) for A in blocks]
 for a,c in [(1,2),(1,3),(3,5)]:
  aa=sum(((1<<len(blocks[i]))-1)<<offs[i] for i in range(3) if a>>i&1);cc=sum(((1<<len(blocks[i]))-1)<<offs[i] for i in range(3) if c>>i&1)
  f=v.hp(v.extend(T,aa,cc));g=v.hp(v.extend(U,aa,cc));v.need((f-g)%8==0,'literal twelve-name blindness');v.need((f-g-functional(D[a,c],[[x%8 for x in q[:4]] for q in profiles]))%8==0,'literal word-functional agreement');literal.append({'child_rows':blocks,'child_profiles':profiles,'physical_query':[aa,cc],'block_query':[a,c],'values':[f,g]})
print(json.dumps({'status':'PASS','scope':'Unbounded block-uniform mod8 blindness for three arbitrary even-order children, each order at least four. Not an exact or unrestricted named-response collision.','symbolic_queries':64,'symbolic_coefficients':512,'coefficient_records':certificates,'literal_controls':literal,'source_sha256':{str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'PROOF.md',ROOT/'verify.py',ROOT.parent/'mod8_transfer'/'PROOF.md',ROOT.parent/'covered_response'/'verify.py']}},sort_keys=True,indent=2))
