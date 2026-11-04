"""Independent four-row audit. Reads frozen data, never imports producer code."""
from pathlib import Path
from itertools import permutations,combinations
from array import array
import hashlib,json,subprocess,sys
D=Path(__file__).resolve().parent
P=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else next((d for d in [D.parent/'session189-nine-noncirculant-probes',D.parent/'nine_noncirculant_probes',D.parent/'nine_noncirculant'] if (d/'PROOF.md').is_file()),None)
if P is None:raise RuntimeError('Frozen producer sibling is unavailable')
def need(ok,why):
 if not ok:raise RuntimeError(why)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(name):return json.loads((P/name).read_text())
def dump(name,data):(D/name).write_text(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n')
Qs=[[192,289,99,71,15,153,306,94,189],[320,193,83,71,43,15,432,61,190],[448,209,195,71,45,15,432,312,62],[384,177,147,135,105,77,15,368,126]]
def converse(rows):
 n=len(rows);return [sum(1<<j for j in range(n) if j!=i and (rows[j]>>i)&1) for i in range(n)]
def extend(rows,a,c):
 n=len(rows);r=rows+[0,0]
 for i in range(n):
  for t,mask in ((n,a),(n+1,c)):
   if (mask>>i)&1:r[i]|=1<<t
   else:r[t]|=1<<i
 r[n]|=1<<(n+1);return r

def valid(rows):
 n=len(rows);return all(0<=r<(1<<n) and not (r>>i)&1 and all(((r>>j)&1)!=((rows[j]>>i)&1) for j in range(i+1,n)) for i,r in enumerate(rows))
def prime(rows):
 n=len(rows)
 return not any(all((rows[v]&s) in (0,s) for v in range(n) if not (s>>v)&1) for s in range(1,(1<<n)-1) if s.bit_count()>1)
def physical(Q,sizes,query,ambient=False):
 labels=[];owner=[]
 for i,k in enumerate(sizes):labels.append(list(range(len(owner),len(owner)+k)));owner.extend([i]*k)
 n=len(owner);out=[0]*n
 for u in range(n):
  for v in range(n):
   if (owner[u]==owner[v] and u<v) or (owner[u]!=owner[v] and (Q[owner[u]]>>owner[v])&1):out[u]|=1<<v
 a,c=(sum(1<<u for u,i in enumerate(owner) if (mask>>i)&1) for mask in query)
 exterior=[]
 if ambient:
  exterior=[n,n+1]
  out=[r|(1<<n) for r in out]+[1<<(n+1),(1<<n)-1]
  a|=1<<(n+1);c|=1<<(n+1)
 return dict(quotient_rows=Q,child_sizes=sizes,child_sets=labels,augmented_rows=extend(out,a,c),physical_query=[a,c],port_names=[len(out),len(out)+1],exterior_names=exterior)
def engine(jobs,filename):
 strings=[]
 for mode,rows,extra in jobs:
  strings.append(' '.join(map(str,[mode,len(rows),*rows,*extra])))
 raw=subprocess.check_output([str(D/'recount')],input=(str(len(jobs))+'\n'+'\n'.join(strings)+'\n').encode())
 (D/filename).write_bytes(raw)
 result=[list(map(int,line.split())) for line in raw.splitlines()]
 need(len(result)==len(jobs),'recount output cardinality')
 return result

def rank(rows):
 work=list(rows);r=0
 for col in range(18):
  pivot=next((i for i in range(r,len(work)) if (work[i]>>col)&1),None)
  if pivot is None:continue
  work[r],work[pivot]=work[pivot],work[r]
  for j in range(len(work)):
   if j!=r and (work[j]>>col)&1:work[j]^=work[r]
  r+=1
 return r

def main():
 pins={'PROOF.md':'0aae8c5a8a93eb5dd9bf5094e969d14417ef159badb5eedd55da028e49e648e2','verify.py':'8ea6433a3625de25cf67a8f82db9540a2695922e3217a06024f167d1661b153e','normal.json':'f5d1d7c19093aa42195a0f019409787c23437d409a7d267b586566fbb14a0d97','optimized.json':'f5d1d7c19093aa42195a0f019409787c23437d409a7d267b586566fbb14a0d97'}
 for name,h in pins.items():need(digest(P/name)==h,'frozen '+name)
 for name,h in read('normal.json')['sha256'].items():need(digest(P/name)==h,'producer receipt pin '+name)
 dependencies=read('DEPENDENCIES.json')
 for r in dependencies:
  p=D.parent/r['path']
  if not p.is_file():p=D.parent/r['path'].replace('session189-','').replace('-audit/','_audit/')
  need(digest(p)==r['sha256'],'dependency pin '+r['path'])
 source_hashes={name:digest(P/name) for name in read('PUBLICATION_FILES.json')}
 a=read('anchor_data.json');d=read('analysis.json');cert=read('linear_certificates.json')
 need(len(a)==len(d)==len(cert)==4,'four models')
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'recount.cpp'),'-o',str(D/'recount')],check=True)
 endpoints=[];anchorjobs=[]
 for mi,Q in enumerate(Qs):
  need(a[mi]['rows']==d[mi]['rows']==cert[mi]['rows']==Q,'fixed physical rows')
  need(valid(Q) and prime(Q),'prime valid fixed quotient')
  joint=[[0]*9 for _ in range(9)]
  for perm in permutations(range(9)):
   if all((Q[u]>>v)&1 for u,v in zip(perm,perm[1:])):joint[perm[0]][perm[-1]]+=1
  need(joint==a[mi]['joint_endpoints'],'full joint endpoint recount')
  S=list(map(sum,joint));E=[sum(row[i] for row in joint) for i in range(9)]
  need(S==a[mi]['S'] and E==a[mi]['E'],'endpoint marginals')
  need(sum(S)==[2251,2239,3027,2479][mi],'Hamilton count')
  need(a[mi]['p']==[(x+y)%2 for x,y in zip(S,E)]==[0]*9,'endpoint parity')
  need(len(set(x.bit_count() for x in Q))>1,'not regular')
  need(a[mi]['nontrivial_modules']==[],'module report')
  endpoints.append(dict(rows=Q,joint=joint,starts=S,ends=E))
  for j in range(9):
   cap=[1]*9;cap[j]=2
   for side in range(2):anchorjobs.append(('A',converse(Q) if side else Q,cap))
  for i in range(9):
   for j in range(-1,9):
    cap=[1]*10
    if j>=0:cap[j]=2
    for side in range(2):
     row=converse(Q) if side else Q;one=row+[0]
     for v in range(9):
      if v==i:one[v]|=512
      else:one[9]|=1<<v
     anchorjobs.append(('R',one,cap))
 print('Endpoint recount passed; recounting anchors.',file=sys.stderr,flush=True)
 ar=iter(engine(anchorjobs,'anchor_exact.txt'));anchor_recount=[]
 for mi,Q in enumerate(Qs):
  start=[[0]*9 for _ in range(9)]
  for j in range(9):
   s,t=next(ar),next(ar)
   for i in range(9):start[i][j]=s[i]-t[i]
  d0=[];d2=[]
  for i in range(9):
   base=next(ar)[0]-next(ar)[0];d0.append(base);row=[]
   for j in range(9):row.append(next(ar)[0]-next(ar)[0])
   d2.append(row)
  need(start==a[mi]['start'] and d0==a[mi]['d0'] and d2==a[mi]['d2'],'all 171 anchor primitives')
  blind=[]
  need(len(a[mi]['patterns'])==512,'512 singleton patterns')
  for V in range(512):
   first=[(d0[i]+2*sum(d2[i][j] for j in range(9) if (V>>j)&1))%4 for i in range(9)]
   t=[(a[mi]['S'][i]-a[mi]['E'][i]+2*sum(start[i][j] for j in range(9) if (V>>j)&1))%4 for i in range(9)]
   second=[2*(t[i]+first[i])%8 for i in range(9)]
   need(a[mi]['patterns'][V]==dict(V=V,first=first,second=second),'complete singleton position pattern')
   if not any(first+second):blind.append(V)
  need(blind==a[mi]['blind']==[[254,511],[511],[511],[54,511]][mi],'exact blind patterns')
  sigma=[(a[mi]['S'][i]-a[mi]['E'][i]+2*sum(start[i]))%4 for i in range(9)]
  need(sigma==a[mi]['sigma']==[0]*9,'second endpoint residues')
  anchor_recount.append(dict(rows=Q,start=start,d0=d0,d2=d2,blind=blind))
 dump('joint_endpoints.json',endpoints);dump('anchor_recount.json',anchor_recount)
 caps=[[1]*11]
 for k in (2,3):
  for i in range(9):v=[1]*11;v[i]=k;caps.append(v)
 for i,j in combinations(range(9),2):v=[1]*11;v[i]=v[j]=2;caps.append(v)
 runjobs=[]
 for mi,Q in enumerate(Qs):
  need(len(d[mi]['queries'])==206 and len(set(tuple(r['query']) for r in d[mi]['queries']))==206,'206 distinct queries')
  for q in d[mi]['queries']:
   for side in range(2):
    graph=extend(converse(Q) if side else Q,*q['query'])
    runjobs.extend(('R',graph,v) for v in caps);runjobs.append(('Z',graph,[9]))
 print('All anchors passed; recounting all 92,288 run and cycle coordinates.',file=sys.stderr,flush=True)
 raw=engine(runjobs,'run_recount.txt');run_values=[r[0] for r in raw];pos=0;affines=[];query_meta=read('query_records.json');qpos=0
 need(len(query_meta)==824,'all physical query metadata')
 for mi,Q in enumerate(Qs):
  equations=[];safe=[]
  for qi,q in enumerate(d[mi]['queries']):
   exact=[run_values[pos:pos+56],run_values[pos+56:pos+112]];pos+=112
   need(exact==[q['source'],q['target']],'exact run/cycle recount '+str((mi,qi)))
   need(all(v%2==0 for side in exact for v in side[10:19]),'each side triple-run coefficient is even')
   delta=[x-y for x,y in zip(*exact)];baseline=(delta[0]+2*sum(delta[1:10])+4*sum(delta[19:55]))%8
   need(baseline in (0,4) and baseline==q['baseline'],'affine baseline')
   bits=[v%2 for v in delta[1:10]]+[(v//2)%2 for v in delta[10:19]]
   need(bits==q['linear_bits'],'all 18 derivatives');is_safe=delta[55]%2==0
   need(is_safe==q['cycle_safe'],'actual cycle safety')
   row=sum(x<<j for j,x in enumerate(bits));eq=(row,baseline//4);equations.append(eq)
   if is_safe:safe.append(eq)
   expected=dict(model_index=mi,query_index=qi,query=q['query'],source_augmented_rows=extend(Q,*q['query']),target_augmented_rows=extend(converse(Q),*q['query']),port_names=[9,10],source_rooted_cycle_count=exact[0][55],target_rooted_cycle_count=exact[1][55])
   need(query_meta[qpos]==expected,'physical source/target query metadata');qpos+=1
  for label,e in [('all_queries',equations),('safe_queries',safe)]:
   r=rank([x for x,y in e]);need(d[mi][label]==dict(rank=r,inconsistent=True,solutions=[]),'retained affine rank/solution statement')
  need(len(safe)==d[mi]['safe_query_count'],'safe query count')
  c=cert[mi]['safe_certificate'];ids=c['contradiction_query_indices']
  need(len(ids)==[7,5,1,4][mi] and len(set(ids))==len(ids),'certificate cardinality')
  need(c['contradiction_actual_queries']==[d[mi]['queries'][i]['query'] for i in ids],'certificate physical query labels')
  inds=c['independent_query_indices'];need(c['rank']==rank([equations[i][0] for i in inds])==18,'listed independent rows')
  need(c['independent_actual_queries']==[d[mi]['queries'][i]['query'] for i in inds] and all(d[mi]['queries'][i]['cycle_safe'] for i in inds),'independent actual safe queries')
  row=rhs=0
  for i in ids:
   need(d[mi]['queries'][i]['cycle_safe'],'certifying actual query safety');row^=equations[i][0];rhs^=equations[i][1]
  need(row==0 and rhs==1,'explicit XOR contradiction')
  zeros=0
  for x in range(1<<18):
   if all(((equations[i][0]&x).bit_count()%2)==equations[i][1] for i in ids):zeros+=1
  need(zeros==0,'complete 18-bit domain')
  if mi==2:
   q=d[mi]['queries'][ids[0]];need(q['query']==[3,132] and q['baseline']==4 and not any(q['linear_bits']) and q['source'][55]%2==q['target'][55]%2==0,'universal Q2 safe query')
  affines.append(dict(model_index=mi,safe_count=len(safe),rank=rank([x for x,y in safe]),certificate_indices=ids,certificate_queries=c['contradiction_actual_queries'],certificate_equations=[equations[i] for i in ids],assignments_checked=1<<18))
 dump('affine_recount.json',affines)
 print('All run/cycle and profile certificates passed; recounting named and literal physical banks.',file=sys.stderr,flush=True)
 banks=read('named_banks.json');phys=read('physical_mod8_records.json');residual=read('residual_records.json');masks=[m for m in range(512) if m.bit_count()<=2];physicaljobs=[];physical_keys=[]
 need(len(banks)==8 and len(phys)==1870 and len(residual)==4,'physical bank cardinalities')
 for bi,b in enumerate(banks):
  Q=converse(Qs[bi//2]) if bi%2 else Qs[bi//2]
  need(b['rows']==Q and b['side']==bi%2 and b['physical_masks']==masks and b['port_names']==[9,10] and b['port_edge']=='9->10' and b['a_convention']=='old vertex points to9' and b['c_convention']=='old vertex points to10','complete named physical labels')
  need(b['name']==['node860','node2005','node4067','node6501'][bi//2],'fixed named label')
  need(len(b['basis_values'])==len(b['coefficients'])==46 and all(len(r)==46 for r in b['basis_values']+b['coefficients']),'complete named arrays')
  for a0 in masks:
   for c0 in masks:physicaljobs.append(('H',extend(Q,a0,c0),[]));physical_keys.append(['named',bi,a0,c0])
 expected_keys=[]
 for mi in range(4):
  for qi in cert[mi]['safe_certificate']['contradiction_query_indices']:
   for side in range(2):
    for si in range(55):expected_keys.append((mi,qi,side,si))
 for r,key in zip(phys,expected_keys):
  mi,qi,side,si=key
  need(tuple(r[k] for k in ('model_index','query_index','side','shape_index'))==key,'every required physical interpolation case exactly once')
  Q=converse(Qs[mi]) if side else Qs[mi];expected=physical(Q,caps[si][:9],d[mi]['queries'][qi]['query']);expected.pop('exterior_names')
  need(all(r[k]==v for k,v in expected.items()),'all interpolation physical metadata')
  need(valid(r['augmented_rows']),'literal interpolation tournament')
  physicaljobs.append(('H',r['augmented_rows'],[]));physical_keys.append(['interpolation',*key])
 for ri,r in enumerate(residual):
  mi,V=[(0,254),(0,254),(3,54),(3,54)][ri]
  need(r['model_index']==mi and r['nontrivial_child_mask']==V and r['block_query']==[1,2] and r['ambient']==bool(ri%2) and r['quotient_rows']==Qs[mi] and len(r['sides'])==2,'exact singleton residual case')
  q=next(x for x in d[mi]['queries'] if x['query']==[1,2]);delta=[s-t for s,t in zip(q['source'],q['target'])]
  need((delta[0]+2*sum(delta[i+1] for i in range(9) if (V>>i)&1))%4==2,'unbounded singleton residual formula')
  for side,s in enumerate(r['sides']):
   Q=converse(Qs[mi]) if side else Qs[mi];expected=physical(Q,[1+((V>>i)&1) for i in range(9)],[1,2],bool(ri%2))
   need(all(s[k]==v for k,v in expected.items()),'complete residual and cyclic ambient physical metadata')
   need(valid(s['augmented_rows']),'actual residual tournament')
   physicaljobs.append(('H',s['augmented_rows'],[]));physical_keys.append(['residual',ri,side])
 vals=[r[0] for r in engine(physicaljobs,'physical_recount.txt')];at=0;full_hashes=[];degrees={str(i):0 for i in range(5)}
 for bi,b in enumerate(banks):
  expected=[vals[at+i*46:at+(i+1)*46] for i in range(46)];at+=2116
  need(expected==b['basis_values'],'every exact named basis value')
  coeff={}
  # Triangular inversion by increasing total subset size, distinct from direct inclusion/exclusion.
  keys=sorted(((a0,c0) for a0 in masks for c0 in masks),key=lambda k:sum(x.bit_count() for x in k))
  index={m:i for i,m in enumerate(masks)}
  for a0,c0 in keys:
   val=expected[index[a0]][index[c0]]
   for (u,v),z in coeff.items():
    if u&a0==u and v&c0==v:val-=z
   coeff[a0,c0]=val;need(val==b['coefficients'][index[a0]][index[c0]],'every named coefficient including zeros')
   degrees[str(a0.bit_count()+c0.bit_count())]+=1
  cube=[0]*(1<<18)
  for (a0,c0),v in coeff.items():cube[a0+512*c0]=v
  for bit in range(18):
   step=1<<bit
   for k in range(1<<18):
    if k&step:cube[k]+=cube[k^step]
  need(all(0<v<=39916800 and v%2==1 for v in cube),'full cube parity and bounds')
  for q in d[bi//2]['queries']:
   aa,cc=q['query'];need(cube[aa+512*cc]==q['target' if bi%2 else 'source'][0],'all exact run/bare cube overlaps including nonbasis query')
  if bi%2:
   for k in range(1<<18):
    aa,cc=k%512,k//512;need(previous_cube[k]==cube[(511-cc)+512*(511-aa)],'full named converse-complement duality')
  else:previous_cube=cube
  arr=array('I',cube)
  if sys.byteorder!='little':arr.byteswap()
  h=hashlib.sha256(arr.tobytes()).hexdigest();full_hashes.append(dict(name=b['name'],side=b['side'],values=1<<18,little_endian_u32_sha256=h))
 need(full_hashes==read('cube_hashes.json'),'all complete cube hashes')
 physicalvals=vals[at:at+1870];at+=1870
 need(physicalvals==[r['H'] for r in phys],'every physical interpolation Hamilton count')
 for first in range(0,1870,55):
  F=physicalvals[first:first+55];out=[F[0]]
  for i in range(9):need((F[1+i]-F[0])%2==0,'double interpolation integrality');out.append((F[1+i]-F[0])//2)
  for i in range(9):
   num=F[10+i]-3*F[1+i]+2*F[0];need(num%12==0,'triple interpolation and evenness');out.append(num//6)
  for k,(i,j) in enumerate(combinations(range(9),2)):
   num=F[19+k]-F[1+i]-F[1+j]+F[0];need(num%4==0,'double-double interpolation integrality');out.append(num//4)
  r=phys[first];q=d[r['model_index']]['queries'][r['query_index']]
  need(out==q['target' if r['side'] else 'source'][:55],'independent physical recovery of all certifying run channels')
 for ri,r in enumerate(residual):
  s,t=vals[at:at+2];at+=2
  need([s,t]==[v['H'] for v in r['sides']] and s-t==r['contrast'] and (s-t)%4==2,'literal residual/ambient counts and gaps')
 need(at==len(vals),'exhaustive physical recount consumption')
 dump('physical_counts.json',[dict(case=k,H=v) for k,v in zip(physical_keys,vals)])
 dump('cube_recount.json',full_hashes)
 need(degrees=={'0':8,'1':144,'2':1224,'3':5184,'4':10368},'full degree retention')
 need(source_hashes=={name:digest(P/name) for name in source_hashes},'producer unchanged')
 outputs=['joint_endpoints.json','anchor_exact.txt','anchor_recount.json','run_recount.txt','affine_recount.json','physical_recount.txt','physical_counts.json','cube_recount.json']
 receipt=dict(status='PASS',scope='Only four fixed prime-nine rows; arbitrary childwise matching profiles; every fixed ambient',models=4,joint_endpoint_coordinates=324,anchor_primitives=684,singleton_patterns=2048,run_and_cycle_coordinates=len(run_values),actual_query_metadata_records=824,certifying_battery_sizes=[7,5,1,4],profile_assignments=4*(1<<18),literal_named_basis_values=16928,named_coefficient_coordinates=16928,named_degree_counts=degrees,full_named_cube_values=8*(1<<18),named_query_crosschecks=1648,converse_complement_cube_pairs=4*(1<<18),physical_interpolation_graphs=1870,residual_physical_sides=8,maximum_literal_order=max(len(r) for _,r,_ in physicaljobs),source_sha256=source_hashes,dependency_sha256={r['path']:r['sha256'] for r in dependencies},audit_source_sha256={n:digest(D/n) for n in ['audit.py','recount.cpp']},audit_data_sha256={n:digest(D/n) for n in outputs})
 print(json.dumps(receipt,sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
