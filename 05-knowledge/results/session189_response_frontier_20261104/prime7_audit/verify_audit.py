#!/usr/bin/env python3
"""Independent review checks. Reads producer data, imports no producer code."""
import functools,itertools,json,hashlib,sys,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
P=ROOT.parent/'session189-profiles'
E=ROOT.parent/'session189-prime7'
if not P.is_dir():P=ROOT.parent/'profiles'
if not E.is_dir():E=ROOT.parent/'prime7'

def need(x,*msg):
 if not x:raise RuntimeError(msg)
def read(name):return json.loads((P/name).read_text())
def opposite(T):return tuple(((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T))
def valid(T):return all(not(r>>i&1) and all(((r>>j)&1)+((T[j]>>i)&1)==1 for j in range(len(T)) if i!=j) for i,r in enumerate(T))
def prime(T):
 n=len(T)
 return not any(all((T[i]&m) in (0,m) for i in range(n) if not(m>>i&1)) for m in range(1,(1<<n)-1) if m.bit_count()>1)
def starts(T,multiplicity):
 n=len(T);neighbors=tuple(tuple(j for j in range(n) if r>>j&1) for r in T)
 @functools.lru_cache(None)
 def go(rem,last):
  if not any(rem):return 1
  total=0
  for j in neighbors[last]:
   if rem[j]:
    nxt=list(rem);nxt[j]-=1;total+=go(tuple(nxt),j)
  return total
 out=[]
 for i in range(n):
  rem=list(multiplicity);rem[i]-=1;out.append(go(tuple(rem),i))
 return out

def ext(T,a,c=None):
 n=len(T);out=[r|(((a>>i)&1)<<n) for i,r in enumerate(T)]
 out.append(((1<<n)-1)^a)
 if c is not None:
  for i in range(n):out[i]|=((c>>i)&1)<<(n+1)
  out[n]|=1<<(n+1);out.append(((1<<n)-1)^c)
 return tuple(out)

def cycle_parity(T,root):
 n=len(T);full=(1<<n)-1;dp=[0]*(1<<n);dp[1<<root]=1<<root
 for mask in range(1<<n):
  ends=dp[mask]
  while ends:
   bit=ends&-ends;ends-=bit;i=bit.bit_length()-1
   options=T[i]&(full^mask)
   while options:
    q=options&-options;options-=q;dp[mask|q]^=q
 closing=sum(1<<i for i,r in enumerate(T) if r>>root&1)
 return (dp[full]&closing).bit_count()%2

def rank(vs):
 piv={}
 for v in vs:
  while v:
   j=v.bit_length()-1
   if j in piv:v^=piv[j]
   else:piv[j]=v;break
 return len(piv)

def audit_endpoints_and_anchors():
 reps=read('quotients.json');need(len(reps)==197,'representative count')
 byrows={tuple(r['rows']):r for r in reps};need(len(byrows)==197,'distinct representatives')
 for row in reps:
  T=tuple(row['rows']);need(valid(T) and prime(T),'invalid/nonprime',T)
  ss=[0]*7;ee=[0]*7
  for p in itertools.permutations(range(7)):
   if all(T[p[i]]>>p[i+1]&1 for i in range(6)):ss[p[0]]+=1;ee[p[-1]]+=1
  need(ss==row['S'] and ee==row['E'],'integer endpoint count',T)
  need([(a+b)%2 for a,b in zip(ss,ee)]==row['p'],'endpoint parity',T)
 recs=read('anchor_exploration.json');need(len(recs)==38,'flat representative count');blind=[];primitive=0
 for rec in recs:
  T=tuple(rec['rows']);U=opposite(T);row=byrows[T]
  need(not any(row['p']),'anchor representative nonflat')
  skew=[[0]*7 for _ in range(7)]
  for j in range(7):
   mult=tuple(2 if i==j else 1 for i in range(7));s=starts(T,mult);u=starts(U,mult)
   for i in range(7):skew[i][j]=s[i]-u[i]
  need(skew==rec['start_run_skews'],'anchor starts',T);primitive+=49
  bases=[];doubles=[]
  for i in range(7):
   tx=ext(T,1<<i);ux=ext(U,1<<i)
   bases.append(sum(starts(tx,(1,)*8))-sum(starts(ux,(1,)*8)))
   dd=[]
   for j in range(7):
    mult=tuple(2 if v==j else 1 for v in range(8));dd.append(sum(starts(tx,mult))-sum(starts(ux,mult)))
   doubles.append(dd)
  need(bases==rec['one_port_base'] and doubles==rec['one_port_double'],'one-port primitives',T);primitive+=56
  bad=[]
  for V,got in enumerate(rec['patterns']):
   d=[(bases[i]+2*sum(doubles[i][j] for j in range(7) if V>>j&1))%4 for i in range(7)]
   t=[(row['S'][i]-row['E'][i]+2*sum(skew[i][j] for j in range(7) if V>>j&1))%4 for i in range(7)]
   a=[(2*x+2*y)%8 for x,y in zip(t,d)]
   need(got==dict(V=V,first=d,second=a),'anchor pattern')
   if not any(d+a):bad.append(V);blind.append((T,V))
  need(bad==rec['blind'],'blind pattern list')
 residual=read('residual_mod4.json');expected={(tuple(r['Q']),r['V']) for r in residual}
 need(expected=={(T,V) for T,V in blind if V!=127},'residual singleton coverage')
 for r in residual:
  masks=[x for x in range(128) if x.bit_count()<=2];idx=masks.index(r['query'][0])*29+masks.index(r['query'][1]);need(r['contrast'][idx]==2,'residual actual witness')
 need(len(blind)==10 and sum(V==127 for T,V in blind)==4,'ten residuals')
 return dict(representatives=197,endpoint_integers=197*14,anchor_primitive_integers=primitive,anchor_patterns=38*128,blind_patterns=10,singleton_residuals=6)

def audit_safe_affine():
 cases=read('affine_banks.json');ambient=read('ambient_analysis.json');masks=[x for x in range(128) if x.bit_count()<=2];coords=list(itertools.product(masks,repeat=2));pos={v:i for i,v in enumerate(coords)}
 subpos=[[j for j,(s,t) in enumerate(coords) if s&a==s and t&c==t] for a,c in coords]
 def actual(coeff):return [sum(coeff[j] for j in js) for js in subpos]
 summaries=[]
 for case,am in zip(cases,ambient):
  T=tuple(case['Q']);U=opposite(T);need(list(T)==am['Q'],'matching quotient order')
  sides=case['full_run_coefficients'];A,B=sides
  difference={key:([a-b for a,b in zip(A[key],B[key])] if key=='K0' else [[a-b for a,b in zip(x,y)] for x,y in zip(A[key],B[key])]) for key in A}
  coeff=[difference['K0'][i]+2*sum(v[i] for v in difference['K2'])+4*sum(v[i] for v in difference['K22']) for i in range(841)]
  base=actual(coeff);cols=[actual([4*x for x in v]) for v in difference['K2']]+[actual([2*x for x in v]) for v in difference['K3']]
  def bits(vals):
   need(all(x%4==0 for x in vals),'not binary contrast');return sum(((x%8)//4)<<i for i,x in enumerate(vals))
  b=bits(base);cv=[bits(v) for v in cols]
  cycle=[cycle_parity(ext(T,a,c),7)^cycle_parity(ext(U,a,c),7) for a,c in coords]
  need(cycle==am['cycle_parity_contrast'],'independent complete cycle-parity bank')
  safe=[i for i,z in enumerate(cycle) if not z];need(safe==am['safe_indices'],'safe indices');sf=sum(1<<i for i in safe)
  local_blind=[];safe_blind=[];bad_odd=[];v=b;previous=0
  for i in range(1<<14):
   gray=i^(i>>1)
   if i:v^=cv[(gray^previous).bit_length()-1]
   if v==0:local_blind.append(gray)
   if not(v&sf):
    safe_blind.append(gray)
    if (gray>>7).bit_count()%2:bad_odd.append(gray)
   previous=gray
  local_blind.sort();safe_blind.sort()
  need(not bad_odd,'uncovered odd-total assignment')
  need(local_blind==am['local_blind_masks'] and safe_blind==am['safe_blind_masks'],'complete residual classification')
  need(rank(cv)==case['rank'] and rank([x&sf for x in cv])==am['safe_rank'],'ranks')
  summaries.append(dict(rows=list(T),safe_queries=len(safe),full_rank=rank(cv),safe_rank=rank([x&sf for x in cv]),local_blind=local_blind,safe_blind=safe_blind,odd_uncovered=[]))
 return summaries

def physical_rows(Q,ns,a,c):
 blocks=[];n=0
 for size in ns:blocks.append(list(range(n,n+size)));n+=size
 out=[0]*(n+2)
 for i,block in enumerate(blocks):
  for u in block:
   for v in block:
    if u<v:out[u]|=1<<v
   for j,other in enumerate(blocks):
    if Q[i]>>j&1:
     for v in other:out[u]|=1<<v
   if a>>i&1:out[u]|=1<<n
   else:out[n]|=1<<u
   if c>>i&1:out[u]|=1<<(n+1)
   else:out[n+1]|=1<<u
 out[n]|=1<<(n+1)
 return out

def packed_covers(T):
 # Backward endpoint recurrence, five coefficients packed in base 2^64.
 # Every intermediate coefficient is at most 13! * binomial(12,4) < 2^64.
 n=len(T);need(n<=13,'packed physical order');width=64;limit=(1<<(5*width))-1
 incoming=[sum(1<<i for i,r in enumerate(T) if r>>j&1) for j in range(n)]
 dp=[[0]*n for _ in range(1<<n)];total=[0]*(1<<n)
 for mask in range(1,1<<n):
  left=mask
  while left:
   bit=left&-left;left-=bit;j=bit.bit_length()-1;previous=mask^bit
   if not previous:value=1
   else:
    value=total[previous]<<width;allowed=incoming[j]&previous
    while allowed:
     b=allowed&-allowed;allowed-=b;value+=dp[previous][b.bit_length()-1]
    value&=limit
   dp[mask][j]=value;total[mask]+=value
  need(total[mask]>>320==0,'packed overflow')
 return [(total[-1]>>(h*width))&((1<<width)-1) for h in range(5)]

def audit_even_identities():
 rows=json.loads((E/'rows.json').read_text());qs=json.loads((E/'even_queries.json').read_text());bank=[list(map(int,s.split())) for s in (E/'even_bank.txt').read_text().splitlines()]
 need(len(qs)==841 and len(bank)==4*841*2 and all(len(v)==110 for v in bank),'even complete bank dimensions')
 specs=[(0,[1,0],8,[]),(1,[0,3],0,[3,5]),(1,[4,68],8,[3,5]),(2,[0,3],0,[0,1]),(2,[4,68],8,[0,1]),(3,[3,1],8,[])]
 for ri,q,b,active in specs:
  idx=2*(ri*841+qs.index(q));s,t=bank[idx:idx+2];actual=[(x-y)%(16 if j<5 else 2) for j,(x,y) in enumerate(zip(s,t))]
  expected=[b,0,0,0,0]+[int(k==2 and i in active and h==1) for k in (2,3,5) for i in range(7) for h in range(1,6)]
  need(actual==expected,'even full-profile identity',ri,q)
 ind=json.loads((E/'independent.json').read_text());models=[tuple(v) for v in ind['models']]
 expected=[ns for ns in itertools.product(range(1,5),repeat=7) if sum(ns)-7<=3]+[tuple(5 if i==j else 1 for i in range(7)) for j in range(7)]
 need(set(models)==set(expected) and len(models)==127,'independent model completeness')
 # Recompute the triangular evaluation, not the producer inversion, from each
 # retained nonnegative kernel and the exact Stirling profile.
 stir=[[0]*6 for _ in range(6)];stir[0][0]=1
 import math
 for n in range(1,6):
  for k in range(1,6):stir[n][k]=k*stir[n-1][k]+stir[n-1][k-1]
 profiles=[[math.factorial(k)*stir[n][k] for k in range(6)] for n in range(6)]
 for rec in ind['records']:
  kernels={tuple(v[:7]):v[7:] for v in rec['interpolated']};need(set(kernels)==set(models),'complete independent kernel')
  for ns,values in zip(models,rec['values']):
   rebuilt=[0]*5
   for ks,kv in kernels.items():
    if all(k<=n for k,n in zip(ks,ns)):
     weight=math.prod(profiles[n][k] for k,n in zip(ks,ns))
     for h in range(5):rebuilt[h]+=weight*kv[h]
   need(rebuilt==values,'full integer physical re-evaluation')
  selected=[(1,)*7]+[tuple(5 if i==j else 1 for i in range(7)) for j in (0,3,6)]+[(1,3,1,1,1,2,1),(2,1,2,1,2,1,1)]
 controls=[]
 for rec in ind['records']:
  Q=tuple(rows[rec['quotient_index']]);Q=Q if rec['side']==0 else opposite(Q)
  for ns in selected:
   T=physical_rows(Q,ns,*rec['query']);values=packed_covers(T);expected_values=rec['values'][models.index(ns)]
   need(values==expected_values,'independent packed physical cover DP')
   controls.append(dict(quotient_index=rec['quotient_index'],query=rec['query'],side=rec['side'],sizes=list(ns),physical_rows=T,cover_counts=values))
 (ROOT/'physical_controls.json').write_text(json.dumps(controls,sort_keys=True,indent=2)+'\n')
 return dict(queries=6,affine_bits=21,cover_coordinates=5,physical_interpolation_cases=1524,shapes=127,additional_independent_physical_graphs=len(controls))

def audit_census():
 reps=read('quotients.json')
 with tempfile.TemporaryDirectory() as tmp:
  exe=str(Path(tmp)/'owner_check')
  subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'owner_check.cpp'),'-o',exe],check=True)
  data=str(len(reps))+'\n'+''.join(' '.join(map(str,r['rows']))+'\n' for r in reps)
  output=subprocess.check_output([exe,str(P/'prime_owners.u16le')],input=data.encode()).decode()
 result=json.loads(output)
 need(result['prime_labels']==973680,'total labeled primes')
 return result

def main():
 result={'status':'PASS','census':audit_census(),'endpoints_and_anchors':audit_endpoints_and_anchors(),'safe_affine':audit_safe_affine(),'even_identities':audit_even_identities()}
 filenames=['quotients.json','anchor_exploration.json','affine_banks.json','ambient_analysis.json','residual_mod4.json','prime_owners.u16le']
 result['input_sha256']={'profiles/'+n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in filenames}
 for n in ['rows.json','even_queries.json','even_bank.txt','independent.json']:
  result['input_sha256']['prime7/'+n]=hashlib.sha256((E/n).read_bytes()).hexdigest()
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
