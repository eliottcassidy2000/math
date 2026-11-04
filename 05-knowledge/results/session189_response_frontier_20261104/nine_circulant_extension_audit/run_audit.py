"""Independent suffix multiset words, rooted cycles, affine certificates and full type domain."""
from pathlib import Path
from itertools import combinations,product
import json,hashlib,subprocess,tempfile
D=Path(__file__).resolve().parent;P=D.parent/'nine_circulant_extension'
if not P.exists():P=D.parent/'session189-nine-circulant-extension'
def need(c,s):
 if not c:raise RuntimeError(s)
def op(R):return [((1<<len(R))-1)^(1<<i)^r for i,r in enumerate(R)]
def aug(R,a,c):
 n=len(R);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(R)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]
def gf_rank(rows):
 piv={}
 for r in rows:
  while r:
   j=r.bit_length()-1
   if j not in piv:piv[j]=r;break
   r^=piv[j]
 return len(piv)
def main():
 data=json.loads((P/'analysis.json').read_text());cert=json.loads((P/'linear_certificates.json').read_text());even=json.loads((P/'even_bank.json').read_text());need([d['rows'] for d in data]==[[480,449,387,263,15,30,60,120,240],[226,452,393,275,39,78,156,312,113]],'two explicitly bounded orientations')
 shapes=[[1]*9]
 for k in [2,3]:
  for i in range(9):s=[1]*9;s[i]=k;shapes.append(s)
 for i,j in combinations(range(9),2):s=[1]*9;s[i]=s[j]=2;shapes.append(s)
 cases=[];case_meta=[]
 for mi,q in enumerate(data):
  need(len(q['queries'])==206 and len({tuple(r['query']) for r in q['queries']})==206,'206 distinct actual query pairs')
  for qi,r in enumerate(q['queries']):
   for side in [0,1]:cases.append((aug(op(q['rows']) if side else q['rows'],*r['query']),shapes));case_meta.append((mi,qi,side))
 even_start=len(cases)
 expected=set(tuple(1+x for x in v) for v in product(range(4),repeat=9) if sum(v)<=3)
 for qi,r in enumerate(even['queries']):
  need(len(r['shapes'])==220 and set(map(tuple,r['shapes']))==expected,'all220 geometric excess shapes')
  for side in [0,1]:cases.append((aug(op(even['rows']) if side else even['rows'],*r['query']),r['shapes']))
 with tempfile.TemporaryDirectory() as tmp:
  binary=Path(tmp)/'runs';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'run_suffix.cpp'),'-o',str(binary)],check=True);inp=str(len(cases))+'\n'+''.join(str(len(R))+' '+str(len(ss))+'\n'+' '.join(map(str,R))+'\n'+''.join(' '.join(map(str,list(s)+[1,1]))+'\n' for s in ss) for R,ss in cases);raw=subprocess.check_output([str(binary)],input=inp.encode()).decode()
 values=[list(map(int,l.split())) for l in raw.splitlines()];need(len(values)==len(cases),'all suffix counters returned');recount=[];equations=[[],[]];safe=[[],[]]
 for (mi,qi,side),v in zip(case_meta,values):
  r=data[mi]['queries'][qi];need(v==r['target' if side else 'source'],'every independent exact55 run and rooted-cycle coordinate');need(all(x%2==0 for x in v[10:19]),'individual triple-return coefficients even');recount.append(dict(model=mi,query_index=qi,query=r['query'],side=side,augmented_rows=cases[len(recount)][0],counts=v))
 for mi,q in enumerate(data):
  for qi,r in enumerate(q['queries']):
   d=[s-u for s,u in zip(r['source'],r['target'])];baseline=(d[0]+2*sum(d[1:10])+4*sum(d[19:55]))%8;need(baseline%4==0 and all(x%2==0 for x in d[10:19]),'divisible affine constant and triple return terms');row=sum((d[i+1]%2)<<i for i in range(9))+sum(((d[i+10]//2)%2)<<(i+9) for i in range(9));rhs=baseline//4;need([row>>j&1 for j in range(18)]==r['linear_bits'] and baseline==r['baseline'],'every affine bit from independently recounted actual values');is_safe=d[-1]%2==0;need(is_safe==r['cycle_safe'],'source/target exact cycle-parity equality');equations[mi].append((row,rhs))
   if is_safe:safe[mi].append(qi)
 need([len(s) for s in safe]==[126,119],'safe actual-query totals');need(all(((row&511).bit_count()&1)==rhs for row,rhs in equations[0]),'E9 actually zeros every retained C0 equation');cert_records=[]
 for mi,c in enumerate(cert):
  C=c['safe_certificate'];ids=C.get('contradiction_query_indices',C['independent_query_indices']);need(len(ids)==[18,6][mi] and len(set(ids))==len(ids) and all(i in safe[mi] for i in ids),'all selected actual certificate equations are distinct and safe');rows=[equations[mi][i] for i in ids]
  if mi==0:need(gf_rank([r for r,b in rows])==18 and all(((r&511).bit_count()&1)==b for r,b in rows),'eighteen independent equations uniquely select E9')
  else:
   lhs=rhs=0
   for r,b in rows:lhs^=r;rhs^=b
   need(lhs==0 and rhs==1,'six explicit safe equations xor to contradiction')
  zeros=[]
  for x in range(1<<18):
   if all(((row&x).bit_count()&1)==rhs for row,rhs in rows):zeros.append(x)
  need(zeros==([511] if mi==0 else []),'complete262144-assignment upper profile domain');need(not any((x>>9).bit_count()%2 for x in zeros),'no uncovered odd-total profile assignment');cert_records.append(dict(model=mi,query_indices=ids,actual_queries=[data[mi]['queries'][i]['query'] for i in ids],equations=rows,rank=gf_rank([r for r,b in rows]),complete_zero_assignments=zeros))
 even_records=[]
 for qi,r in enumerate(even['queries']):
  vals=values[even_start+2*qi:even_start+2*qi+2];need(vals[0][:-1]==r['source_counts'] and vals[1][:-1]==r['target_counts'],'all880 independently recounted geometric run coordinates');delta=[a-b for a,b in zip(vals[0][:-1],vals[1][:-1])];constant=sum((-2)**(sum(s)-9)*v for s,v in zip(r['shapes'],delta))%16;derivatives=[]
  for i in range(9):shape=[1]*9;shape[i]=2;derivatives.append(delta[r['shapes'].index(shape)]%2)
  need(constant==8*qi and derivatives==[0,0,0,0,1,0,0,1,0],'two E9 mod16 equations with exact physical derivative positions');even_records.append(dict(query=r['query'],shapes=r['shapes'],source_counts=vals[0][:-1],target_counts=vals[1][:-1],geometric_contrast_mod16=constant,q2_derivatives=derivatives))
 for e in range(512):
  gaps=[(r['geometric_contrast_mod16']+8*sum(((e>>i)&1)*v for i,v in enumerate(r['q2_derivatives'])))%16 for r in even_records];need(sorted(gaps)==[0,8],'every q2-bit assignment has exactly one gap8');need(all(sorted((u*g)%16 for g in gaps)==[0,8] for u in range(1,16,2)),'common Hamilton and ambient odd-unit safety')
 outputs={'run_recount.json':recount,'certificate_recount.json':cert_records,'even_recount.json':even_records}
 for name,obj in outputs.items():
  f=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if f.exists():need(f.read_text()==txt,'independent exact artifact stable')
  else:f.write_text(txt)
 return dict(status='PASS',actual_query_pairs=412,exact_run_and_cycle_coordinates=824*56,selected_safe_equations=[18,6],profile_assignments_checked=2*(1<<18),even_shapes=220,exact_even_run_coordinates=880,even_q2_assignments_checked=512,source_pins={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['run_audit.py','run_suffix.cpp']},output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs},producer_proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
