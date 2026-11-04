from pathlib import Path
from itertools import product
import hashlib,json,subprocess,sys
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(n):return hashlib.sha256((D/n).read_bytes()).hexdigest()
def main():
 logs={}
 for script in ['explore.py','certificates.py','even_check.py','physical_even.py','physical_mod8.py']:
  logs[script]=json.loads(subprocess.check_output([sys.executable,str(D/script)]))
 from explore import aug,op,solve
 from physical_even import build
 data=json.loads((D/'analysis.json').read_text());certs=json.loads((D/'linear_certificates.json').read_text());need(len(data)==2,'exact two quotient representatives');query_records=[];cubes=[]
 for mi,(q,cert) in enumerate(zip(data,certs)):
  need(len(q['queries'])==206,'all retained selected actual queries');equations=[];safe=[]
  for qi,r in enumerate(q['queries']):
   a,c=r['query'];need(0<=a<512 and 0<=c<512,'actual query masks');S,U=r['source'],r['target'];need(len(S)==len(U)==56,'complete55 run and1 cycle channels');d=[s-u for s,u in zip(S,U)];baseline=(d[0]+2*sum(d[1:10])+4*sum(d[19:55]))%8;need(baseline%4==0 and baseline==r['baseline'],'complete affine constant');need(all(x%2==0 for x in d[10:19]),'return-segment parity');bits=[x%2 for x in d[1:10]]+[(x//2)%2 for x in d[10:19]];need(bits==r['linear_bits'],'all18 affine derivatives');need(r['cycle_safe']==(d[55]%2==0),'exact cycle-safety test');eq=(sum(b<<i for i,b in enumerate(bits)),baseline//4);equations.append(eq)
   if r['cycle_safe']:safe.append(eq)
   query_records.append(dict(model_index=mi,query_index=qi,query=r['query'],source_augmented_rows=aug(q['rows'],a,c),target_augmented_rows=aug(op(q['rows']),a,c),port_names=[9,10],source_rooted_cycle_count=S[-1],target_rooted_cycle_count=U[-1]))
  need(solve(equations)==q['all_queries'] and solve(safe)==q['safe_queries'],'exact affine-system solution')
  need(len(safe)==[126,119][mi],'all safe actual queries retained')
  certificate=cert['safe_certificate'];ids=certificate.get('contradiction_query_indices',certificate['independent_query_indices']);chosen=[]
  for i in ids:
   r=q['queries'][i];need(r['cycle_safe'],'certificate actual queries all safe');chosen.append((sum(b<<j for j,b in enumerate(r['linear_bits'])),r['baseline']//4))
  zeros=[x for x in range(1<<18) if all(((row&x).bit_count()%2)==rhs for row,rhs in chosen)]
  need(zeros==([511] if mi==0 else []),'all262144 profile assignments checked against concrete certificate');cubes.append(dict(model_index=mi,assignments_checked=1<<18,zero_assignments=zeros,odd_total_zero_assignments=[x for x in zeros if (x>>9).bit_count()%2]))
 even=json.loads((D/'even_bank.json').read_text());need(even['rows']==data[0]['rows'],'same physical carousel orientation');expected_shapes={tuple(1+x for x in e) for e in product(range(4),repeat=9) if sum(e)<=3}
 for qi,r in enumerate(even['queries']):
  need(set(map(tuple,r['shapes']))==expected_shapes and len(r['shapes'])==220,'complete geometric mod16 run domain');need(r['geometric_contrast_mod16']==8*qi and r['q2_derivative_parities']==[0,0,0,0,1,0,0,1,0],'local even-profile identities')
 for e in range(512):
  gaps=[(r['geometric_contrast_mod16']+8*sum(bit*((e>>i)&1) for i,bit in enumerate(r['q2_derivative_parities'])))%16 for r in even['queries']];need(sorted(gaps)==[0,8],'every even-profile q2-bit assignment has a nonzero probe')
  for odd in range(1,16,2):need([(odd*g)%16 for g in gaps]==gaps,'normalization and ambient odd-unit safety')
 physical_count=0
 for name in ['physical_even_records.json','physical_mod8_records.json']:
  records=json.loads((D/name).read_text())
  for r in records:
   if name=='physical_even_records.json':query=even['queries'][r['query_index']]['query']
   else:query=data[r['model_index']]['queries'][r['query_index']]['query']
   expected=build(r['quotient_rows'],r['child_sizes'],*query)
   need(all(r[k]==v for k,v in expected.items()),'retained actual physical rows/masks/child sets');T=r['augmented_rows'];need(all(not(row>>i&1) and all(((row>>j)&1)+((T[j]>>i)&1)==1 for j in range(len(T)) if j!=i) for i,row in enumerate(T)),'physical tournament');need(r['H']>0 and r['H']%2==1,'Hamilton oddness on physical controls');physical_count+=1
 need(physical_count==3520,'all literal interpolation graphs')
 (D/'query_records.json').write_text(json.dumps(query_records,separators=(',',':'))+'\n');(D/'profile_cube_checks.json').write_text(json.dumps(cubes,indent=2)+'\n')
 names=['PROOF.md','explore.py','runs.cpp','certificates.py','even_check.py','even_runs.cpp','physical_h.cpp','physical_even.py','physical_mod8.py','analysis.json','run_values.txt','linear_certificates.json','even_bank.json','even_values.txt','physical_even_records.json','physical_even_values.txt','physical_mod8_records.json','physical_mod8_values.txt','query_records.json','profile_cube_checks.json']
 print(json.dumps(dict(status='PASS',scope='Two displayed nine-name prime circulant classes, arbitrary matching children and fixed ambient; not a prime9 census',actual_query_pairs=412,exact_run_and_cycle_counts=412*2*56,certifying_safe_queries=[18,6],profile_assignments_checked=2*(1<<18),complete_even_shapes=220,exact_even_run_counts=880,literal_physical_graphs=3520,physical_control_summaries={k:v for k,v in logs.items() if k.startswith('physical_')},sha256={n:sha(n) for n in names}),sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
