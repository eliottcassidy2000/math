"""Independent physical completion counts and direct signed-Stirling tensor inversion."""
from pathlib import Path
from itertools import product,combinations
import json,math,hashlib,subprocess,tempfile
D=Path(__file__).resolve().parent;P=D.parent/'nine_circulant_extension'
if not P.exists():P=D.parent/'session189-nine-circulant-extension'
def need(c,s):
 if not c:raise RuntimeError(s)
def op(R):return [((1<<len(R))-1)^(1<<i)^r for i,r in enumerate(R)]
def build(Q,sizes,a,c):
 starts=[];n=0
 for s in sizes:starts.append(n);n+=s
 sets=[list(range(o,o+s)) for o,s in zip(starts,sizes)];R=[0]*(n+2)
 for i,S in enumerate(sets):
  for u in S:
   R[u]=sum(1<<v for v in S if u<v)|sum(1<<v for j,A in enumerate(sets) if Q[i]>>j&1 for v in A)
   if a>>i&1:R[u]|=1<<n
   else:R[n]|=1<<u
   if c>>i&1:R[u]|=1<<(n+1)
   else:R[n+1]|=1<<u
 R[n]|=1<<(n+1);aa=sum(1<<v for i,S in enumerate(sets) if a>>i&1 for v in S);cc=sum(1<<v for i,S in enumerate(sets) if c>>i&1 for v in S)
 return dict(quotient_rows=Q,child_sizes=sizes,child_sets=sets,augmented_rows=R,physical_query=[aa,cc],port_names=[n,n+1])
def main():
 data=json.loads((P/'analysis.json').read_text());even=json.loads((P/'even_bank.json').read_text());cert=json.loads((P/'linear_certificates.json').read_text());standard=[[1]*9]
 for k in [2,3]:
  for i in range(9):s=[1]*9;s[i]=k;standard.append(s)
 for i,j in combinations(range(9),2):s=[1]*9;s[i]=s[j]=2;standard.append(s)
 groups={};allrecords=[];expected_groups=set()
 for mi,c in enumerate(cert):
  C=c['safe_certificate'];ids=C.get('contradiction_query_indices',C['independent_query_indices'])
  for qi in ids:
   for side in [0,1]:expected_groups.add(('mod8',mi,qi,side))
 for qi in range(2):
  for side in [0,1]:expected_groups.add(('even',0,qi,side))
 for kind,name in [('mod8','physical_mod8_records.json'),('even','physical_even_records.json')]:
  records=json.loads((P/name).read_text())
  for r in records:
   mi=r.get('model_index',0);qi=r['query_index'];side=r['side'];key=(kind,mi,qi,side);need(key in expected_groups,'every interpolation group is a needed explicit certificate query');query=data[mi]['queries'][qi]['query'] if kind=='mod8' else even['queries'][qi]['query'];shapes=standard if kind=='mod8' else even['queries'][qi]['shapes'];sizes=shapes[r['shape_index']];Q=data[mi]['rows'] if not side else op(data[mi]['rows']);expected=build(Q,sizes,*query);need(all(r[k]==v for k,v in expected.items()),'all exact quotient rows, child sizes/sets, augmented arcs, physical masks and port IDs');need(tuple(sizes) not in groups.setdefault(key,{}),'unique interpolation shape per group');groups[key][tuple(sizes)]=len(allrecords);allrecords.append(dict(kind=kind,model_index=mi,query_index=qi,side=side,shape_index=r['shape_index'],**expected,producer_H=r['H']))
 need(set(groups)==expected_groups and len(allrecords)==3520,'all52 complete groups and3520 physical graphs')
 for key,points in groups.items():
  shapes=standard if key[0]=='mod8' else even['queries'][key[2]]['shapes'];need(set(points)==set(map(tuple,shapes)),'complete downward-closed shape domain')
 with tempfile.TemporaryDirectory() as tmp:
  binary=Path(tmp)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'literal.cpp'),'-o',str(binary)],check=True);inp=str(len(allrecords))+'\n'+''.join(str(len(r['augmented_rows']))+'\n'+' '.join(map(str,r['augmented_rows']))+'\n' for r in allrecords);raw=subprocess.check_output([str(binary)],input=inp.encode()).decode()
 vals=list(map(int,raw.split()));need(len(vals)==3520,'every independent physical Hamilton count returned')
 for r,v in zip(allrecords,vals):need(v==r.pop('producer_H'),'exact independent physical Hamilton count');r['H']=v
 # Signed Stirling numbers of the first kind invert the entire transform directly.
 first=[[0]*5 for _ in range(5)];second=[[0]*5 for _ in range(5)];first[0][0]=second[0][0]=1
 for n in range(1,5):
  for k in range(1,n+1):first[n][k]=first[n-1][k-1]-(n-1)*first[n-1][k];second[n][k]=second[n-1][k-1]+k*second[n-1][k]
 recovered=[];channels={}
 for key,points in sorted(groups.items()):
  kind,mi,qi,side=key;shapes=standard if kind=='mod8' else even['queries'][qi]['shapes'];K={}
  for shape in shapes:
   num=0
   for ns in product(*(range(1,k+1) for k in shape)):
    weight=math.prod(first[k][n] for k,n in zip(shape,ns));num+=weight*allrecords[points[ns]]['H']
   denominator=math.prod(math.factorial(k) for k in shape);need(num%denominator==0,'exact direct tensor inverse division');K[tuple(shape)]=num//denominator;need(K[tuple(shape)]>=0,'recovered combinatorial count nonnegative')
  for shape in shapes:
   reconstructed=sum(math.prod(math.factorial(k)*second[n][k] for n,k in zip(shape,ks))*K[ks] for ks in product(*(range(1,n+1) for n in shape)));need(reconstructed==allrecords[points[tuple(shape)]]['H'],'every physical count reconstructed by complete forward tensor transform')
  vals=[K[tuple(s)] for s in shapes];expected=data[mi]['queries'][qi]['target' if side else 'source'][:55] if kind=='mod8' else even['queries'][qi]['target_counts' if side else 'source_counts'];need(vals==expected,'all exact certifying run counts independently recovered');channels[key]=K;recovered.append(dict(kind=kind,model_index=mi,query_index=qi,side=side,query=data[mi]['queries'][qi]['query'] if kind=='mod8' else even['queries'][qi]['query'],shapes=shapes,recovered_counts=vals,physical_record_indices=[points[tuple(s)] for s in shapes]))
 for qi in range(2):
  S=channels['even',0,qi,0];T=channels['even',0,qi,1];gap=sum((-2)**(sum(shape)-9)*(S[shape]-T[shape]) for shape in S)%16;need(gap==8*qi,'independent physical geometric contrast')
  for i in range(9):shape=tuple(2 if j==i else 1 for j in range(9));need((S[shape]-T[shape])%2==int(i in [4,7]),'independent physical derivative label positions')
 outputs={'literal_interpolation_records.json':allrecords,'signed_stirling_recovery.json':recovered}
 for name,obj in outputs.items():
  file=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if file.exists():need(file.read_text()==txt,'complete independent interpolation data stable')
  else:file.write_text(txt)
 return dict(status='PASS',literal_physical_graphs=len(allrecords),physical_mod8_graphs=2640,physical_mod16_graphs=880,interpolation_groups=len(groups),exact_recovered_counts=sum(len(r['recovered_counts']) for r in recovered),maximum_augmented_order=max(len(r['augmented_rows']) for r in allrecords),source_pins={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['interpolation_audit.py','literal.cpp']},output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs},producer_proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
