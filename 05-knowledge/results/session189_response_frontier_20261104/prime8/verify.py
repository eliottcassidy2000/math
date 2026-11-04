from pathlib import Path
import json,subprocess,hashlib
from prepare_quotients import prime
ROOT=Path(__file__).resolve().parent

def need(c,m):
 if not c:raise RuntimeError(m)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run_python(name):return subprocess.check_output(['python3',str(ROOT/name)],cwd=ROOT).decode().strip()
def compile_cpp(name):subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(ROOT/(name+'.cpp')),'-o',str(ROOT/name)],check=True)
def relabel(T,p):return [sum(((T[p[i]]>>p[j])&1)<<j for j in range(len(T))) for i in range(len(T))]

def main():
 logs=[run_python('make_cover.py'),run_python('prepare_quotients.py'),run_python('anchor_explore.py')]
 C=json.loads((ROOT/'cover.json').read_text())['levels'];need([len(c['representatives']) for c in C]==[1,1,1,2,4,12,56,456,6880],'complete canonical level sizes');need(C[0]['representatives']==[[]],'empty base')
 total=0
 for n in range(1,9):
  previous=C[n-1]['representatives'];current=C[n]['representatives'];records=C[n]['extension_records'];need(len(records)==len(previous)*(1<<(n-1)),'all rooted extensions');seen=set()
  for k,(owner,packed) in enumerate(records):
   parent=k>>(n-1);mask=k&((1<<(n-1))-1);T=previous[parent];G=[T[i]|((1-((mask>>i)&1))<<(n-1)) for i in range(n-1)]+[mask];p=[(packed>>(3*i))&7 for i in range(n)];need(sorted(p)==list(range(n)) and packed>>(3*n)==0,'actual permutation encoding');need(0<=owner<len(current) and relabel(G,p)==current[owner],'edgewise extension isomorphism');seen.add(owner);total+=1
  need(seen==set(range(len(current))),'every retained representative has extension witness')
 R=json.loads((ROOT/'quotients.json').read_text());P=[(i,T) for i,T in enumerate(C[-1]['representatives']) if prime(T)];need(len(P)==4008 and [(r['cover_index'],r['rows']) for r in R]==P,'full prime filter, not a chosen subfamily')
 compile_cpp('endpoints');inp=str(len(R))+'\n'+''.join(' '.join(map(str,r['rows']))+'\n' for r in R);end=subprocess.check_output([str(ROOT/'endpoints')],input=inp.encode()).decode();(ROOT/'endpoint_independent.txt').write_text(end);E=[list(map(int,l.split())) for l in end.splitlines()];need(len(E)==4008,'all independent endpoint rows')
 for r,e in zip(R,E):need(e[0]==1 and e[1:9]==r['S'] and e[9:]==r['E'] and r['p']==[(x+y)%2 for x,y in zip(r['S'],r['E'])],'64,128 exact permutation endpoint coordinates')
 zero=[r for r in R if not any(r['p'])];need([r['rows'] for r in zero]==[[4,17,42,131,45,75,31,119],[16,5,41,67,46,139,55,95]],'exact two residuals');need(all('sigma' not in r for r in R if any(r['p'])),'do not retain partial sigma fields')
 known=[238,4,240,118,195,146,34,74];op=[255^(1<<i)^r for i,r in enumerate(known)];maps=[[1,6,5,7,2,4,3,0],[0,3,4,2,7,5,6,1]];need(relabel(known,maps[0])==zero[1]['rows'] and relabel(op,maps[1])==zero[0]['rows'],'known obstruction converse-pair maps')
 compile_cpp('anchor_check');inp='2\n'+''.join(' '.join(map(str,r['rows']))+'\n' for r in zero);raw=subprocess.check_output([str(ROOT/'anchor_check')],input=inp.encode()).decode();(ROOT/'anchor_independent.txt').write_text(raw);ind=[list(map(int,l.split())) for l in raw.splitlines()];A=json.loads((ROOT/'anchor_exploration.json').read_text());need(len(A)==2 and all(not a['blind'] for a in A),'all512 singleton-position cases')
 for a,iv in zip(A,ind):need(iv==sum(a['start_run_skews'],[])+a['one_port_base']+sum(a['one_port_double'],[]),'all272 exact independent anchor primitives')
 specifications=[[(0,5,0,20),(0,7,1,104),(1,0,0,124)],[(0,0,1,22),(0,2,0,40),(1,7,0,62)]]
 for a,spec in zip(A,specifications):
  need(spec[0][2]^spec[1][2]^spec[2][2]==1 and spec[0][3]^spec[1][3]^spec[2][3]==0,'three affine bits sum to one')
  for kind,i,b,mask in spec:
   field='second' if kind else 'first';scale=4 if kind else 2
   for pattern in a['patterns']:need(pattern[field][i]//scale==(b^((pattern['V']&mask).bit_count()%2)),'explicit three-query affine identity')
 physical_log=run_python('physical_checks.py');logs.append(physical_log);F=json.loads((ROOT/'physical_records.json').read_text());need(F['status']=='PASS' and F['physical_cases']==6548 and F['physical_order_max']==20,'physical result scope')
 for p in F['anchor_probes']+F['ambient_probes']:
  for side,T in enumerate(p['augmented_rows']):
   n=len(T);a,c=p['physical_queries'][side];need(p['port_names']==[n-2,n-1] and 0<=a<1<<(n-2) and 0<=c<1<<(n-2),'physical port metadata bounds');need(T[n-2]==(((1<<(n-2))-1)^a)|(1<<(n-1)) and T[n-1]==(((1<<(n-2))-1)^c),'physical retained port rows');need(all(((T[i]>>(n-2))&1)==((a>>i)&1) and ((T[i]>>(n-1))&1)==((c>>i)&1) for i in range(n-2)),'physical retained input rows')
   children=p['child_sets'];used=sorted(v for block in children for v in block);need(used==list(range(len(used))) and [len(block)>1 for block in children]==[bool(p['V']>>i&1) for i in range(8)],'physical child sets and singleton bits');exterior=p.get('exterior_names',[]);need(sorted(used+exterior)==list(range(n-2)),'physical old-name partition')
   Q=A[p['ri']]['rows'];Q=Q if side==0 else [255^(1<<i)^r for i,r in enumerate(Q)];need(all(((T[u]>>v)&1)==((Q[i]>>j)&1) for i in range(8) for j in range(8) if i!=j for u in children[i] for v in children[j]),'physical quotient arcs')
 data=['cover.json','quotients.json','anchor_exploration.json','endpoint_independent.txt','anchor_independent.txt','physical_records.json'];sources=['PROOF.md','verify.py','make_cover.py','prepare_quotients.py','core.py','anchor_explore.py','endpoints.cpp','anchor_check.cpp','physical_h.cpp','physical_checks.py']
 def pins(files):return {f:dict(bytes=(ROOT/f).stat().st_size,sha256=sha(ROOT/f)) for f in files}
 print(json.dumps(dict(status='PASS',scope='All prime-eight quotient converses; arbitrary nonempty children, equal-profile target replacements, every fixed ambient. Explicit source-relative original-module root-family consequence; no unrestricted reconstruction.',extension_records_checked=total,order_eight_rooted_extensions=58368,canonical_classes=6880,prime_classes=4008,endpoint_zero_classes=2,independent_exact_endpoint_coordinates=64128,independent_exact_anchor_primitives=272,singleton_position_patterns=512,three_query_affine_certificates=specifications,full_named_physical_response_values=262144,physical_cases=6548,cyclic_ambient_controls=8,source_pins=pins(sources),data_pins=pins(data)),sort_keys=True,indent=2))
if __name__=='__main__':main()
