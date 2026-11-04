from pathlib import Path
from array import array
from itertools import product
import hashlib,json,subprocess,sys
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(name):return hashlib.sha256((D/name).read_bytes()).hexdigest()
def rank(rows):
 piv={}
 for v in rows:
  while v:
   j=v.bit_length()-1
   if j not in piv:piv[j]=v;break
   v^=piv[j]
 return len(piv)
def main():
 for name in ('explore.py','build_banks.py','physical_checks.py'):subprocess.check_output([sys.executable,str(D/name)])
 from explore import mods
 from physical_checks import build
 data=json.loads((D/'families.json').read_text());need(len(data)==18,'exact stated finite family');need(len({tuple(q['rows']) for q in data})==18,'distinct displayed physical orientations')
 for q in data:
  rows=q['rows'];need(all(not(r>>i&1) and all(((r>>j)&1)+((rows[j]>>i)&1)==1 for j in range(9) if i!=j) for i,r in enumerate(rows)),'actual quotient tournament')
  need(q['nontrivial_modules']==mods(rows),'complete proper nontrivial module test')
  need(len(q['patterns'])==512,'all singleton patterns')
 orbit0=[0,3,5,10,12,15];orbit1=[1,4,7,8,11,14];transports=[]
 for ids in (orbit0,orbit1):
  src=data[ids[0]]
  for target in ids:
   dst=data[target];u=next(u for u in (1,2,4,5,7,8) if sorted(d*u%9 for d in src['construction']['steps'])==sorted(dst['construction']['steps']))
   p=[i*u%9 for i in range(9)];need(all(((src['rows'][i]>>j)&1)==((dst['rows'][p[i]]>>p[j])&1) for i in range(9) for j in range(9)),'physical isomorphism edge-by-edge');transports.append(dict(source=src['name'],target=dst['name'],source_to_target=p))
 for i,q in enumerate(data[:16]):
  need(q['prime']==(i in orbit0+orbit1),'exact12 prime cases')
  need(not any(q['p']) and not any(q['sigma']),'circulant endpoint flatness')
  need(q['anchor_blind_patterns']==([0,511] if q['prime'] else [0,73,146,219,292,365,438,511]),'exact blind patterns')
  if not q['prime']:need(q['nontrivial_modules']==[73,146,292],'composite residue-class modules')
 for V in range(512):
  v=lambda j:(V>>(j%9))&1
  for i in range(9):
   need(data[0]['patterns'][V]['first'][i]==2*(v(i+4)^v(i+5)),'carousel first-anchor cycle')
   need(data[1]['patterns'][V]['second'][i]==4*(v(i-1)^v(i+1)),'second circulant second-anchor cycle')
  w=data[17]['patterns'][V];need(w['first'][8]==2*(v(0)^v(7)) and w['second'][2]==4*(1^v(0)^v(7)),'W9 fixed two-query contradiction')
 need(data[16]['p']==[0,0,0,0,1,0,0,0,1] and data[16]['prime'] and data[17]['prime'],'displayed critical cases')
 need(sum(data[0]['S'])==3267 and sum(data[1]['S'])==3357,'different Hamilton counts separate circulant classes')
 matrices=[]
 for i in (0,1,2,6):
  q=data[i];first=[sum((x%2)<<j for j,x in enumerate(row)) for row in q['one_port_double']];second=[sum(((x+y)%2)<<j for j,(x,y) in enumerate(zip(r,s))) for r,s in zip(q['one_port_double'],q['start_run_skews'])];matrices.append(dict(name=q['name'],first_rank=rank(first),second_rank=rank(second),combined_rank=rank(first+second)))
 need([(m['first_rank'],m['second_rank'],m['combined_rank']) for m in matrices]==[(8,0,8),(6,8,8),(6,0,6),(6,0,6)],'anchor matrix ranks')
 (D/'transports.json').write_text(json.dumps(transports,indent=2)+'\n')
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'endpoint_check.cpp'),'-o',str(D/'endpoint_check')],check=True)
 inp=str(len(data))+'\n'+''.join(' '.join(map(str,q['rows']))+'\n' for q in data);raw=subprocess.check_output([str(D/'endpoint_check')],input=inp.encode()).decode();(D/'endpoint_values.txt').write_text(raw)
 for q,line in zip(data,raw.splitlines()):need(list(map(int,line.split()))==sum(q['joint_endpoints'],[]),'independent permutation joint endpoints')
 banks=json.loads((D/'named_banks.json').read_text());need(len(banks)==8,'all source and target complete named banks');cube_hashes=[];counts={d:0 for d in range(5)}
 for b in banks:
  masks=b['physical_masks'];need(masks==[m for m in range(512) if m.bit_count()<=2],'complete bidegree2 basis');cube=[0]*(1<<18)
  for i,a in enumerate(masks):
   for j,c in enumerate(masks):cube[a+512*c]=b['coefficients'][i][j];counts[a.bit_count()+c.bit_count()]+=1
  for bit in range(18):
   step=1<<bit
   for block in range(0,1<<18,2*step):
    for j in range(block+step,block+2*step):cube[j]+=cube[j-step]
  for i,a in enumerate(masks):
   for j,c in enumerate(masks):need(cube[a+512*c]==b['basis_values'][i][j],'Möbius reconstruction all retained basis values')
  need(all(v>0 and v%2 and v<=39916800 for v in cube),'entire reconstructed cube Hamilton bounds and parity')
  arr=array('I',cube)
  if sys.byteorder!='little':arr.byteswap()
  cube_hashes.append(dict(name=b['name'],side=b['side'],full_values=1<<18,little_endian_u32_sha256=hashlib.sha256(arr.tobytes()).hexdigest()))
 for left,right in zip(banks[::2],banks[1::2]):
  if left['name'].startswith('circulant'):
   pos={m:i for i,m in enumerate(left['physical_masks'])};gap=left['basis_values'][pos[1]][pos[2]]-right['basis_values'][pos[1]][pos[2]];need(gap%4==2,'singleton residual gap2mod4')
 records=json.loads((D/'physical_records.json').read_text())
 for r in records:
  for side,s in enumerate(r['sides']):need({k:v for k,v in s.items() if k!='H'}==build(r['quotient_rows'],r['source_child_rows'],*r['block_query'],side,r['ambient']),'physical row and mask transport')
  need(r['contrast']==r['sides'][0]['H']-r['sides'][1]['H'] and r['contrast']%r['modulus']==r['expected_residue'],'actual aligned physical contrast')
 (D/'cube_hashes.json').write_text(json.dumps(cube_hashes,indent=2)+'\n')
 names=['PROOF.md','explore.py','anchors.cpp','endpoint_check.cpp','build_banks.py','physical_h.cpp','physical_checks.py','families.json','anchor_values.txt','endpoint_values.txt','transports.json','named_values.txt','named_banks.json','cube_hashes.json','physical_records.json','physical_values.txt']
 print(json.dumps(dict(status='PASS',scope='16 cyclic step choices plus2 explicitly displayed prime9 quotients; not a prime9 census',models=18,prime_circulants=12,unit_multiplier_families=2,complete_joint_endpoint_coordinates=18*81,singleton_patterns=18*512,anchor_ranks=matrices,complete_named_bank_coordinates=8*2116,named_coefficient_counts=counts,reconstructed_full_values=8*(1<<18),literal_physical_sides=2*len(records),sha256={name:digest(name) for name in names}),sort_keys=True,separators=(',',':')))
if __name__=='__main__':main()
