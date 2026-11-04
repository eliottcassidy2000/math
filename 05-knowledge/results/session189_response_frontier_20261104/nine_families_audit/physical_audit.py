"""Independent full named banks and all-cover safety on literal matched-profile models."""
from pathlib import Path
from array import array
import subprocess,tempfile,json,hashlib,math,sys
D=Path(__file__).resolve().parent;P=D.parent/'nine_families'
if not P.exists():P=D.parent/'session189-nine-families'
def need(c,s):
 if not c:raise RuntimeError(s)
def op(R):return [((1<<len(R))-1)^(1<<i)^r for i,r in enumerate(R)]
def augment(R,a,c):
 n=len(R);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(R)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]
def sub(Q,children):
 rows=[];sets=[]
 for C in children:sets.append(list(range(len(rows),len(rows)+len(C))));rows.extend([0]*len(C))
 for i,C in enumerate(children):
  for j,r in enumerate(C):rows[sets[i][j]]=(r<<sets[i][0])|sum(1<<v for k,S in enumerate(sets) if Q[i]>>k&1 for v in S)
 return rows,sets

def small_profile(R):
 n=len(R);dp=[[0]*n for _ in range(1<<n)];total=[0]*(1<<n);B=64
 for S in range(1,1<<n):
  for v in range(n):
   if S>>v&1:
    M=S^(1<<v);dp[S][v]=1 if not M else (total[M]<<B)+sum(dp[M][u] for u in range(n) if M>>u&1 and R[u]>>v&1)
  total[S]=sum(dp[S])
 return [(total[-1]>>(B*j))&((1<<B)-1) for j in range(n)]
def main():
 banks=json.loads((P/'named_banks.json').read_text());family=json.loads((P/'families.json').read_text());byname={r['name']:r for r in family};need(len(banks)==8,'eight full physical named banks');masks=[m for m in range(512) if m.bit_count()<=2];cases=[]
 for b in banks:
  expected=byname[b['name']]['rows'];expected=op(expected) if b['side'] else expected;need(b['rows']==expected and b['physical_masks']==masks and b['port_names']==[9,10],'exact bank physical rows and ports')
  cases.extend(augment(expected,a,c) for a in masks for c in masks)
 records=json.loads((P/'physical_records.json').read_text());actual=[];meta=[];childprof={}
 for index,r in enumerate(records):
  Q=r['quotient_rows'];need(Q==byname[r['name']]['rows'],'correct finite family quotient');need(sum((len(C)>1)<<i for i,C in enumerate(r['source_child_rows']))==r['nontrivial_child_mask'],'actual singleton pattern matches claim')
  for dual in [False,True]:
   a,c=r['block_query'];a,c=(511-c,511-a) if dual else (a,c)
   for side in [0,1]:
    s=r['sides'][side];children=s['target_child_rows']
    for C,S in zip(children,r['source_child_rows']):
     for X in [C,S]:childprof.setdefault(tuple(X),small_profile(X))
     need(childprof[tuple(C)]==childprof[tuple(S)],'complete corresponding child profiles')
    old,sets=sub(op(Q) if side else Q,children);A=sum(1<<v for i,S in enumerate(sets) if a>>i&1 for v in S);C=sum(1<<v for i,S in enumerate(sets) if c>>i&1 for v in S);ext=[]
    if r['ambient']:
     n=len(old);old,_=sub([2,4,1],[old,[0],[0]]);A|=1<<(n+1);C|=1<<(n+1);ext=[n,n+1]
    rows=augment(old,A,C)
    if not dual:need(rows==s['augmented_rows'] and sets==s['child_sets'] and [A,C]==s['physical_query'] and ext==s['exterior_names'] and [len(old),len(old)+1]==s['port_names'],'every retained physical arc, set, port and exterior bit')
    actual.append(rows);meta.append(dict(record=index,dual=dual,side=side,augmented_rows=rows,child_sets=sets,physical_query=[A,C],exterior_names=ext,port_names=[len(old),len(old)+1]))
 # Independent full physical transport of the representative singleton probe.
 transports=json.loads((P/'transports.json').read_text());transport_cases=[];transport_meta=[]
 for t in transports:
  source=byname[t['source']]['rows'];target=byname[t['target']]['rows'];perm=t['source_to_target'];p=perm+[9,10];A=1<<perm[0];C=1<<perm[1]
  for side in [0,1]:
   R=op(source) if side else source;S=op(target) if side else target;ra=augment(R,1,2);sa=augment(S,A,C);need(all(((ra[i]>>j)&1)==((sa[p[i]]>>p[j])&1) for i in range(11) for j in range(11)),'entire augmented physical transport keeps both ports fixed')
   transport_cases.append(sa);transport_meta.append(dict(source=t['source'],target=t['target'],side=side,source_to_target=p,transported_query=[A,C],augmented_rows=sa))
 need(math.factorial(20)<1<<64,'exact Hamilton no-overflow bound')
 with tempfile.TemporaryDirectory() as tmp:
  binary=Path(tmp)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'literal.cpp'),'-o',str(binary)],check=True)
  def run(graphs,profile=False):
   inp=str(len(graphs))+'\n'+''.join(str(len(R))+'\n'+' '.join(map(str,R))+'\n' for R in graphs);raw=subprocess.check_output([str(binary)]+(['profiles'] if profile else []),input=inp.encode()).decode();vals=[list(map(int,line.split())) for line in raw.splitlines()];need(len(vals)==len(graphs),'every actual literal count returned');return vals
  values=run(cases);physical=run(actual,True);transport_values=run(transport_cases)
 newbanks=[];counts=[0]*5;hashes=[];at=0;ix={m:i for i,m in enumerate(masks)}
 for b in banks:
  val=[x[0] for x in values[at:at+2116]];at+=2116;matrix=[val[46*i:46*(i+1)] for i in range(46)];need(matrix==b['basis_values'],'all independent exact named basis values');coef=[]
  for a in masks:
   row=[]
   for c in masks:
    result=0;u=a
    while True:
     v=c
     while True:
      result+=(-1)**(a.bit_count()+c.bit_count()-u.bit_count()-v.bit_count())*matrix[ix[u]][ix[v]]
      if not v:break
      v=(v-1)&c
     if not u:break
     u=(u-1)&a
    row.append(result);counts[a.bit_count()+c.bit_count()]+=1
   coef.append(row)
  need(coef==b['coefficients'],'all complete named coefficients including zero cubic/quartic coordinates');cube=[0]*(1<<18)
  for i,a in enumerate(masks):
   for j,c in enumerate(masks):cube[a|(c<<9)]=coef[i][j]
  for bit in range(18):
   step=1<<bit
   for block in range(0,1<<18,2*step):
    for j in range(block+step,block+2*step):cube[j]+=cube[j-step]
  arr=array('I',cube)
  if sys.byteorder!='little':arr.byteswap()
  hashes.append(dict(name=b['name'],side=b['side'],full_values=1<<18,little_endian_u32_sha256=hashlib.sha256(arr.tobytes()).hexdigest()));newbanks.append(dict(name=b['name'],side=b['side'],rows=b['rows'],physical_masks=masks,port_names=[9,10],basis_values=matrix,coefficients=coef))
 need(hashes==json.loads((P/'cube_hashes.json').read_text()),'all reconstructed complete cube hashes')
 for i,(m,p) in enumerate(zip(meta,physical)):
  r=records[m['record']];n=len(m['augmented_rows']);m['H']=p[0];m['complete_profile_mod8']=[p[0]%8,p[1],p[2]]+[0]*(n-3)
  if not m['dual']:need(p[0]==r['sides'][m['side']]['H'],'all producer physical H counts')
 for j,r in enumerate(records):
  p=meta[4*j:4*j+4];delta=p[0]['H']-p[1]['H'];dual=p[2]['H']-p[3]['H'];need(delta==r['contrast'] and delta%r['modulus']==r['expected_residue'],'physical claimed anchor residue');need(dual==-delta,'fixed-exterior core-only converse-complement duality over integers')
  gap=[a-b for a,b in zip(p[0]['complete_profile_mod8'],p[1]['complete_profile_mod8'])];need(all(x%r['modulus']==0 for x in gap[1:]),'complete higher-cover safety for the literal anchor or singleton residual')
 for m,value in zip(transport_meta,transport_values):
  b=next(b for b in banks if b['name']==m['source'] and b['side']==m['side']);need(value[0]==b['basis_values'][ix[1]][ix[2]],'literal transported query value');m['H']=value[0]
 outputs={'independent_named_banks.json':newbanks,'cube_hashes.json':hashes,'physical_and_dual_controls.json':meta,'transport_controls.json':transport_meta}
 for name,obj in outputs.items():
  file=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if file.exists():need(file.read_text()==txt,'independent full bank and physical data stable')
  else:file.write_text(txt)
 return dict(status='PASS',full_named_banks=8,exact_named_basis_values=len(cases),complete_named_coefficients=len(cases),named_degree_counts=counts,reconstructed_complete_cube_values=8*(1<<18),producer_physical_sides=2*len(records),additional_dual_physical_sides=2*len(records),literal_singleton_transport_sides=len(transport_cases),complete_mod8_profiles=len(meta),maximum_augmented_order=max(map(len,actual)),source_pins={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['physical_audit.py','literal.cpp']},output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs},producer_proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
