from pathlib import Path
import json,hashlib,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'prime_eight_uniform_anchor'
def need(c,m):
 if not c:raise RuntimeError(m)
def tt(n):return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]
def op(T):return [((1<<len(T))-1)^(1<<i)^r for i,r in enumerate(T)]
def blow(Q,blocks):
 offs=[];T=[]
 for B in blocks:offs.append(len(T));T.extend([0]*len(B))
 for i,B in enumerate(blocks):
  for u,r in enumerate(B):T[offs[i]+u]=(r<<offs[i])|sum(((1<<len(C))-1)<<offs[j] for j,C in enumerate(blocks) if Q[i]>>j&1)
 return T,offs
def ports(T,a,c):
 n=len(T);R=T.copy()+[1<<(n+1),0]
 for v in range(n):
  for bits,p in[(a,n),(c,n+1)]:
   if bits>>v&1:R[v]|=1<<p
   else:R[p]|=1<<v
 return R
def uniform(core,blocks,a,c,outer=None):
 offs=[];n=0
 for B in blocks:offs.append(n);n+=len(B)
 aa=sum(((1<<len(B))-1)<<offs[i] for i,B in enumerate(blocks) if a>>i&1);cc=sum(((1<<len(B))-1)<<offs[i] for i,B in enumerate(blocks) if c>>i&1)
 if outer:
  core,off=blow(outer,[core]+[[0] for _ in outer[1:]]);inc=sum(1<<off[i] for i in range(1,len(outer)) if outer[i]&1);aa|=inc;cc|=inc
 return ports(core,aa,cc)
def main():
 src=json.loads((SRC/'normal.json').read_text());Q=src['quotient'];qs=[(1,0),(255,254)];cases=[];expected=[];labels=[]
 for r in src['large_physical_controls']:
  T,_=blow(Q,r['blocks']);U,_=blow(op(Q),r['blocks'])
  for a,c in qs:
   for side,R in enumerate([T,U]):cases.append(uniform(R,r['blocks'],a,c,r['ambient']));labels.append([r['name'],a,c,side])
  expected+=r['values']
 blocks=[tt(2)]*8;T,_=blow(Q,blocks);U,_=blow(op(Q),blocks)
 for r in src['full_named_batteries']:
  for j,(a,c) in enumerate(r['queries']):
   for side,R in enumerate([T,U]):cases.append(ports(R,a,c));labels.append([r['name'],j,side]);expected.append(r['values'][j][side])
 blocks=[tt(3)]+[tt(2)]*7;T,_=blow(Q,blocks);U,_=blow(op(Q),blocks)
 for a,c in qs:
  for side,R in enumerate([T,U]):cases.append(uniform(R,blocks,a,c));labels.append(['family_baseline',a,c,side])
 expected+=src['root_family_baseline_values']
 for r in src['root_family_controls']:
  for a,c in qs:cases.append(uniform(r['target_rows'],blocks,a,c));labels.append(['root_family',r['selected_bits'],a,c])
  expected+=r['values']
 lines=[]
 for R in cases:lines.extend([str(len(R)),' '.join(map(str,R))])
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'physical';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(exe)],check=True);got=list(map(int,subprocess.check_output([str(exe)],input=('\n'.join(lines)+'\n').encode()).split()))
 need(got==expected,'every exact large/ambient/nonuniform/family physical value matches')
 print(json.dumps({'status':'PASS','physical_count':len(got),'maximum_order':max(map(len,cases)),'records':[[label,v,hashlib.sha256(repr(R).encode()).hexdigest()] for label,v,R in zip(labels,got,cases)],'source_receipt_sha256':hashlib.sha256((SRC/'normal.json').read_bytes()).hexdigest(),'audit_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['physical.py','physical.cpp']}},sort_keys=True,indent=2))
if __name__=='__main__':main()
