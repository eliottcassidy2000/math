from pathlib import Path
import hashlib,json,subprocess,tempfile,sys
D=Path(__file__).resolve().parent;P=D.parent/'prime8'
if not P.exists():P=D.parent/'session189-prime8'
def need(x,s):
 if not x:raise RuntimeError(s)
def op(Q):return [((1<<len(Q))-1)^(1<<i)^r for i,r in enumerate(Q)]
def extend(T,a,c):
 n=len(T);R=T+[0,0];R[n]|=1<<(n+1)
 for i in range(n):
  if a>>i&1:R[i]|=1<<n
  else:R[n]|=1<<i
  if c>>i&1:R[i]|=1<<(n+1)
  else:R[n+1]|=1<<i
 return R
def model(Q,V,a,c,ambient):
 B=[];T=[]
 for i in range(8):B.append(list(range(len(T),len(T)+1+(V>>i&1))));T.extend([0]*len(B[-1]))
 for i,bs in enumerate(B):
  for u in bs:
   for v in bs:
    if u<v:T[u]|=1<<v
   for j,cs in enumerate(B):
    if Q[i]>>j&1:
     for v in cs:T[u]|=1<<v
 aa=sum(1<<v for i,bs in enumerate(B) if a>>i&1 for v in bs);cc=sum(1<<v for i,bs in enumerate(B) if c>>i&1 for v in bs)
 if ambient:
  n=len(T);T+=[0,0]
  for v in range(n):T[v]|=1<<n;T[n+1]|=1<<v
  T[n]|=1<<(n+1);aa|=1<<(n+1);cc|=1<<(n+1)
 return extend(T,aa,cc),aa,cc,B
def subsets(a):
 s=a
 while True:
  yield s
  if not s:break
  s=(s-1)&a

def main():
 R=json.loads((P/'physical_records.json').read_text());cases=[];expected=[];checks=[];phys=[];basis=R['basis_queries'];need(len(basis)==1369,'complete named basis')
 for ri,record in enumerate(R['records']):
  for side in range(2):
   Q=record['rows'] if side==0 else op(record['rows']);B=record['sides'][side];values={tuple(q):v for q,v in zip(basis,B['values'])};co=[]
   for a,c in basis:
    co.append(sum((-1)**(a.bit_count()+c.bit_count()-u.bit_count()-v.bit_count())*values[u,v] for u in subsets(a) for v in subsets(c)));cases.append(extend(Q,a,c));expected.append(values[a,c]);checks.append(['named_basis',ri,side,a,c])
   need(co==B['coefficients'],'all named coefficients including cubic quartic zeros');f=[0]*65536
   for (a,c),v in zip(basis,co):f[a|(c<<8)]=v
   for b in range(16):
    for mask in range(65536):
     if mask>>b&1:f[mask]+=f[mask^(1<<b)]
   need(f==B['complete_query_values_a_low_c_high'],'complete full physical cube retained exactly')
 for kind in ['anchor_probes','ambient_probes']:
  for p in R[kind]:
   rr=[]
   for side in range(2):
    Q=R['records'][p['ri']]['rows'];Q=Q if side==0 else op(Q);T,a,c,B=model(Q,p['V'],*p['query'],kind=='ambient_probes')
    need([a,c]==p['physical_queries'][side],'retained complete physical masks including ambient incidences');cases.append(T);expected.append(p['values'][side]);checks.append([kind,p['ri'],p['V'],side,a,c]);rr.append(dict(side=side,augmented_rows=T,core_child_physical_sets=B,actual_masks=[a,c],H=p['values'][side]))
   need((p['values'][0]-p['values'][1])%p['modulus']==p['predicted'],'actual modular separator');phys.append(dict(kind=kind,ri=p['ri'],V=p['V'],sides=rr))
 for h in R['nonbasis_literal_values']:
  Q=R['records'][h['ri']]['rows'];Q=Q if h['side']==0 else op(Q);a,c=h['query'];cases.append(extend(Q,a,c));expected.append(h['value']);checks.append(['nonbasis',h['ri'],h['side'],a,c])
 need(len(cases)==6548,'every physical query counted');need(max(map(len,cases))==20,'safe20! bound')
 with tempfile.TemporaryDirectory(prefix='prime8-physical-audit-') as td:
  e=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'literal.cpp'),'-o',str(e)],check=True)
  data=str(len(cases))+'\n'+''.join(str(len(T))+' '+' '.join(map(str,T))+'\n' for T in cases);actual=list(map(int,subprocess.check_output([str(e)],input=data.encode()).split()));need(actual==expected,'complete independent backward recurrence physical values')
 for name,obj in [('physical_probe_recount.json',phys),('physical_counts.json',dict(case_keys=checks,values=actual))]:
  s=json.dumps(obj,separators=(',',':'))+'\n';p=D/name
  if p.exists():need(p.read_text()==s,'audit physical artifact stable')
  else:p.write_text(s)
 return {'status':'PASS','independent_literal_counts':len(cases),'complete_named_basis_coordinates':5476,'complete_retained_full_cube_values':262144,'anchor_patterns':512,'ambient_patterns':8,'literal_order_max':20,'source_pins':{f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['physical_audit.py','literal.cpp']},'producer_physical_sha256':hashlib.sha256((P/'physical_records.json').read_bytes()).hexdigest(),'output_pins':{f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['physical_probe_recount.json','physical_counts.json']}}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
