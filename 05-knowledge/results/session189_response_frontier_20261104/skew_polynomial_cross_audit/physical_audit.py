"""Independent compression recurrence, literal controls, and every small named coefficient."""
from pathlib import Path
import json,hashlib,math,tempfile,subprocess,collections
from polynomial_audit import inflate,op
D=Path(__file__).resolve().parent;P=D.parent/'skew_polynomial'
if not P.exists():P=D.parent/'session189-skew-polynomial'
def need(c,s):
 if not c:raise RuntimeError(s)
def aug(R,a,c):
 n=len(R);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(R)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]
def cyclic(R):
 n=len(R);return [r|(1<<n) for r in R]+[1<<(n+1),(1<<n)-1]
def main():
 selected=json.loads((P/'selected.json').read_text());byid={q['node']:q for q in selected};weighted=json.loads((P/'selected.out').read_text());byw={q['node']:w for q,w in zip(selected,weighted)};physical=json.loads((P/'physical_queries.json').read_text());queries=physical['queries'];oldvals=list(map(int,(P/'physical_values.txt').read_text().split()));need(len(queries)==len(oldvals)==208,'complete208 physical-query records');used=set()
 for case in physical['cases']:
  q=byid[case['node']];m=case['m'];g=case['group'];ambient=case['ambient'];need(g==q['certificate_group'],'same selected group');old=m+16;types=m+8;group=[m+2*(g-1),m+2*(g-1)+1];need(group==case['group_names'],'fixed physical group names');expected_terms=[(m-1,1)]+([(2,3)] if m==5 else []);need([(t['r'],t['weight']) for t in case['terms']]==expected_terms,'correct odd binomial weights')
  for t in case['terms']:
   need(len(t['queries'])==16,'complete finite-difference corners')
   for corner,x in enumerate(t['queries']):
    aa,cc=divmod(corner,4);a=(aa&1)|((aa>>1&1)<<1);c=((1<<t['r']) if cc&1 else 0)|(sum(1<<v for v in group) if cc&2 else 0);ca=a;cy=((1<<t['r']) if cc&1 else 0)|((1<<(m+g-1)) if cc&2 else 0);cap=[1]*m+[2]*8
    if ambient:a|=1<<(old+1);c|=1<<(old+1);ca|=1<<(types+1);cy|=1<<(types+1);cap.extend([1,1])
    need(x['a']==a and x['c']==c and x['sign']==(-1)**(4-aa.bit_count()-cc.bit_count()),'every actual physical mask and finite-difference sign')
    for side,key in [(0,'source'),(1,'target')]:
     Q=op(q['rows']) if side else q['rows'];actual,_=inflate(Q,[m]+[2]*8);compressed,_=inflate(Q,[m]+[1]*8)
     if ambient:actual=cyclic(actual);compressed=cyclic(compressed)
     actual=aug(actual,a,c);compressed=aug(compressed,ca,cy);r=queries[x[key]];used.add(x[key]);need(r==dict(compressed_rows=compressed,capacities=cap+[1,1],physical_rows=actual,a_mask=a,c_mask=c),'entire source/quotient-target physical query and compressed representation')
 need(used==set(range(208)),'every retained query belongs to the declared battery')
 compression=json.loads((P/'compression_controls.json').read_text());need(len(compression)==100,'all100 independent literal controls')
 for r in compression:need(inflate(r['rows'],r['capacities'])[0]==r['physical_rows'],'complete actual compression graph')
 bank=json.loads((P/'small_bank.json').read_text());small=json.loads((P/'small_controls.json').read_text());cases=[]
 for b in bank:
  n=len(b['rows']);masks=[0]+[1<<i for i in range(n)]+[m for m in range(1<<n) if m.bit_count()==2];need(b['masks']==masks,'complete canonical physical named basis');cases.extend(aug(b['rows'],a,c) for a in masks for c in masks)
 case_count=len(cases);cases.extend(r['physical_rows'] for r in compression);need(max(map(len,cases))<=20 and math.factorial(20)<1<<64,'literal integer bound');need(math.factorial(max(len(q['physical_rows']) for q in queries))<1<<127,'large compressed exact integer bound')
 with tempfile.TemporaryDirectory() as tmp:
  direct=Path(tmp)/'literal';compressed=Path(tmp)/'compressed'
  for src,dst in [('literal.cpp',direct),('compressed_suffix.cpp',compressed)]:subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/src),'-o',str(dst)],check=True)
  inp=str(len(cases))+'\n'+''.join(str(len(R))+'\n'+' '.join(map(str,R))+'\n' for R in cases);literal=list(map(int,subprocess.check_output([str(direct)],input=inp.encode()).split()));need(len(literal)==len(cases),'every literal count returned')
  rc=[(q['compressed_rows'],q['capacities']) for q in queries]+[(q['rows'],q['capacities']) for q in compression];inp=str(len(rc))+'\n'+''.join(str(len(R))+'\n'+' '.join(map(str,R))+'\n'+' '.join(map(str,cap))+'\n' for R,cap in rc);compressed_values=list(map(int,subprocess.check_output([str(compressed)],input=inp.encode()).split()));need(len(compressed_values)==308,'all independently compressed counts')
 need(compressed_values[:208]==oldvals,'all208 large actual query counts by independent completion method')
 control=[]
 for r,H,Hc in zip(compression,literal[case_count:],compressed_values[208:]):need(H==Hc==r['H'],'every direct physical control equals exact compressed count');control.append(dict(**r,independent_H=H))
 grouped=[]
 for c in physical['cases']:
  gap=sum(t['weight']*sum(x['sign']*(compressed_values[x['source']]-compressed_values[x['target']]) for x in t['queries']) for t in c['terms']);k=byw[c['node']]['skews'][c['group']];want=2*(3*k[0]+k[1]) if c['ambient'] else 2*k[0];need(gap==want and gap>0,'full physical grouped quotient contrast and fixed ambient transfer');grouped.append(dict(node=c['node'],m=c['m'],ambient=c['ambient'],actual_contrast=gap,expected=want))
 at=0;newbank=[];degrees=collections.Counter()
 for b in bank:
  masks=b['masks'];L=len(masks);v=literal[at:at+L*L];at+=L*L;need(v==b['values'],'all independently recounted named basis values');idx={m:i for i,m in enumerate(masks)};co=[]
  for a in masks:
   for c in masks:
    result=0;u=a
    while True:
     w=c
     while True:
      result+=(-1)**(a.bit_count()+c.bit_count()-u.bit_count()-w.bit_count())*v[idx[u]*L+idx[w]]
      if not w:break
      w=(w-1)&c
     if not u:break
     u=(u-1)&a
    co.append(result);degrees[a.bit_count()+c.bit_count()]+=1
  need(co==b['coefficients'],'every named coefficient including zero cubic and quartic entries');newbank.append(dict(id=b['id'],rows=b['rows'],masks=masks,values=v,coefficients=co))
 need(at==15219,'complete small named coordinate count')
 for c in small['cases']:
  source,target,local=c['graph_ids'];m=c['m'];Q=c['Q'];caps=[m]+([2,2] if len(Q)==3 else [1]*4);S,_=inflate(Q,caps);U,_=inflate(op(Q),caps);V=S.copy()
  for i in range(m):V[i]=(V[i]&~((1<<m)-1))|op(inflate([0],[m])[0])[i]
  if c['ambient']:S,U,V=map(cyclic,[S,U,V])
  need([S,U,V]==[newbank[i]['rows'] for i in c['graph_ids']],'all small source, quotient target and local-converse physical graphs')
 def coeff(g,a,c):
  b=newbank[g];L=len(b['masks']);ix={m:i for i,m in enumerate(b['masks'])};return b['coefficients'][ix[a]*L+ix[c]]
 small_results=[];sw=json.loads((P/'small_weighted.out').read_text())
 for c,k in zip(small['cases'],sw):
  s,t,l=c['graph_ids'];m=c['m'];weights=[(m-1,1)]+([(2,3)] if m==5 else [])
  for g,G in enumerate(c['groups'],1):
   ds=sum(w*sum(coeff(s,3,(1<<r)|(1<<i))-coeff(t,3,(1<<r)|(1<<i)) for i in G) for r,w in weights);dl=sum(w*sum(coeff(s,3,(1<<r)|(1<<i))-coeff(l,3,(1<<r)|(1<<i)) for i in G) for r,w in weights);pol=k['skews'][g];pred=2*(3*pol[0]+pol[1]) if c['ambient'] else 2*pol[0];need(ds==dl==pred,'full named quotient/local restoration and grouped ambient formula');small_results.append(dict(source_graph=s,group=G,quotient_contrast=ds,local_contrast=dl,polynomial_prediction=pred))
 outputs={'independent_small_banks.json':newbank,'physical_query_recount.json':dict(queries=[dict(**q,H=H) for q,H in zip(queries,compressed_values)],cases=grouped),'compression_recount.json':control,'small_grouped_controls.json':small_results}
 for name,obj in outputs.items():
  p=D/name;txt=json.dumps(obj,separators=(',',':'))+'\n'
  if p.exists():need(p.read_text()==txt,'independent exact physical artifacts stable')
  else:p.write_text(txt)
 return dict(status='PASS',large_actual_queries=208,maximum_large_physical_order=max(len(q['physical_rows']) for q in queries),independent_literal_compression_controls=100,complete_small_named_graphs=len(newbank),complete_small_named_coordinates=at,named_degree_counts=dict(sorted(degrees.items())),source_pins={f:hashlib.sha256((D/f).read_bytes()).hexdigest() for f in ['physical_audit.py','literal.cpp','compressed_suffix.cpp','polynomial_audit.py']},output_pins={name:hashlib.sha256((D/name).read_bytes()).hexdigest() for name in outputs},producer_proof_sha256=hashlib.sha256((P/'PROOF.md').read_bytes()).hexdigest())
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
