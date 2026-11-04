from pathlib import Path
from itertools import product
import json,subprocess
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(Q):return [((1<<len(Q))-1)^(1<<i)^r for i,r in enumerate(Q)]
def aug(Q,a,c):
 n=len(Q);return [r|(((a>>i)&1)<<n)|(((c>>i)&1)<<(n+1)) for i,r in enumerate(Q)]+[(((1<<n)-1)^a)|(1<<(n+1)),((1<<n)-1)^c]
def subsets(m):
 s=m
 while True:
  yield s
  if not s:break
  s=(s-1)&m

def main():
 families=json.loads((D/'anchor_data.json').read_text());names=['node860','node2005','node4067','node6501'];families=[dict(q,name=name) for q,name in zip(families,names)];models=[];cases=[];masks=[m for m in range(512) if m.bit_count()<=2]
 for name in names:
  q=next(r for r in families if r['name']==name)
  for side in range(2):
   rows=op(q['rows']) if side else q['rows'];models.append(dict(name=name,side=side,rows=rows,physical_masks=masks,port_names=[9,10],port_edge='9->10',a_convention='old vertex points to9',c_convention='old vertex points to10'))
   cases.extend(aug(rows,a,c) for a,c in product(masks,repeat=2))
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical_h.cpp'),'-o',str(D/'physical_h')],check=True)
 inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases)
 raw=subprocess.check_output([str(D/'physical_h')],input=inp.encode()).decode();(D/'named_values.txt').write_text(raw);values=list(map(int,raw.split()));need(len(values)==len(cases),'all exact basis counts');pos={m:i for i,m in enumerate(masks)};at=0
 for m in models:
  m['basis_values']=[values[at+46*i:at+46*(i+1)] for i in range(46)];at+=2116
  m['coefficients']=[[sum((-1)**(a.bit_count()+c.bit_count()-u.bit_count()-v.bit_count())*m['basis_values'][pos[u]][pos[v]] for u in subsets(a) for v in subsets(c)) for c in masks] for a in masks]
 (D/'named_banks.json').write_text(json.dumps(models,separators=(',',':'))+'\n')
 summary=[]
 for left,right in zip(models[::2],models[1::2]):
  gap=[(a,c,left['basis_values'][i][j]-right['basis_values'][i][j]) for i,a in enumerate(masks) for j,c in enumerate(masks) if left['basis_values'][i][j]!=right['basis_values'][i][j]]
  best=min(gap,key=lambda r:((abs(r[2])&-abs(r[2])).bit_length(),r[0].bit_count()+r[1].bit_count(),r[:2]))
  summary.append(dict(name=left['name'],best_basis_gap=best,nonzero_basis_gaps=len(gap)))
 print(json.dumps(dict(status='PASS',complete_named_banks=8,named_basis_values=len(cases),coefficient_coordinates=len(cases),summaries=summary),indent=2))
if __name__=='__main__':main()
