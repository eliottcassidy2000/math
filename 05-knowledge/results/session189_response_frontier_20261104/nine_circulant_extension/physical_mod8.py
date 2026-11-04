from pathlib import Path
from itertools import combinations
import json,subprocess
from physical_even import build,op,need
D=Path(__file__).resolve().parent

def main():
 data=json.loads((D/'analysis.json').read_text());certs=json.loads((D/'linear_certificates.json').read_text());records=[];cases=[];groups=[]
 shapes=[[1]*9]
 for k in (2,3):
  for i in range(9):s=[1]*9;s[i]=k;shapes.append(s)
 for i,j in combinations(range(9),2):s=[1]*9;s[i]=s[j]=2;shapes.append(s)
 for mi,(q,cert) in enumerate(zip(data,certs)):
  c=cert['safe_certificate'];ids=c.get('contradiction_query_indices',c['independent_query_indices'])
  for qi in ids:
   query=q['queries'][qi]
   for side in range(2):
    Q=op(q['rows']) if side else q['rows'];start=len(records)
    for si,sizes in enumerate(shapes):
     r=build(Q,sizes,*query['query']);r.update(model_index=mi,query_index=qi,side=side,shape_index=si);records.append(r);cases.append(r['augmented_rows'])
    groups.append((mi,qi,side,start))
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical_h.cpp'),'-o',str(D/'physical_h')],check=True);inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases);raw=subprocess.check_output([str(D/'physical_h')],input=inp.encode()).decode();(D/'physical_mod8_values.txt').write_text(raw);values=list(map(int,raw.split()));need(len(values)==len(cases),'all selected certificate physical models')
 for r,H in zip(records,values):r['H']=H
 for mi,qi,side,start in groups:
  F=values[start:start+55];k=[F[0]]
  for i in range(9):need((F[i+1]-F[0])%2==0,'single double-run exact divisibility');k.append((F[i+1]-F[0])//2)
  for i in range(9):num=F[10+i]-3*F[1+i]+2*F[0];need(num%6==0,'triple-run exact divisibility');k.append(num//6)
  for index,(i,j) in enumerate(combinations(range(9),2)):
   num=F[19+index]-F[1+i]-F[1+j]+F[0];need(num%4==0,'two-double-run exact divisibility');k.append(num//4)
  need(k==data[mi]['queries'][qi]['target' if side else 'source'][:55],'full selected certificate run bank recovered from physical names')
 (D/'physical_mod8_records.json').write_text(json.dumps(records,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',selected_query_pairs=len(groups)//2,physical_graphs=len(cases),exact_recovered_run_counts=len(cases),maximum_augmented_order=max(map(len,cases)))))
if __name__=='__main__':main()
