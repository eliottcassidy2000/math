from pathlib import Path
import json,subprocess
from physical_helpers import build,op,need
D=Path(__file__).resolve().parent

def main():
 data=json.loads((D/'analysis.json').read_text());anchors=json.loads((D/'anchor_data.json').read_text());records=[];cases=[]
 for mi,V in [(0,254),(3,54)]:
  q=data[mi];a,c=1,2;r=next(r for r in q['queries'] if r['query']==[a,c]);delta=[s-t for s,t in zip(r['source'],r['target'])];need((delta[0]+2*sum(delta[i+1] for i in range(9) if V>>i&1))%4==2,'complete singleton/nontrivial mod4 run formula')
  for ambient in [False,True]:
   rec=dict(model_index=mi,quotient_rows=q['rows'],nontrivial_child_mask=V,block_query=[a,c],ambient=ambient,sides=[])
   for side in range(2):
    Q=op(q['rows']) if side else q['rows'];s=build(Q,[2 if V>>i&1 else 1 for i in range(9)],a,c);s['exterior_names']=[]
    if ambient:
     n=len(s['augmented_rows'])-2;oldmask=(1<<n)-1;T=[(row&oldmask)|(1<<n) for row in s['augmented_rows'][:-2]]+[1<<(n+1),oldmask];aa,cc=[v|(1<<(n+1)) for v in s['physical_query']];nn=n+2;s['augmented_rows']=[row|(((aa>>i)&1)<<nn)|(((cc>>i)&1)<<(nn+1)) for i,row in enumerate(T)]+[(((1<<nn)-1)^aa)|(1<<(nn+1)),((1<<nn)-1)^cc];s['physical_query']=[aa,cc];s['port_names']=[nn,nn+1];s['exterior_names']=[n,n+1]
    rec['sides'].append(s);cases.append(s['augmented_rows'])
   records.append(rec)
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'physical_h.cpp'),'-o',str(D/'physical_h')],check=True);inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases);raw=subprocess.check_output([str(D/'physical_h')],input=inp.encode()).decode();(D/'residual_values.txt').write_text(raw);values=list(map(int,raw.split()));need(len(values)==8,'every residual physical side')
 for j,r in enumerate(records):
  for side in range(2):r['sides'][side]['H']=values[2*j+side]
  r['contrast']=values[2*j]-values[2*j+1];need(r['contrast']%4==2,'literal local or ambient singleton residual gap')
 (D/'residual_records.json').write_text(json.dumps(records,separators=(',',':'))+'\n');print(json.dumps(dict(status='PASS',literal_physical_sides=8,maximum_augmented_order=max(map(len,cases)),singleton_residuals=2)))
if __name__=='__main__':main()
