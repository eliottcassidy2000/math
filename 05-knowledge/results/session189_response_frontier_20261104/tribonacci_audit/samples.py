#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,hashlib
from controls import need,valid,birth,nearjoin
P=Path(__file__).resolve().parent

def decode(n,code,op):
 r=[0]*n;b=n*(n-1)//2-1
 for j in range(1,n):
  for i in range(j):
   a=((code>>b)&1)^op;b-=1
   if a:r[i]|=1<<j
   else:r[j]|=1<<i
 return r

def main():
 packed=[int(x,16)for x in(P/'witnesses.hex').read_text().splitlines()];targets=[69,12875,12877,38629,38635,702189];records=[];requests=[]
 for v in targets:
  z=packed[(v-69)//2];code=z&((1<<36)-1);t=z>>36&15;n=z>>40&15;op=z>>44&1;mode=z>>45&3;k=z>>47&15;r=decode(n,code,op);need(valid(r),'physical witness base');out=r[:];root=t
  if mode==1:out=nearjoin(r,t,k)
  else:
   for _ in range(k):out=birth(out,root);root=len(out)-1
  need(valid(out),'physical output tournament');records.append({'H':v,'base_order':n,'base_rows':r,'root':t,'mode':mode,'height':k,'output_order':len(out),'output_rows':out});requests.append((out,0))
 text=''.join(str(len(r))+' '+str(t)+' '+' '.join(map(str,r))+'\n'for r,t in requests);(P/'sample_graphs.in').write_text(text);data=subprocess.check_output([str(P/'ham')],input=text,text=True);(P/'sample_counts.txt').write_text(data);counts=[list(map(int,l.split()))for l in data.splitlines()]
 for record,count in zip(records,counts):need(count[0]==record['H']and count[3]==1,'direct full physical H and strongness');record['direct_H']=count[0];record['strong']=True
 (P/'sample_controls.json').write_text(json.dumps(records,sort_keys=True,indent=2)+'\n');receipt={'status':'PASS','direct_physical_values':targets,'maximum_sample_order':max(r['output_order']for r in records),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['sample_graphs.in','sample_counts.txt','sample_controls.json']}};(P/'sample_receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
