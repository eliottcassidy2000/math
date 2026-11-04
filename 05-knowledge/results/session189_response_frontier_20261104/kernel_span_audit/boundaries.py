#!/usr/bin/env python3
from pathlib import Path
import json,subprocess,hashlib
from audit import physical,need,control
P=Path(__file__).resolve().parent

def main():
 cases=[('singleton_out',[0],0),('singleton_in',[0],1),('mixed_group_triangle',[2,0],2)];graphs=[]
 for _,r,i in cases:graphs.extend([physical(3,r,i),physical(3,r,i,True)])
 text=str(len(graphs))+'\n'+''.join(str(len(r))+'\n'+' '.join(map(str,r))+'\n'for r in graphs);data=subprocess.check_output([str(P/'literal_polynomial')],input=text,text=True);(P/'boundary_banks.json').write_text(data);banks=json.loads(data);records=[]
 for j,(name,r,into)in enumerate(cases):
  left,right=banks[2*j:2*j+2];L=len(left['masks']);ix={a:i for i,a in enumerate(left['masks'])};dif=[a-b for a,b in zip(left['coefficients'],right['coefficients'])];quartics=[dif[ix[3]*L+ix[(1<<2)|(1<<(3+i))]]for i in range(len(r))];nonzero=next((p,v)for p,v in enumerate(dif)if v);need(nonzero,'generic named detectability persists')
  if len(r)==1:need(quartics==[0],'m3 selected restricted battery can vanish')
  else:need(quartics==[4,-4]and sum(quartics)==0,'mixing incidence sides cancels restricted quartics')
  moments=control(r,into)['raw_moments'];records.append({'name':name,'m':3,'exterior_rows':r,'into':into,'raw_moments':moments,'restricted_quartics':quartics,'generic_nonzero_coefficient':{'a_mask':left['masks'][nonzero[0]//L],'c_mask':left['masks'][nonzero[0]%L],'difference':nonzero[1]}})
 (P/'boundary_controls.json').write_text(json.dumps(records,sort_keys=True,indent=2)+'\n');out={'status':'PASS','graphs':len(banks),'full_named_coefficients':sum(len(b['coefficients'])for b in banks),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'artifacts':{f:hashlib.sha256((P/f).read_bytes()).hexdigest()for f in ['boundary_banks.json','boundary_controls.json']}};(P/'boundary_receipt.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
