from pathlib import Path
import json
D=Path(__file__).parent
def sub(m):
 s=m
 while True:
  yield s
  if not s:break
  s=(s-1)&m

def main():
 models=json.load(open(D/'models.json'));it=iter(map(int,(D/'named_values.txt').read_text().split()));out=[]
 for m in models:
  n=next(it);l=next(it);ms=[v for v in range(1<<n) if v.bit_count()<=2]
  if n!=len(m['rows']) or l!=len(ms):raise RuntimeError('named bank dimensions')
  pos={v:i for i,v in enumerate(ms)};values=[[next(it) for _ in ms] for _ in ms];coeff=[]
  for a in ms:
   row=[]
   for c in ms:row.append(sum((-1)**(a.bit_count()+c.bit_count()-u.bit_count()-v.bit_count())*values[pos[u]][pos[v]] for u in sub(a) for v in sub(c)))
   coeff.append(row)
  out.append(dict(**m,port_edge='x->y',a_convention='old vertex points to x',c_convention='old vertex points to y',physical_masks=ms,values=values,coefficients=coeff))
 if next(it,None) is not None:raise RuntimeError('trailing physical counts')
 json.dump(out,open(D/'named_banks.json','w'),separators=(',',':'))
 print('complete named response coordinates',sum(len(r['physical_masks'])**2 for r in out))
if __name__=='__main__':main()
