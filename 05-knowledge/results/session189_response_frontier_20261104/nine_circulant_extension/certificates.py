from pathlib import Path
import json
D=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def certificate(bank,safe):
 piv={};contradictions=[];basis=[]
 for i,q in enumerate(bank):
  if safe and not q['cycle_safe']:continue
  row=sum(b<<j for j,b in enumerate(q['linear_bits']));rhs=q['baseline']//4;used={i}
  while row:
   k=row.bit_length()-1
   if k not in piv:piv[k]=(row,rhs,used);basis.append(i);break
   a,b,ids=piv[k];row^=a;rhs^=b;used^=ids
  if not row and rhs:contradictions.append(sorted(used))
 out=dict(rank=len(piv),independent_query_indices=basis,independent_actual_queries=[bank[i]['query'] for i in basis])
 if contradictions:
  ids=min(contradictions,key=lambda ids:(len(ids),ids));out.update(contradiction_query_indices=ids,contradiction_actual_queries=[bank[i]['query'] for i in ids])
 else:
  need(len(piv)==18,'expected unique solution');x=0
  for k,(row,rhs,_) in sorted(piv.items()):x|=(rhs^((row&x).bit_count()%2))<<k
  out['unique_solution']=x
 return out

def main():
 data=json.loads((D/'analysis.json').read_text());out=[]
 for q in data:
  c=certificate(q['queries'],True)
  if 'contradiction_query_indices' in c:
   row=rhs=0
   for i in c['contradiction_query_indices']:
    r=q['queries'][i];need(r['cycle_safe'],'safe actual queries in contradiction');row^=sum(b<<j for j,b in enumerate(r['linear_bits']));rhs^=r['baseline']//4
   need(row==0 and rhs==1,'explicit affine contradiction')
  out.append(dict(rows=q['rows'],safe_certificate=c))
 need(out[0]['safe_certificate'].get('unique_solution')==511,'carousel unique E9 type');need('contradiction_query_indices' in out[1]['safe_certificate'],'C1 no type');(D/'linear_certificates.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(status='PASS',carousel_safe_basis_queries=18,noncarousel_safe_contradiction_queries=len(out[1]['safe_certificate']['contradiction_query_indices']))))
if __name__=='__main__':main()
