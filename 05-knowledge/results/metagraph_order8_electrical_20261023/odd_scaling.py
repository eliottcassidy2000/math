"""Exact positive-odd black-weight scaling test on the full order-eight graph.
Imports the audited filtration engine; no changes to that engine or its inputs.
Default compares the compact frozen sidecar; --write regenerates it.
"""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('dyadic',HERE/'check.py');dyadic=importlib.util.module_from_spec(spec);spec.loader.exec_module(dyadic)
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
 folder=HERE.parent/'metagraph_order8_20261023';nr=dyadic.raw_read(folder/'merged_nodes.tsv');er=dyadic.raw_read(folder/'merged_edges.tsv');n=len(list(csv.DictReader(nr.decode().splitlines(),delimiter='\t')));L=[];B=[]
 for r in csv.DictReader(er.decode().splitlines(),delimiter='\t'):
  i,j,b,k=(int(r[x]) for x in ('source','target','blue','black'))
  dyadic.need(k%2==0,'paired black')
  if b:L.append((i,j,b))
  if k:B.append((i,j,k//2))
 base=json.loads(dyadic.raw_read(HERE/'receipt.json'))['order8'];cases=[]
 for scale in (1,3,5,7):
  before=dyadic.CHECKS
  out=base if scale==1 else dyadic.calculate(n,L,[(i,j,scale*w) for i,j,w in B])
  multiplicities={1:out['exact_valuation_one'],2:out['exact_valuation_two'],3:out['exact_valuation_three']}
  layers=[{'level':1,'dimension':out['critical_two_rank'],'rank':out['half_form_rank'],'next_dimension':out['four_divisible_factors']},{'level':2,'dimension':out['four_divisible_factors'],'rank':out['quarter_form_rank'],'next_dimension':out['eight_divisible_factors']},{'level':3,'dimension':out['eight_divisible_factors'],'rank':out['exact_valuation_three'],'next_dimension':out['sixteen_divisible_factors']}]
  for row in out['higher_layers']:
   multiplicities[row['level']]=row['rank'];layers.append({k:row[k] for k in ('level','dimension','rank','next_dimension')})
  exact=[(e,m) for e,m in sorted(multiplicities.items()) if m]
  dyadic.need(out['complete_two_primary'],'finite filtration did not complete');dyadic.need(sum(multiplicities.values())==out['critical_two_rank'],'invariant factor multiplicity total')
  dyadic.need(out['critical_two_rank']==base['critical_two_rank'] and out['half_form_rank']==base['half_form_rank'],'odd scaling changed first layer')
  case={'black_scale':scale,'complete':out['complete_two_primary'],'layers':layers,'two_primary_elementary_divisors':[[2**e,m] for e,m in exact],'determinant_v2':sum(e*m for e,m in exact),'computed_checks':dyadic.CHECKS-before,'full_calculation_sha256':hashlib.sha256(json.dumps(out,sort_keys=True,separators=(',',':')).encode()).hexdigest()};cases.append(case);print(json.dumps(case),flush=True)
 result={'status':'FINITE-EXACT odd scaling b=1,3,5,7 at n=8; no all-odd scaling assertion','definition':'G_b = blue Laplacian + 2*b*paired-black Laplacian','input_sha256_lf':{'nodes':hashlib.sha256(nr).hexdigest(),'edges':hashlib.sha256(er).hexdigest(),'engine':hashlib.sha256(dyadic.raw_read(HERE/'check.py')).hexdigest(),'base_receipt':hashlib.sha256(dyadic.raw_read(HERE/'receipt.json')).hexdigest()},'cases':cases,'checks':dyadic.CHECKS}
 raw=(json.dumps(result,indent=2)+'\n').encode()
 if args.write:(HERE/'odd_scaling.json').write_bytes(raw)
 else:dyadic.need(raw==dyadic.raw_read(HERE/'odd_scaling.json'),'frozen scaling sidecar mismatch')
 print('sidecar_sha256',hashlib.sha256(raw).hexdigest())
if __name__=='__main__':main()
