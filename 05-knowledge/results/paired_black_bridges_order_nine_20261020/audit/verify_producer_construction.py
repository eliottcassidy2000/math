"""Read-only re-execution of the producer constructor, with independent checks."""
import importlib.util,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent

def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
ind=load('independent',ROOT/'audit/verify_score_construction.py')
producer=load('producer',ROOT/'construct_black_quadrilaterals.py')
claims=json.loads((ROOT/'black_quadrilateral_construction_checks.json').read_text())
assert claims['checked_orders']==len(claims['rows'])==101
for row in claims['rows']:
 n=row['n'];a,p,q=producer.construct(n)
 binary=[[bool(a[i]>>j&1) for j in range(n)] for i in range(n)]
 result=ind.verify(binary,p,q)
 assert result['alpha']==row['a'] and result['beta']==row['b'] and result['four_c3_values']==row['triangle_counts']
(ROOT/'audit/producer_construction_audit.json').write_text(json.dumps({'checked_orders':len(claims['rows']),'orders':[r['n'] for r in claims['rows']],'all_checks':'passed'},indent=2)+'\n')
print('Independent checks of all 101 producer construction outputs passed.')
