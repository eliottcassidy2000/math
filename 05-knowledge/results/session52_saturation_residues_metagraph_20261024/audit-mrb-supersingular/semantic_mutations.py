#!/usr/bin/env python3
from pathlib import Path
import tempfile, subprocess, sys, json
ROOT=Path(__file__).resolve().parent
source=(ROOT/'actual-manifest-copy/check.py').read_text()
mutations={
 'residue_preserving_wrong_ogg_prime':('O=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]','O=[2,3,5,7,11,13,107,19,23,29,31,41,47,59,71]'),
 'wrong_curve_with_same_claimed_invariants':('trace=-sum(legendre(x*x*x+4*x*x+2*x,p) for x in range(p))','trace=-sum(legendre(x*x*x+x,p) for x in range(p))'),
 'wrong_divisor_intersection_sign':('[(1,1),(3,-1),(5,-1),(15,1)]','[(1,1),(3,-1),(5,-1),(15,-1)]'),
 'unit_density_substitution':('unit_C=scaled[1]-scaled[3]-scaled[5]+scaled[15]','unit_C=8*scaled[1]/15'),
 'all_j_count_omits_exceptional_j_zero':('if count==p+1: js.append(j)','if count==p+1 and j!=0: js.append(j)'),
 'drop_nonfundamental_order_contribution':('h4=forms(-4*p)','h4=forms(-4*p) if p%4==1 else forms(-p)'),
}
rows=[]
with tempfile.TemporaryDirectory(prefix='mrb-semantic-audit-') as temp:
 for name,(old,new) in mutations.items():
  count=source.count(old)
  if count<1: raise RuntimeError('missing semantic mutation '+name)
  path=Path(temp)/(name+'.py'); path.write_text(source.replace(old,new))
  for opt in (False,True):
   r=subprocess.run([sys.executable]+(['-O'] if opt else [])+[str(path)],capture_output=True,text=True)
   if r.returncode==0 or 'RuntimeError:' not in r.stderr: raise RuntimeError('semantic mutation escaped '+name)
   reason=r.stderr.strip().splitlines()[-1]
   rows.append({'mutation':name,'optimized':opt,'replacements':count,'reason':reason})
   print('REJECTED',name,'optimized='+str(opt),reason)
(ROOT/'semantic_mutations.json').write_text(json.dumps(rows,indent=2)+'\n')
print('PASS: twelve semantic-mutation executions rejected')
