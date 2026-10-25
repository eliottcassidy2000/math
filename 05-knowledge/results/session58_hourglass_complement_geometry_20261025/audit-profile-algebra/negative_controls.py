import importlib.util,json
from pathlib import Path
p=Path(__file__).with_name('audit.py');s=importlib.util.spec_from_file_location('independent',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
out={}
for name,condition in [('false_rank34_claim',32==34),('false_single_pair_collision_claim',34<34),('false_nonzero_integral_claim',m.exact_integral()!='0')]:
 try:m.require(condition,name)
 except RuntimeError:out[name]='rejected'
 else:raise RuntimeError('Negative control accepted: '+name)
print(json.dumps(out,indent=2))
