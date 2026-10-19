from pathlib import Path
from copy import deepcopy
import json
import exact_replay_engine as R
P=Path(__file__).resolve().parent;r=json.loads((P/'local_repair.json').read_text());p=(0,28,0,0,12,12,13,125,125)
def cap(x):return next(c for c in x['certificates'] if c['kind']=='ancestor_caps')
checks=[]
for name,f in [('wrong_convention',lambda c:c.update(coefficient_convention='normalized_mu')),('wrong_slack',lambda c:c['caps'][0].__setitem__(3,99)),('wrong_parent',lambda c:c['caps'][0].__setitem__(1,c['caps'][0][2])),('negative_multiplier',lambda c:c['caps'][0].__setitem__(4,-1)),('duplicate_cap',lambda c:c['caps'].append(c['caps'][0])),('false_margin',lambda c:c.update(max_margin=1))]:
 x=deepcopy(r);f(cap(x))
 try:R.replay(p,x)
 except (RuntimeError,KeyError):checks.append(name)
 else:raise RuntimeError('invalid certificate accepted '+name)
R.replay(p,r);print(json.dumps({'status':'PASS','rejected':checks,'positive_control':'passed'}))
