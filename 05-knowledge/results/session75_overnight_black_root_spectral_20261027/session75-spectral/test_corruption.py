#!/usr/bin/env python3
import copy,json,pathlib
import check
root=pathlib.Path(__file__).resolve().parent
c=json.loads((root/'certificate.json').read_text())
trials=[]
x=copy.deepcopy(c['cases'][0]);x['subset_lower_bounds'][1]=1;trials.append(('nonzero_singleton',x))
x=copy.deepcopy(c['cases'][1]);x['minimum']=2;trials.append(('wrong_target_minimum',x))
x=copy.deepcopy(c['cases'][1]);x['subset_lower_bounds'][-1]=4;trials.append(('false_lower_bound',x))
x=copy.deepcopy(c['cases'][0]);x['witness']['children'][0]['vertex']=x['witness']['children'][1]['vertex'];trials.append(('duplicate_named_leaf',x))
rejected=[]
for name,bad in trials:
 try:check.check_bellman(bad)
 except RuntimeError:rejected.append(name)
 else:raise RuntimeError('missed corruption: '+name)
print(json.dumps({'status':'PASS','rejected':rejected},sort_keys=True,indent=2))
