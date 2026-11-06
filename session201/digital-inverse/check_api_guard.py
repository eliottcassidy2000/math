"""Release-only check of the original-index API guard; internal h=1 is legal."""
import json
from verify import fixed_inverse

try:
    fixed_inverse((1,1,0,0),2,vmin=1)
except RuntimeError as error:
    if str(error)!='original index lower bound must be at least two':raise
else:raise RuntimeError('out-of-domain original singleton was accepted')

hit,stats=fixed_inverse((1,1,0,0),3,vmin=2)
if hit!=[2]:raise RuntimeError('valid source or internal scalar singleton was rejected')
print(json.dumps({'status':'PASS','vmin_1':'rejected','valid_v_2':hit,'statistics':stats},sort_keys=True,indent=2))
