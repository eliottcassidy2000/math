"""Independent base permutation check, without a lane engine."""
import json
from independent import count,need
rows=(158,156,8,16,36,335,159,60,223)
r=count(rows,{0,6})
need(r=={'H':143,'W':142,'U':2004,'delta':1},'three-clock base')
pairs=[(a,b) for a in range(9) for b in range(a+1,9) if a not in (0,6) and b not in (0,6) and all(bool(rows[a]>>v&1)==bool(rows[b]>>v&1) for v in range(9) if v not in (a,b))]
need(not pairs,'unselected twins')
out={'rows':rows,'selected':[0,6],'base':r,'unselected_twin_pairs':pairs}
r4=count((670,668,8,528,36,847,671,60,735,148),{0,6})
need(r4=={'H':543,'W':542,'U':8640,'delta':1},'four-clock base')
out['four_clock_base']=r4
print(json.dumps(out,sort_keys=True,indent=2))
