from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import sys,json
from explore import *
import literal
ns={'state':literal.state}
a,_,names=family(0,0);got=states(0,0)
for i,name in enumerate(names):
 want=ns['state'](a,i)
 if got[name]!=want:raise ValueError((name,got[name],want))
print(json.dumps({'status':'PASS','roots':got},indent=2))
