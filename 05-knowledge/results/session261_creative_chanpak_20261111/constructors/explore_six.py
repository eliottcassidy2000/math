from module_tools import *
from math import gcd
import json
seen5={canon(tour(5,b)) for b in range(1024)};seen=set()
for A in seen5:
 for b in range(32):
  T=tuple(A[i]|((1<<5) if b>>i&1 else 0) for i in range(5))+(sum(1<<i for i in range(5) if not b>>i&1),)
  seen.add(canon(T))
D={}
for T in sorted(seen):D.setdefault(covers(T),[]).append(T)
out=[]
for q,ls in D.items():
 if gcd(*q)>1 or len(ls)<2:continue
 rs=[sorted(rootprofile(T,v) for v in range(6)) for T in ls]
 if len({str(x) for x in rs})>1:
  z={'Q':q,'modules':ls,'root_profiles':rs};out.append(z);print(z,flush=True)
Path(__file__).with_name('six_module_collisions.json').write_text(json.dumps({'classes':len(seen),'profiles':len(D),'collisions':out},indent=2))
print('END',len(seen),len(D),len(out))
