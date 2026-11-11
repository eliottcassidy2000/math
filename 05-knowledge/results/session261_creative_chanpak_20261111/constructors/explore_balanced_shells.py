from module_tools import *
import json

def shell(T):
 n=len(T);R=tuple(sum(1<<i for i in range(n) if T[i]>>j&1) for j in range(n));out=[]
 for i in range(n):out.append(T[i]|(((1<<n)-1)<<n))
 for i in range(n):out.append((R[i]<<n)|(1<<(2*n)))
 out.append((1<<n)-1)
 return tuple(out)
out=[]
for n in range(1,5):
 cs={canon(tour(n,b)) for b in range(1<<(n*(n-1)//2))}
 for T in sorted(cs):
  Z=shell(T)
  for O in range(1,1<<n):
   for apex in (0,1):
    P=O|(O<<n)|(apex<<(2*n));h,w,u=m.data(Z,P)
    if h-w==1:
     d={'core':T,'selected_pair_mask':O,'apex_selected':apex,'H':h,'U':u,'rows':Z,'O':P};out.append(d);print(d,flush=True)
Path(__file__).with_name('balanced_shells.json').write_text(json.dumps(out,indent=2))
