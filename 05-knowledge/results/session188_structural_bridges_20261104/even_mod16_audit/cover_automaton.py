from itertools import product
import json

def cover_kernel(a,c,reverse):
 rows=[0]*5
 for i in range(3):
  rows[i]|=1<<((i+(2 if reverse else 1))%3)
  for mask,p in [(a,3),(c,4)]:
   if mask>>i&1:rows[i]|=1<<p
   else:rows[p]|=1<<i
 rows[3]|=1<<4
 levels=[{} for _ in range(18)]
 for v in range(5):
  e=[0]*3;ports=0
  if v<3:e[v]=1
  else:ports=1<<(v-3)
  levels[1][tuple(e)+(ports,v,1)]=1
 series=[{} for _ in range(6)]
 for length in range(1,18):
  for state,z in levels[length].items():
   e=state[:3];ports,last,h=state[3:]
   if ports==3:series[h][e]=series[h].get(e,0)+z
   if length==17:continue
   for v in range(5):
    ee=list(e);pp=ports
    if v<3:
     if ee[v]==5:continue
     ee[v]+=1
    else:
     bit=1<<(v-3)
     if pp&bit:continue
     pp|=bit
    # Joining without a cut requires an actual directed quotient/port arc.
    # Cutting is always legal, including same-child labels and reversed arcs.
    for hh in ([h] if rows[last]>>v&1 else [])+([h+1] if h<5 else []):
     key=tuple(ee)+(pp,v,hh);levels[length+1][key]=levels[length+1].get(key,0)+z
 return series
