from nested_core_check import graph, orbit_chart
from pathlib import Path
import json,sys
C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
A=graph(C,S)
# Mode 1 adds the marker only. Mode 2 adds the remaining 50 reference contacts.
mode=int(sys.argv[1]) if len(sys.argv)>1 else 1
if mode:
    refs=json.loads(Path("rigid_51_contacts.json").read_text())
    for h,y in refs:
        if mode==1 and (h,y)!=(2,7):continue
        x,z=115,58+56*h+y
        A[x].add(z);A[z].add(x)
N=3136
edges=[(x-58,y-58) for x in range(58,58+N) for y in A[x] if x<y<58+N]
print(56,len(edges))
for x,y in edges:print(x,y)
