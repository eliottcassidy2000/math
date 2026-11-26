from nested_core_check import orbit_chart, graph, legal
import sys,json
from pathlib import Path
mode=sys.argv[1] if len(sys.argv)>1 else "marked"
C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
A=graph(C,S)
E=[] if mode=="bare" else [(2,7)] if mode=="marked" else json.loads(Path(__file__).with_name("rigid_51_contacts.json").read_text())
for h,y in E:
    w=58+56*h+y
    assert legal(A,115,w)
    A[115].add(w)
    A[w].add(115)
print(len(A))
for Z in A:print(len(Z),*sorted(Z))
