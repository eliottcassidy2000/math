"""Prepare optional search inputs; put rigid_51_contacts.json in the working directory."""
import sys
from replay import orbit_chart, graph, IDENTITY, SLOTS

C,S=orbit_chart('period_fourteen_site.txt',4,((10,45,2),(27,49,0)))
if sys.argv[1]=="transport":
    marks={u:1+x[1]//4 for x,E in S.items() for u in [x,*E]}
    for h in IDENTITY:
        if h not in SLOTS:print(*([-1]*168));continue
        print(*C[h,0],*C[h,28],*(marks.get((h,y),0) for y in IDENTITY))
elif sys.argv[1]=="grid":
    A=graph(C,S)
    x,y=58+56+1,58+2*56+7
    A[x].add(y);A[y].add(x)
    E=[(i-58,j-58) for i in range(58,3194) for j in sorted(A[i]) if i<j<3194]
    print(56,len(E))
    for i,j in E:print(i,j)
else:
    raise SystemExit("expected transport or grid")
