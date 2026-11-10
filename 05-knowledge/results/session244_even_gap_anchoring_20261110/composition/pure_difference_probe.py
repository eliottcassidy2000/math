import json
O=((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1))
Z=((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1))
J=((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))
M=dict(O=O,Z=Z,J=J)
def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)) for i in range(4))
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
states={I:''};allrows={};reports=[]
for n in range(1,13):
 new={}
 for a,w in states.items():
  for c in M:new.setdefault(mul(M[c],a),w+c)
 states=new
 found=[]
 for a,w in states.items():
  if w[-1]!='Z':continue
  side=a[1:];key=(side,a[0]);allrows.setdefault(key,w)
 for (side,r),w in allrows.items():
  for j in range(4):
   t=tuple(x+2*(i==j) for i,x in enumerate(r))
   if (side,t) in allrows:found.append(dict(row=r,coordinate=j,low=w,high=allrows[(side,t)],side=side))
 reports.append(dict(length=n,distinct_matrices=len(states),all_terminal_rows=len(allrows),pairs=found[:30],num_pairs=len(found)))
 print(json.dumps(reports[-1]),flush=True)
 if found:break
