"""Exploratory finite pure-D matrix search. No all-word completeness claim."""
import json,sys
O=lambda a:(a[0],(0,0,0,0),tuple(2*x for x in a[2]),a[3])
Z=lambda a:(tuple(x-y+z for x,y,z in zip(a[0],a[1],a[2])),a[3],tuple(2*x-y for x,y in zip(a[2],a[3])),a[3])
J=lambda a:(tuple(x-y+z for x,y,z in zip(a[0],a[1],a[2])),tuple(x-y for x,y in zip(a[0],a[1])),tuple(x+y for x,y in zip(a[0],a[2])),a[0])
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
state={I:''};rows={}
for n in range(1,int(sys.argv[1] if len(sys.argv)>1 else 13)+1):
 nxt={}
 for a,w in state.items():
  for c,F in (('O',O),('Z',Z),('J',J)):
   b=F(a);nxt.setdefault(b,w+c)
   # Terminal admissibility belongs to a witness, not the chosen matrix representative.
   if c=='Z':rows.setdefault((b[1:],b[0]),w+c)
 state=nxt
 counts=[0]*4;escapes=[]
 for (side,r),w in rows.items():
  for j in range(4):
   t=tuple(x+2*(i==j) for i,x in enumerate(r))
   if (side,t) in rows:
    counts[j]+=1
    if any(x%3 for x in side[-1]):escapes.append(dict(low=w,high=rows[(side,t)],coordinate=j,matrix=[r]+list(side)))
 print(json.dumps(dict(length=n,matrices=len(state),pairs_by_coordinate=counts,divisor_escapes=escapes[:20],escape_count=len(escapes))),flush=True)
 if escapes:break
