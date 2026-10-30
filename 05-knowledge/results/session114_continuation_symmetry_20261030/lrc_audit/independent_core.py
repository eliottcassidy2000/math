from fractions import Fraction as F
from math import floor
import json
C=tuple(range(1,12))+(13,)
checks=0
def gate(x):
 global checks
 checks+=1
 if not x:raise RuntimeError('independent gate '+str(checks))
def distance(v,t):
 x=v*t-floor(v*t)
 return min(x,1-x)
def reconstruct(r):
 cuts=sorted({F(0),F(1)}|{(k+s)/v for v in C for k in range(v) for s in (r,1-r)})
 points=[(t,t) for t in cuts if min(distance(v,t) for v in C)>=r]
 cells=[(a,b) for a,b in zip(cuts,cuts[1:]) if min(distance(v,(a+b)/2) for v in C)>=r]
 segments=[]
 for a,b in sorted(points+cells):
  if segments and a<=segments[-1][1]:segments[-1]=(segments[-1][0],max(b,segments[-1][1]))
  else:segments.append((a,b))
 return segments,len(cuts)
r=F(1,14)
segments,ncuts=reconstruct(r)
expected=[(F(15,182),F(13,154)),(F(29,70),F(41,98)),(F(57,98),F(41,70)),(F(141,154),F(167,182))]
expected+= [(F(a,14),F(a,14)) for a in (1,3,5,9,11,13)]
gate(segments==sorted(expected))
# Verify the complete affine triangles through their vertices. Every speed
# must remain on one branch; each safety inequality is affine on that branch.
triangles=[]
for a,b in segments:
 if a==b:continue
 c=min((F(k,12) for k in (1,5,7,11)),key=lambda c:abs(c-(a+b)/2))
 owners=[]
 for v in C:
  branch=floor(v*a)
  gate(floor(v*b)==branch and floor(v*c)==branch)
  for t,radius in ((a,r),(b,r),(c,F(1,12))):
   gate(radius<=v*t-branch<=1-radius)
 owners=([v for v in C if distance(v,a)==r and (v*a)%1==r],[v for v in C if distance(v,b)==r and (v*b)%1==1-r])
 gate(len(owners[0])==len(owners[1])==1)
 triangles.append(dict(left=str(a),right=str(b),center=str(c),owners=owners))
# Direct owner-based slopes, with no use of producer formula functions.
slopes=[]
for v in list(range(1,97))+[997,1009,1012,1200]:
 V=tuple(sorted(set(C+(v,))))
 if v%12:
  M=F(1,12);ts=[F(a,12) for a in (1,5,7,11)]
  base={1:13,5:5,7:7,11:11}
  predicted=2*sum((F(1,max(u,v) if v%12==res else u) for res,u in base.items()),F(0))
 elif v<=24:
  M=F(1,14);ts=[F(a,14) for a in (1,3,5,9,11,13)];predicted=F(23324,6435)
 else:
  m=v//12;M=F(m,v+5);ts=[F(5*m+2,v+5),1-F(5*m+2,v+5)];predicted=2*(F(1,5)+F(1,v))
 B=F(0)
 for t in ts:
  gate(min(distance(u,t) for u in V)==M)
  plus=[u for u in V if (u*t)%1==M];minus=[u for u in V if (u*t)%1==1-M]
  gate(bool(plus) and bool(minus))
  B+=F(1,max(plus))+F(1,max(minus))
 gate(B==predicted)
 slopes.append(dict(tail=v,slope=str(B)))
print(json.dumps(dict(status='PASS',checks=checks,core_walls=ncuts,core_segments=[[str(a),str(b)] for a,b in segments],affine_triangles=triangles,slopes=slopes),indent=2))
