#!/usr/bin/env python3
from fractions import Fraction as F
import json
from pathlib import Path
mul=lambda A,B:tuple(tuple(sum(A[i][k]*B[k][j]for k in range(2))for j in range(2))for i in range(2))
sub=lambda A,B:tuple(tuple(A[i][j]-B[i][j]for j in range(2))for i in range(2))
add=lambda A,B:tuple(tuple(A[i][j]+B[i][j]for j in range(2))for i in range(2))
scale=lambda c,A:tuple(tuple(c*x for x in row)for row in A)
tr=lambda A:tuple(zip(*A))
comm=lambda A,B:sub(mul(A,B),mul(B,A))
e=((0,1),(0,0));f=((0,0),(1,0));h=((1,0),(0,-1));J=((0,1),(-1,0));zero=((0,0),(0,0))
assert comm(h,e)==scale(2,e)
assert comm(h,f)==scale(-2,f)
assert comm(e,f)==h
for X in (e,f,h):assert add(mul(tr(X),J),mul(J,X))==zero
assert mul(e,e)==zero and e[0][0]==e[1][0]==0
# Bracket profile: multiplicity weights +/-1 and adjoint weights +/-2,0.
assert sum((133,2*56,2))==247 and sum((133,2*56,3))==248
assert 2*(133+56+1)==380
# Centralizer equations for X=[[a,b],[c,-a]],Y=[[u,v],[w,-u]]:
# [X,Y]=[[bw-cv,2(av-bu)],[2(cu-aw),cv-bw]].
for a,b,c in ((1,0,0),(0,1,0),(0,0,1),(1,2,3)):
 X=((a,b),(c,-a))
 for u,v,w in ((0,0,0),(a,b,c),(2*a,2*b,2*c),(1,1,1)):
  Y=((u,v),(w,-u))
  expected=((b*w-c*v,2*(a*v-b*u)),(2*(c*u-a*w),c*v-b*w))
  assert comm(X,Y)==expected
result={'paired_lower_bound':247,'proved_minimum':248,'isotropic_restriction_bound':380,'centralizer_dimensions_at_candidates':[2,3],'multiplicity_module_dimension':2,'sl2_relations_verified':True,'symplectic_multiplicity_action_verified':True,'automatic_grading_dimensions':{'-2':1,'-1':56,'0':134,'1':56,'2':1},'same_H57_symplectic_control':{'levi':28*57,'semidirect':28*57+57,'completion':29*59},'scope':'Exact matrix/dimension controls only; abstract proof in UNGRADED_METRIC_COMPLETION.md.'}
Path(__file__).with_name('ungraded_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
