"""Same guarded section-height carrier in three Q-distinct K3 surfaces."""
import json,sympy as S
from named_family import family
T=S.symbols('T');out=[];target=S.Matrix([[4,-2,0],[-2,4,2],[0,2,4]])
for i in [273,302,4239]:
 _,pts=family(i,T);G=4*S.eye(3);intersections=[]
 for j in range(3):
  for k in range(j):
   dx=S.Poly(pts[j][0]-pts[k][0],T);dy=S.Poly(pts[j][1]-pts[k][1],T);finite=S.gcd(dx,dy).degree();infinity=min(4-dx.degree(),6-dy.degree());I=finite+infinity;G[j,k]=G[k,j]=2-I
   intersections.append({'sections':[k,j],'finite':int(finite),'infinity':int(infinity)})
 if G!=target:raise ValueError('common carrier failed')
 U=S.diag(1,1,-1);A=U.T*G*U
 if A!=2*S.Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]]):raise ValueError('A3 identification failed')
 out.append({'curve':i,'height_gram':[list(map(int,G.row(k)))for k in range(3)],'determinant':int(G.det()),'intersections':intersections})
print(json.dumps(out,indent=2))
