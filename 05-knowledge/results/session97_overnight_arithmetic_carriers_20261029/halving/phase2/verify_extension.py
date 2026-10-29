"""Independent native-model rational verifier. Standard library only."""
import copy,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
D=Path(__file__).resolve().parent;checks=0;special={'O':0,'T':0}
def need(b):
 global checks
 checks+=1
 if not b:raise ValueError('extension rejected')
def rank(masks):
 rows=[[m>>i&1 for i in range(30)]for m in masks];rank=0
 for j in range(30):
  k=next((k for k in range(rank,len(rows))if rows[k][j]),None)
  if k is None:continue
  rows[rank],rows[k]=rows[k],rows[rank]
  for l in range(len(rows)):
   if l!=rank and rows[l][j]:rows[l]=[a^b for a,b in zip(rows[l],rows[rank])]
  rank+=1
 return rank

def verify(curve,z):
 need(z['curve_id']==curve['id']);a1,a2,a3,a4,a6=map(F,curve['ainvs']);b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
 delta=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6;need(delta==F(curve['discriminant']))
 def on(P):
  if P is None:return True
  x,y=P;return y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
 def add(P,R):
  need(on(P)and on(R))
  if P is None:return R
  if R is None:return P
  x,y=P;u,v=R
  if x==u and y+v+a1*x+a3==0:return None
  if P==R:m=(3*x*x+2*a2*x+a4-a1*y)/(2*y+a1*x+a3)
  else:need(x!=u);m=(v-y)/(u-x)
  nu=y-m*x;xx=m*m+a1*m-a2-x-u;yy=-(m+a1)*xx-nu-a3
  need(on((xx,yy)));return xx,yy
 def native(q):
  x,Y=map(F,q);return x,Y-(a1*x+a3)/2
 def rowscheck(rows,pts):
  masks=[]
  for row in rows:
   p,r=row['p'],row['root'];need(p>2 and all(p%d for d in range(2,math.isqrt(p)+1)))
   red=lambda v:v.numerator*pow(v.denominator,-1,p)%p
   need(red(delta)!=0);need((4*r**3+red(b2)*r*r+2*red(b4)*r+red(b6))%p==0)
   squares={x*x%p for x in range(1,p)};mask=0
   for i,(x,y) in enumerate(pts):
    if x.denominator%p==0:
     # On an integral smooth model, negative v_p(x) is reduction to O.
     need(y.denominator%p==0);bit=0;special['O']+=1
    else:
     need(y.denominator%p!=0);v=(red(x)-r)%p
     if not v:
      need((2*red(y)+red(a1)*r+red(a3))%p==0)
      v=(3*r*r+red(b2)*pow(2,-1,p)*r+red(b4)*pow(2,-1,p))%p;special['T']+=1
     need(v!=0);bit=int(v not in squares)
    mask|=bit<<i
   need(mask==row['bits']);masks.append(mask)
  return rank(masks)
 pts=[tuple(map(F,p))for p in curve['points']];M=[[F(i==j)for j in range(30)]for i in range(30)]
 for t,step in enumerate(z['steps']):
  need(step['step']==t+1);need(rowscheck(step['rows_before'],pts)==step['rank_before'])
  support=step['support'];replace=step['replace_index'];need(replace in support and len(set(support))==len(support)and all(0<=i<30 for i in support))
  P=None
  for i in support:P=add(P,pts[i])
  need(P==native(step['sum_completed_square']));H=native(step['half_completed_square']);need(on(H));need(add(H,H)==P)
  row=[sum(M[i][j]for i in support)/2 for j in range(30)]
  need(row==list(map(F,step['new_basis_row'])))
  pts[replace]=H;M[replace]=row
 need(pts==[native(p)for p in z['final_points_completed_square']]);need(M==[list(map(F,r))for r in z['final_basis_matrix']])
 finalrank=rowscheck(z['final_character_rows'],pts);need(finalrank==z['final_character_rank']==30);need(len(z['steps'])<=20)
 p=z['two_torsion_exclusion_prime'];need(p>2 and all(p%d for d in range(2,math.isqrt(p)+1)));need(delta.numerator%p!=0)
 need(all((4*x**3+int(b2)*x*x+2*int(b4)*x+int(b6))%p for x in range(p)))
 # Exact determinant of the transported basis, independently from step count.
 mat=copy.deepcopy(M);det=F(1)
 for j in range(30):
  k=next((k for k in range(j,30)if mat[k][j]),None);need(k is not None)
  if k!=j:mat[k],mat[j]=mat[j],mat[k];det=-det
  pivot=mat[j][j];det*=pivot
  for k in range(j+1,30):
   scale=mat[k][j]/pivot;mat[k]=[a-scale*b for a,b in zip(mat[k],mat[j])]
 need(det==F(1,2**len(z['steps'])));need(all(v.denominator in (1,2)for row in M for v in row))
 transport=json.loads((D/f'inverse_transport_{curve["id"]}.json').read_text());need(transport['curve_id']==curve['id'])
 N=[list(map(F,r))for r in transport['original_basis_in_final_basis']];need(all(v.denominator==1 for row in N for v in row))
 for i in range(30):
  for j in range(30):need(sum(N[i][k]*M[k][j]for k in range(30))==F(i==j))
 # Each legal half replacement also explicitly recovers the previous replaced
 # point as 2H minus the other supported points, so old L is contained in new L.
 return {'curve':curve['id'],'steps':len(z['steps']),'index':2**len(z['steps']),'binary_rank':finalrank,'quotient':f'(Z/2Z)^{len(z["steps"])}','two_saturated':True,'max_final_prime':max(r['p']for r in z['final_character_rows'])}
source=D.parent/'halving_certificate.json';packet=json.loads(source.read_text());out=[];rejected=0
for entry in packet['curves']:
 curve=entry['curve'];z=json.loads((D/f'extension_{curve["id"]}.json').read_text());need(z['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest());out.append(verify(curve,z))
 for mode in ('half','row','matrix'):
  bad=copy.deepcopy(z)
  if mode=='half':bad['steps'][0]['half_completed_square'][1]=str(F(bad['steps'][0]['half_completed_square'][1])+1)
  if mode=='row':bad['final_character_rows'][0]['bits']^=1
  if mode=='matrix':bad['final_basis_matrix'][0][0]=str(F(bad['final_basis_matrix'][0][0])+1)
  try:verify(curve,bad)
  except ValueError:rejected+=1
  else:raise ValueError('hostile accepted')
need(rejected==6)
print(json.dumps({'results':out,'checks':checks,'hostiles_rejected':rejected,'special_character_evaluations_including_hostile_replays':special,'scope':'2-saturation of these subgroups only, not global rank or odd-primary saturation'},sort_keys=True))
