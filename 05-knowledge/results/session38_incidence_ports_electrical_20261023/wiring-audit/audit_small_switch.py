"""Independent Bareiss / exhaustive-binary-kernel audit of n5 invisible switch."""
import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent;D=P.parent/'electrical'
if not D.exists():D=P.parent/'session38-electrical'
r=json.loads((D/'receipt.json').read_text());control=r['switches'][0]['invisible_control'];N=10

def need(q,m):
 if not q:raise ValueError(m)
def lap(blue,black):
 A=[[0]*N for _ in range(N)]
 for u,v,w in blue+[(u,v,2*w) for u,v,w in black]:
  if u!=v:A[u][u]+=w;A[v][v]+=w;A[u][v]-=w;A[v][u]-=w
 return [row[:-1] for row in A[:-1]]
def determinant(A):
 A=[r[:] for r in A];n=len(A);prev=1;sign=1
 for k in range(n-1):
  if A[k][k]==0:
   j=next(j for j in range(k+1,n) if A[j][k]);A[k],A[j]=A[j],A[k];sign=-sign
  pivot=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=pivot*A[i][j]-A[i][k]*A[k][j];need(num%prev==0,'Bareiss division');A[i][j]=num//prev
   A[i][k]=0
  prev=pivot
 return sign*A[-1][-1]
def half(A):
 n=len(A);basis=[]
 for x in itertools.product(range(2),repeat=n):
  y=[sum(a*b for a,b in zip(row,x)) for row in A]
  if all(v%2==0 for v in y):basis.append((x,y))
 q=[[sum(a*b for a,b in zip(x,Ay))//2%2 for y,Ay in basis] for x,Ax in basis]
 need(all(sum(a*b for a,b in zip(x,Ay))%2==0 for x,Ax in basis for y,Ay in basis),'half integral')
 return basis,q
A=lap(control['blue_edges'],control['original_paired_black_edges']);B=lap(control['blue_edges'],control['changed_paired_black_edges'])
a,qa=half(A);b,qb=half(B)
need([x for x,y in a]==[x for x,y in b],'binary harmonic spaces');need(qa==qb,'entire half form')
need(len(a)==8,'dimension3 kernel');need(sum(not any(row) for row in qa)==1,'half form nondegenerate')
trees=[determinant(A),determinant(B)];need(trees==[41096,36680],'independent spanning-tree determinants')
O={'status':'PASS','method':'Every binary harmonic vector (512 candidates) and all64 pairings; independent fraction-free Bareiss determinants','harmonic_dimension':3,'entire_half_form_equal':True,'half_form_nondegenerate':True,'spanning_tree_counts':trees,'complete_two_primary_both':'(Z/2)^3','strongest_scope':'A synthetic same-margin colored graph, not a second tournament realization.'}
(P/'independent_switch_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
