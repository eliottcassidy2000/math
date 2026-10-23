"""Independent finite controls for proved incidence transport and all-size hostiles."""
from itertools import combinations
from pathlib import Path
import json,argparse
HERE=Path(__file__).resolve().parent
checks=0
def need(x,m):
 global checks
 checks+=1
 if not x:raise AssertionError(m)
def incidences(A):
 n=len(A);cv=[0]*n;ce=[[0]*n for _ in range(n)];tri={}
 for i,j,k in combinations(range(n),3):
  z=int(A[i][j]==A[j][k]==A[k][i]);tri[i,j,k]=z
  for u in (i,j,k):cv[u]+=z
  for u,v in ((i,j),(i,k),(j,k)):ce[u][v]+=z;ce[v][u]+=z
 return cv,ce,tri
def E(A):
 c,e,_=incidences(A)
 return tuple(sorted(tuple(sorted((c[v],e[u][v]) for v in range(len(A)) if A[u][v])) for u in range(len(A))))
def ham_subsets(A):
 n=len(A);full=(1<<n)-1;outs=[sum(A[u][v]<<v for v in range(n)) for u in range(n)];dp=[[0]*n for _ in range(1<<n)]
 for u in range(n):dp[1<<u][u]=1
 for mask in range(1,full+1):
  for u,x in enumerate(dp[mask]):
   if x:
    options=outs[u]&~mask
    while options:
     b=options&-options;v=b.bit_length()-1;options-=b;dp[mask|b][v]+=x
 h=list(map(sum,dp));return h[full],sum(h[S]*h[full^S] for S in range(1,full))
def main():
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();O={'scope':'PROVED formulas with exhaustive finite checks n3..6, exact n15 hostiles','transport':[]};path_hostile=None
 for n in range(3,7):
  free=[(i,j) for i,j in combinations(range(n),2) if j!=i+1];groups={}
  for mask in range(1<<len(free)):
   A=[[0]*n for _ in range(n)]
   for i in range(n-1):A[i][i+1]=1
   for k,(i,j) in enumerate(free):A[i][j]=mask>>k&1;A[j][i]=1-A[i][j]
   for i in range(n-1):A[i+1][i]=0
   B=[row[:] for row in A]
   for i,j in free:B[i][j]=1-B[i][j];B[j][i]=1-B[j][i]
   cv,ce,tri=incidences(A);dv=[0]*n;de=[[0]*n for _ in range(n)];touched=0
   for (i,j,k),C in tri.items():
    adj=[(u,v) for u,v in ((i,j),(j,k)) if v==u+1]
    if not adj:delta=0
    elif len(adj)==2:delta=1-2*C
    else:
     u,v=adj[0];w=next(x for x in (i,j,k) if x not in (u,v));delta=A[u][w]-A[v][w]
    touched+=bool(adj)
    for u in (i,j,k):dv[u]+=delta
    for u,v in ((i,j),(i,k),(j,k)):de[u][v]+=delta;de[v][u]+=delta
   bv,be,bt=incidences(B);need(touched==(n-2)**2,'wrong local update count');need([cv[u]+dv[u] for u in range(n)]==bv,'vertex incidence transport');need([[ce[u][v]+de[u][v] for v in range(n)] for u in range(n)]==be,'edge incidence transport')
   for u in range(n):need(cv[u]==sum(ce[u][v] for v in range(n) if A[u][v]),'out balance')
   source,target=E(A),E(B)
   if source in groups and groups[source][0]!=target and path_hostile is None:
    old=groups[source];path_hostile={'n':n,'fixed_path':'0->1->...->n-1','free_order':free,'masks':[old[1],mask],'source_rows':[old[2],[sum(x<<v for v,x in enumerate(row)) for row in A]],'partner_signatures':[old[0],target]}
   groups[source]=(target,mask,[sum(x<<v for v,x in enumerate(row)) for row in A])
  O['transport'].append({'n':n,'fixed_path_masks':1<<len(free),'affected_triangles_per_mask':(n-2)**2})
 need(path_hostile is not None,'missing path dependence hostile');O['source_signature_does_not_supply_path']=path_hostile
 S1=[14,13,12,11,5,6,8];S2=[14,13,12,4,10,9,8];graphs=[[[int((j-i)%15 in S) for j in range(15)] for i in range(15)] for S in (S1,S2)];need(E(graphs[0])==E(graphs[1]),'circulant E mismatch');records=[]
 for S,A in zip((S1,S2),graphs):
  c,e,_=incidences(A);H,F2=ham_subsets(A);v4=[]
  for u in range(15):v4.append(sum(sorted(sum(A[x][y] for y in B) for x in B)==[1,1,2,2] for B in combinations(range(15),4) if u in B))
  need(set(c)=={28},'hostile vertex triangles');need(set(v4)=={184},'hostile rooted strong-four');need(all(sorted(e[u][v] for v in range(15) if A[u][v])==[3,3,3,4,4,5,6] for u in range(15)),'hostile edge profile');records.append({'connection_set':S,'H':H,'F2':F2,'score':7,'vertex_triangles':28,'rooted_strong4':184,'outgoing_edge_triangle_counts':[3,3,3,4,4,5,6]})
 need([r['H'] for r in records]==[196264905,196277235],'hostile path counts');need(records[0]['H']!=records[1]['H'],'nonisomorphism not witnessed');O['circulant_n15']=records;O['checks']=checks
 O=json.loads(json.dumps(O));dest=HERE/'structural_receipt.json'
 if args.write:dest.write_text(json.dumps(O,indent=2)+'\n')
 else:need(json.loads(dest.read_text())==O,'frozen structural receipt differs')
 print(json.dumps(O,indent=2))
if __name__=='__main__':main()
