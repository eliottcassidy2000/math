from itertools import permutations,product
import json
A=[[0,0,0,0,1,1,1,0],[1,0,0,0,0,1,1,0],[1,1,0,0,0,0,1,0],[1,1,1,0,0,0,0,0],[0,1,1,1,0,0,0,0],[0,0,1,1,1,0,0,0],[0,0,0,1,1,1,0,1],[1,1,1,1,1,1,0,0]]
paths=[(3,0,6,7,1,5,4,2),(3,1,6,7,0,5,4,2),(2,0,6,7,1,5,4,3),(2,1,6,7,0,5,4,3)]
g=(4,0,2,3,1,5,6,7)
checks=0

def need(c,msg):
 global checks
 checks+=1
 if not c:raise RuntimeError(msg)
H=0;aut=[]
for p in permutations(range(8)):
 H+=all(A[u][v] for u,v in zip(p,p[1:]))
 if all(A[u][v]==A[p[u]][p[v]] for u,v in product(range(8),repeat=2)):aut.append(p)
need(H==409 and aut==[tuple(range(8))],'original actual core path count and automorphisms')

def partner(B,P):
 E={frozenset(e) for e in zip(P,P[1:])};n=len(B)
 return [[0 if i==j else B[i][j] if frozenset((i,j)) in E else B[j][i] for j in range(n)] for i in range(n)]

hosts=[]
for r,s in list(product(range(9),repeat=2))+[(32,0),(0,32)]:
 n=8+r+s;pre=tuple(range(8,8+r));post=tuple(range(8+r,n));order=pre+tuple(range(8))+post;rank={v:i for i,v in enumerate(order)}
 B=[[0 if u==v else A[u][v] if u<8 and v<8 else int(rank[u]<rank[v]) for v in range(n)] for u in range(n)]
 P=[pre+p+post for p in paths];G=g+tuple(range(8,n))
 need(all(B[u][v]+B[v][u]==1 for u in range(n) for v in range(u+1,n)),'actual padded tournament')
 deg=list(map(sum,B));need(all(deg.count(deg[u])==1 for u in range(8,n)),'new owners uniquely determined by degree')
 need(all(s<=deg[u]<=s+7 for u in range(8)),'core degree interval disjoint')
 K=[]
 for p in P:
  need(all(B[u][v] for u,v in zip(p,p[1:])),'actual Hamiltonian path')
  need(p[:r]==pre and (not s or p[-s:]==post),'literal prefix and suffix')
  need(p[r+2:r+4]==(6,7) and p[r+5:r+7]==(5,4),'shifted fixed shared arcs')
  K.append(partner(B,p))
  C=list(map(list,zip(*B)));R=partner(C,p[::-1])
  need(R==list(map(list,zip(*K[-1]))),'literal converse partner square')
 for i,j in [(0,1),(2,3)]:
  need(all(K[i][u][v]==K[j][G[u]][G[v]] for u,v in product(range(n),repeat=2)),'extended actual target isomorphism')
  need(G[P[i][r]]==P[j][r] and G[P[i][r+7]]==P[j][r+7],'core endpoints fixed by target transporter')
  need([a for a in range(n) if P[i][a]!=P[j][a]]==[r+1,r+4],'constant literal support')
 compatible=[(i,j) for i in (0,1) for j in (2,3) if P[i][r+1:r+7]==P[j][r+1:r+7]]
 need(compatible==[(0,2),(1,3)],'two retained states; no one-slot edge')
 need(len(set(P))==4 and len({tuple(map(tuple,z)) for z in K})==4,'four distinct path and target-matrix resources')
 hosts.append({'r':r,'s':s,'n':n,'slot_budget':2})
print(json.dumps({'status':'PASS','checks':checks,'core_hamiltonian_paths':H,'core_automorphisms':len(aut),'hosts':hosts,'scope':'Actual family at every n>=8 by proof; finite direct tests83 hosts; numeric higher-order potentials not assumed.'},indent=2))
