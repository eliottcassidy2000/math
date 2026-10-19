import itertools,json,sympy as s

def mat(mask,n):
 A=s.zeros(n)
 for j,(u,v) in enumerate(itertools.combinations(range(n),2)):
  if mask>>j&1:A[u,v]=1
  else:A[v,u]=1
 return A

def cycles(A,k):
 n=A.rows;total=0
 for S in itertools.combinations(range(n),k):
  a=min(S)
  for tail in itertools.permutations([x for x in S if x!=a]):
   path=(a,)+tail
   total+=all(A[path[j],path[(j+1)%k]] for j in range(k))
 return total
out=[]
for mask in [73,83]:
 A=mat(mask,6);B=s.zeros(7);B[:6,:6]=A
 for i in range(6):B[6,i]=1
 row={'mask':mask,'charpoly6':[int(x) for x in A.charpoly().all_coeffs()],'charpoly7':[int(x) for x in B.charpoly().all_coeffs()],'cycles6':cycles(A,6),'subcycles6in7':cycles(B,6),'cycles7':cycles(B,7)}
 for C in [A,B]:
  for k in [3,4,5]:
   if int(s.trace(C**k))!=k*cycles(C,k):raise RuntimeError('low order')
 out.append(row)
if out[0]['charpoly6']!=out[1]['charpoly6'] or out[0]['charpoly7']!=out[1]['charpoly7'] or [r['subcycles6in7'] for r in out]!=[0,1]:raise RuntimeError('hostile mismatch')
print(json.dumps({'status':'PASS','witnesses':out},indent=2))
