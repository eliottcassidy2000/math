from itertools import combinations, permutations
def H(A):
 if not A:return 1
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,d in enumerate(row):
   if not d:continue
   todo=A[i]&~S
   while todo:
    z=todo&-todo;todo^=z;D[S|z][z.bit_length()-1]+=d
 return sum(D[-1])

def induced(T,rem):
 return [sum(((T[i]>>j)&1)<<k for k,j in enumerate(rem)) for i in rem]
T=[120,113,115,54,32,64,24]
v=6
if (H(induced(T,list(range(v)))),H(T))!=(5,35):raise ValueError("incorrect literal lift")
data={}
for l in (3,5,7):
 vals=[]
 for s in permutations(range(v),l-1):
  cyc=(v,)+s
  if all(T[i]>>j&1 for i,j in zip(cyc,cyc[1:]+cyc[:1])):
   vals.append(H(induced(T,[u for u in range(v) if u not in s])))
 data[l]=(len(vals),sum(vals))
if data!={3:(4,8),5:(6,6),7:(1,1)}:raise ValueError(data)
print("Exact path counts 5 -> 35; apex cycle length: (number, complement-weighted sum)",data)
for n in range(1,7):
 edges=list(combinations(range(n),2));hits=0
 for code in range(1<<len(edges)):
  A=[0]*n
  for k,(i,j) in enumerate(edges):
   if not code>>k&1:i,j=j,i
   A[i]|=1<<j
  if H(A)==35:
   hits+=sum(H(induced(A,[i for i in range(n) if i!=a]))==5 for a in range(n))
 if hits:raise ValueError((n,hits))
 print("No one-vertex 5 -> 35 lift at output order",n)
