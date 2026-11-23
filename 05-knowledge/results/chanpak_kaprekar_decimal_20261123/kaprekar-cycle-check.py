from itertools import combinations, permutations

def height(A):
 n=len(A);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S in range(1,len(d)):
  for i,c in enumerate(d[S]):
   if not c:continue
   r=A[i]&~S
   while r:
    b=r&-r;r^=b;d[S|b][b.bit_length()-1]+=c
 return sum(d[-1])

def flip(A,i,j):
 A=list(A);A[i]^=1<<j;A[j]^=1<<i;return tuple(A)

def relabel(A,p):
 return tuple(sum((A[i]>>j&1)<<p[j] for j in range(len(A))) for i in sorted(range(len(A)),key=lambda x:p[x]))

def check_rows(A):
 assert all(not (a>>i&1) and 0<=a<1<<len(A) for i,a in enumerate(A))
 assert all((A[i]>>j&1)+(A[j]>>i&1)==1 for i,j in combinations(range(len(A)),2))

def walk(A,E,h,p):
 check_rows(A)
 B=A;N=[height(B)]
 for i,j in E:
  B=flip(B,i,j);N.append(height(B))
 assert N==h and B==relabel(A,p)
 print('order',len(A),'exact heights',N,'labelled return',B==A)

def check_classes():
 E=tuple(combinations(range(5),2))
 def mask(A):return sum((A[i]>>j&1)<<k for k,(i,j) in enumerate(E))
 def can(A):return min(mask(relabel(A,p)) for p in permutations(range(5)))
 seen=set()
 for m in range(1<<len(E)):
  A=[0]*5
  for k,(i,j) in enumerate(E):
   if m>>k&1:A[i]|=1<<j
   else:A[j]|=1<<i
  if height(A)==9:seen.add(can(A))
 assert seen=={can((6,4,8,19,7)),can((2,4,9,19,7))}
 print('exact five-vertex H=9 class check',sorted(seen))

if __name__=='__main__':
 check_classes()
 walk((486,484,488,499,487,448,384,256,0),
      ((2,8),(2,4),(7,8),(6,8),(2,4),(6,8),(2,8)),
      [9,81,115,63,37,27,45,9],(0,1,2,3,4,5,6,8,7))
 walk((250,252,249,208,224,200,128,0),
      ((0,7),(0,3),(0,6),(3,5),(0,3),(0,7),(3,5),(0,6)),
      [9,81,83,109,63,73,27,45,9],tuple(range(8)))
