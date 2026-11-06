import json
F=json.load(open(__file__.replace('digit_descent.py','forms.json')))
I=[[int(i==j)for j in range(4)]for i in range(4)]
B=[[[1,0,0,0],[0,2,0,0],[0,0,0,0],[0,0,0,1]],[[1,1,-1,0],[0,2,0,-1],[0,0,0,1],[0,0,0,1]]]
def tr(A):return list(map(list,zip(*A)))
def mm(A,B):return [[sum(A[i][k]*B[k][j]for k in range(4))for j in range(4)]for i in range(4)]
def add(A,B,scale=1):return [[A[i][j]+scale*B[i][j]for j in range(4)]for i in range(4)]
for side in ['left','right','both']:
 hits=[]
 for name,M in F.items():
  for i in range(2):
   for j in range(2 if side=='both'else 1):
    A=mm(tr(B[i]),M)if side!='right'else mm(M,B[i])
    if side=='both':A=mm(A,B[j])
    for target,N in F.items():
     for scale in [2,4]:
      D=add(A,N,-scale)
      if all(D[a][b]==0 for a in range(4)for b in range(4)if(a,b)!=(3,3)):
       hits.append((name,i,j,scale,target,D[3][3]))
 print(side,len(hits),hits[:100])
def proj(M,p,q):
 A=[[1,0,0],[0,1,0],[0,0,p],[0,0,1]];B0=[[1,0,0],[0,1,0],[0,0,q],[0,0,1]]
 return [[sum(A[k][i]*M[k][l]*B0[l][j]for k in range(4)for l in range(4))for j in range(3)]for i in range(3)]
for side in ['left','right','both']:
 hits=[]
 for name,M in F.items():
  for i in range(2):
   for j in range(2 if side=='both'else 1):
    A=mm(tr(B[i]),M)if side!='right'else mm(M,B[i])
    if side=='both':A=mm(A,B[j])
    for p in range(2):
     for q in range(2):
      AA=proj(A,p,q)
      for target,N in F.items():
       NN=proj(N,p,q)
       for scale in [2,4]:
        D=[[AA[a][b]-scale*NN[a][b]for b in range(3)]for a in range(3)]
        if all(D[a][b]==0 for a in range(3)for b in range(3)if(a,b)!=(2,2)):
         hits.append((name,i,j,p,q,scale,target,D[2][2]))
 print('parity',side,len(hits),hits[:100])
