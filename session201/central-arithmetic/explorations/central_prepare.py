from itertools import product
import json
I=[[int(i==j)for j in range(4)]for i in range(4)]
O=[[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]]
Z=[[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]]
J=[[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]]
Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]]
M=dict(O=O,Z=Z,J=J)
def mm(A,B):return[[sum(A[i][k]*B[k][j]for k in range(4))for j in range(4)]for i in range(4)]
W=['','J']
for le in range(2,10):W+=['J'+''.join(w)+'J'for w in product('OZJ',repeat=le-2)]
forms=[];seen=set()
for w in W:
 A=I
 for c in w:A=mm(M[c],A)
 A=mm(Q,A);key=tuple(map(tuple,A))
 if key not in seen:forms.append((w,A));seen.add(key)
json.dump(forms,open(__file__.replace('central_prepare.py','central_forms.json'),'w'))
with open(__file__.replace('central_prepare.py','central_forms.h'),'w')as f:
 f.write('const int FORMCOUNT='+str(len(forms))+';\n')
 f.write('const int length[]={'+','.join(str(len(w))for w,A in forms)+'};\n')
 f.write('const int roots[]={'+','.join(str(w.count('J'))for w,A in forms)+'};\n')
 P=[[1,0,2,1],[0,0,1,0],[1,1,2,2],[0,0,1,1]];u=[1,1,2,2]
 bounds=[]
 for w,A in forms:
  C=mm(mm(list(map(list,zip(*P))),A),P)
  bounds.append(min(4*C[i][j]//(u[i]*u[j])for i in range(4)for j in range(4)))
 f.write('const int lower4[]={'+','.join(map(str,bounds))+'};\n')
 f.write('const long long forms[][4][4]='+str([A for w,A in forms]).replace('[','{').replace(']','}')+';\n')
print(len(forms))
