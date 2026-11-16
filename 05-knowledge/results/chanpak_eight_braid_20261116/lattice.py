from check_rank_fan import *
from itertools import product
from sympy.matrices.normalforms import hermite_normal_form
import json
def mul(A0,B0):
 need(int(A0.sum(axis=1).max())*int(B0.max())<2**63)
 return A0@B0
P=R0.copy();q=np.array(v0,dtype=np.int64)
for a in V[1:]:P=mul(A[a],P);q=mul(A[a],q)
K=Matrix(mul(L,P))
_,j=K.to_DM().rref();C=K[:,list(j)];P0=P[:,list(j)]
_,i=C.T.to_DM().rref()
den=C[list(i),:].inv();SS=O*Matrix(P0)
def app(U,w):
 for a in w:U=mul(A[a],U)
 return U
T=[];words=[]
C0=C[list(i),:]
for w in product(range(8),repeat=3):
 R=w+S
 # Exact products below are tiny (max 8 more letters, at most profile envelope ~4^8 times initial).
 Y=Matrix(mul(L,app(P0,R)))
 Z=Y[list(i),:]
 if Z.det()==0:continue
 T0=Z*den
 if Y!=C*den*Z:continue
 T.append(T0);words.append(w)
print('number exact return blocks',len(T),'pivot rows',i,flush=True)
H0=hermite_normal_form(K[list(i),:])
cnt=0
while True:
 out=[a*H0 for a in T]
 Z=Matrix.hstack(H0,*out);need(all(x.q==1 for x in Z));H1=hermite_normal_form(Matrix(Z.rows,Z.cols,list(map(int,Z))))
 cnt+=1
 if H1==H0:break
 H0=H1
print('HNF',H0,'iterations',cnt,flush=True)
B1=[H0.inv()*a*H0 for a in T];z1=H0.inv()*Matrix(mul(L,q))[list(i),:]
f1=C[0,:]*den*H0
need(all(x.q==1 for x in z1) and all(x.q==1 for x in f1) and all(x.q==1 for M1 in B1 for x in M1))
z=list(map(int,z1));f=list(map(int,f1));ms=[list(map(int,b)) for b in B1]
need(len(T)==168 and cnt==1 and z==[-43,-5,4,53] and f==[32,14,10,27])
print('starting coordinates',z,'HCOV',f,flush=True)
with open('return_automaton.json','w') as fp:json.dump({'words':words,'start':z,'readout':f,'matrices':ms},fp)
