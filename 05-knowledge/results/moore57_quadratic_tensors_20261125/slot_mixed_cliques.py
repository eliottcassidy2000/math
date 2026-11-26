from slot_interface_check import read_board
import sys
from itertools import combinations
from collections import Counter,defaultdict
_,R=read_board(sys.argv[1]);m=R[0][7];C=[]
L=defaultdict(int)
for i,D in enumerate(R[1:],1):
 for y,t in D.items():
  if any(u==v for u,v in zip(t,m)):continue
  v=len(C);e=tuple(enumerate((i,*t)))
  C.append((i,y,t,e))
  for c in e:L[c]|=1<<v
def each(s):
 while s:
  b=s&-s;s-=b;yield b.bit_length()-1
N=[]
for i,(_,_,_,e) in enumerate(C):
 z=0
 for c in e:z|=L[c]
 N.append(z^(1<<i))
maximal=set()
def grow(r,p,x):
 if not p:
  if not x:maximal.add(tuple(sorted(r)))
  return
 for w in each(p):
  grow(r+(w,),p&N[w],x&N[w]);p^=1<<w;x|=1<<w
for i,(_,_,_,e) in enumerate(C):
 for c in e:
  for j in each(L[c]>> (i+1)):
   j+=i+1
   for k in each(((N[i]&N[j])&~L[c])>>(j+1)):
    k+=j+1
    # All maximal extensions of this triangle (only forward among its added sites).
    grow((i,j,k),N[i]&N[j]&N[k],0)
T=Counter(map(len,maximal))
print(T)
K=max(sorted(maximal),key=len);print([C[i][:3] for i in K])
print("original identities",K)
assert T=={3:6948,4:88248,5:822}
assert K==(0,8,19,313,1438)
