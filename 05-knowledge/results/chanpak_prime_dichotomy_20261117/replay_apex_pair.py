from itertools import product
from slot_profiles import end_profiles

def test(L,w):
 R=tuple(sum(1<<j for j in range(i+1,L))|((1-w[i])<<L) for i in range(L))+(sum(v<<i for i,v in enumerate(w)),)
 h,ee=end_profiles(R); full=len(h)-1
 dual=tuple(sum((R[j]>>i&1)<<j for j in range(L+1)) for i in range(L+1))
 _,ss=end_profiles(dual)
 E=ee[full]
 for j in range(L):
  U=sum(ss[s][j]*ee[full^s][L]+ss[s][L]*ee[full^s][j] for s in range(full+1))
  found=(U+E[j]+E[L]==h[full])
  goal=(j==L-2 and all(w[:-1]) and L>=3)
  if found!=goal:raise ValueError((L,w,j))
  if 0<j<L-1:
   f=sum(w[i]<<(j-i-1) for i in range(j))
   z=sum((1-w[i])<<(i-j-1) for i in range(j+1,L))
   if E[j]!=f+w[j] or U!=z+1-w[j-1]+f+w[j+1]:raise ValueError((w,j,U,E))
for L in range(2,9):
 for b in product(range(2),repeat=L-2):test(L,(1,*b,0))
 print('all',L)
