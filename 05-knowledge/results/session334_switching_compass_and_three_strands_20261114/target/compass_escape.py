from cert import *
M=288
F=('1','1',('0',24),('2',24),'3','1','1')
def comp(t):
 if isinstance(t,tuple):return (str(3-int(t[0])),t[1])
 return str(3-int(t))
# Add codes 1,3, since complements of the homogeneous runs occur.
T,P0=build(M)
P={}
for d in "0123":
 powers=[np.eye(36,dtype=np.int64)]
 for n in range(1,30):powers.append(powers[-1]@T[int(d)]%M)
 assert np.all(powers[5]==powers[29])
 P[d]=(powers,5,24)
Z=set()
def toggle(t):
 if isinstance(t,tuple): return (str(int(t[0])^1),t[1])
 return str(int(t)^1)
def cut(l,r):
 w=r+tuple(map(comp,l))
 for u in (w,tuple(map(comp,w))):
  Z.add(int(raw(T,P,u,M,1)))
  Z.add(int(raw(T,P,tuple(map(toggle,u)),M,0)))
for i in range(7):cut(F[:i],F[i:])
for i,d in ((2,'0'),(3,'2')):
 for a,b in product(range(1,29),repeat=2):
  if max(a,b)<5 or (a+b)%24:continue
  cut(F[:i]+((d,a),),((d,b),)+F[i+1:])
print(len(Z),sorted(Z))
assert 3 not in Z and len(Z)==49
