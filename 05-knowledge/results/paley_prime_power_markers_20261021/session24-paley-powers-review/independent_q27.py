"""Independent q=27 arithmetic and all-three-subset stabilizer census."""
import itertools,json
from collections import Counter
# F3[z]/(z^3+2z^2+1), separate arithmetic without importing the producer.
p=3;q=27
D=[(n%3,(n//3)%3,n//9) for n in range(q)]
def enc(v): return sum((a%3)*3**i for i,a in enumerate(v))
def add(x,y):return enc([a+b for a,b in zip(D[x],D[y])])
def mul(x,y):
 c=[0]*5
 for i,a in enumerate(D[x]):
  for j,b in enumerate(D[y]): c[i+j]+=a*b
 for i in (4,3):
  c[i-3]-=c[i];c[i-1]-=2*c[i]
 return enc(c[:3])
def power(x,k):
 out=1
 for _ in range(k):out=mul(out,x)
 return out
perms=[];labels=[]
for j in range(3):
 for a in range(1,27):
  for b in range(27):
   perm=tuple(add(mul(a,power(x,3**j)),b) for x in range(27))
   if len(set(perm))!=27:raise RuntimeError('not permutation')
   perms.append(perm);labels.append((a,b,j))
if len(set(perms))!=2106:raise RuntimeError('duplicate group map')
counts=Counter()
for s in itertools.combinations(range(q),3):
 ss=set(s)
 stabilizers=[k for k,f in enumerate(perms) if {f[x] for x in s}==ss]
 counts[len(stabilizers)]+=1
s=(0,1,5)
hostile=[(labels[k],tuple(f[x] for x in s)) for k,f in enumerate(perms) if {f[x] for x in s}==set(s)]
if counts!={1:2106,3:702,18:117}:raise RuntimeError(str(counts))
if sum(label[2]==0 for label,im in hostile)!=1 or len(hostile)!=3:raise RuntimeError('hostile')
print(json.dumps({'modulus_low_to_high':[1,0,2,1],'q':27,'all_three_subset_stabilizer_histogram':dict(sorted(counts.items())),'free_signed_triple_orbits':counts[1]//2106,'affine_regular_semilinear_hostile':{'S':s,'t_polynomial':'2+z','stabilizers_and_images':hostile}},indent=2,sort_keys=True))
