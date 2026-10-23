from itertools import permutations,product
from collections import Counter
from fractions import Fraction
from math import comb
import json
G=list(permutations(range(3)));one=(0,1,2)
def mul(a,b):return tuple(a[b[i]] for i in range(3))
def inv(a):return tuple(a.index(i) for i in range(3))
def parity(a):return sum(a[i]>a[j] for i in range(3) for j in range(i+1,3))%2
def char(a):return (1,(-1)**parity(a),sum(a[i]==i for i in range(3))-1)
dims=(1,1,2);checks=0
def eq(a,b):
 global checks
 checks+=1
 if a!=b:raise RuntimeError((a,b))
# Full character convolution table and all group arguments, independently direct.
for r,s,g in product(range(3),range(3),G):
 lhs=sum(Fraction(char(x)[r]*char(mul(inv(x),g))[s],6) for x in G)
 eq(lhs,Fraction(char(g)[r],dims[r]) if r==s else 0)
hist=Counter()
for a,b,c in product(G,repeat=3):
 faces=(a,b,c,mul(mul(a,b),inv(c)))
 hist[sum(x==one for x in faces)]+=1
poly=[Fraction(comb(4,k)*(5**(4-k)+5*(-1)**(4-k)),6) for k in range(5)]
eq([hist[k] for k in range(5)],poly)
# arbitrary real central inversion symmetric weights, including signed hostile values
weighted=[]
for weights in [(2,3,5),(0,1,-2),(7,-3,4)]:
 def w(g):return weights[0 if g==one else 1 if parity(g) else 2]
 direct=sum(Fraction(w(a)*w(b)*w(c)*w(mul(mul(a,b),inv(c))),216) for a,b,c in product(G,repeat=3))
 coeff=[sum(Fraction(w(g)*char(g)[r],6*dims[r]) for g in G) for r in range(3)]
 spectral=sum(dims[r]**2*coeff[r]**4 for r in range(3))
 eq(direct,spectral);weighted.append({'class_weights':weights,'direct_and_spectral':str(direct),'a_coefficients':list(map(str,coeff))})
print(json.dumps({'status':'PASS','checks':checks,'S3_identity_face_histogram':dict(sorted(hist.items())),'rooted_configurations':216,'arbitrary_central_weight_controls':weighted},indent=2))
