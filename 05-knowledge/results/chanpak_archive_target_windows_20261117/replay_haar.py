from slot_profiles import *
D=(20,121,114,245,32,1,49,119)
a,b=profile(D,77),profile(D,146)
if (a,b)!=((51,1512,705,534,15,22,-163),(51,1512,705,534,19,28,-155)):raise ValueError("tiles")
sentinel=profile((2,0),1),profile((0,),0)
for A,B,expected in ((a,b,-2064993),(b,a,-2064991)):
 rows=(1,0,1,0)
 for p in (sentinel[0],A,B,sentinel[1]):rows=append(rows,p)
 if rows!=(2601,11348406,expected,0):raise ValueError(rows)
 # Strength is tested on the entire literal exterior, not on its scalar profile.
 t=ordinal((2,0),ordinal(D,ordinal(D,(0,))))
 Q=1|((77 if A==a else 146)<<2)|((146 if A==a else 77)<<10)
 if not strong(external(t,Q,1)):raise ValueError("strong")
 print(rows)
print("pure mixed determinant",a[4]*b[5]-b[4]*a[5])
from math import isqrt
p=11348406*4-2064991
if p!=43328633 or any(p%d==0 for d in range(2,isqrt(p)+1)):
    raise ValueError("prime read")
print("prime",p)
