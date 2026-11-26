
from fractions import Fraction
from math import ceil

D=[b for b in range(1,57) if 56%b==0 and (b-1)%5==0]
assert D==[1,56]
rows=[]
for a in range(4,54):
    u=max(a*a+1,2*a+1,2*a*(57-a)-(a-1)-a*(a-1)//2-(a*a-4*a))
    assert u<=3024
    old=ceil(Fraction(a*u*(6695-u),130*(53+a)))
    new=ceil(Fraction(a*u*(103*3024-46*u),6048*(53+a)))
    rows.append((old,new,a,u))
assert min(rows)==(1407,1408,4,415)
assert min(x[1] for x in rows)==1408
print("First rows (old bound, exterior bound, min loss degree, support):",rows[:5])
print("Compatible holonomy block counts:",D)
S=[((57-m)*(m+1)+1)//2 for m in range(3,28)]
assert S.count(108)==1 and S[0]==108
assert 2*(53**2+1)>54*56 and 108*56==6048
