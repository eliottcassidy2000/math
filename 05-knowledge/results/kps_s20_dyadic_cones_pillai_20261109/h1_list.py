# One-step finiteness (THM-5109 (E)): a sheet-compatible dyadic transitive cone step from value h needs
# h = 2^t - 3^a with 1 <= h <= 2^a (because a >= j and 2^j = B >= h).  For t >= 28 Ellison gives
# 2^t - 3^a > 2^t e^(-t/10) > 2^(0.8557 t) > 2^(0.631 t) > 2^a, so only t <= 27 can occur.  We list t <= 60 anyway.
import math
sol=[]
for t in range(1,61):
    a=0
    while 3**a < 2**t:
        h=2**t-3**a
        if 1<=h<=2**a: sol.append((t,a,h))
        a+=1
print('solutions (t,a,h) of 1 <= 2^t-3^a <= 2^a, t<=60:',sol)
# Ellison sanity: for 28<=t<=60 check 2^t-3^a > 2^t*exp(-t/10) for all a with 3^a<2^t
bad=[(t,a) for t in range(28,61) for a in range(0,int(t*math.log(2)/math.log(3))+1) if 3**a<2**t and not (2**t-3**a > 2**t*math.exp(-t/10))]
print('Ellison bound violations for 28<=t<=60:',bad)
