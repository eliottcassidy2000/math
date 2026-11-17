from fractions import Fraction as F
pairs = [(F(1,2), F(504,1000)), (F(1,12),F(86,1000)),
         (F(2,3),F(670,1000)),(F(4,7),F(575,1000)),
         (F(8,17),F(475,1000)),(F(16,51),F(318,1000))]
if not all(a**99 < b**100 for a,b in pairs):
    raise RuntimeError("invalid rational upper bound")
x,b,a0,a1,a2,a3 = (q for p,q in pairs)
l,l2,l3=(x**i/(1-x**i) for i in (1,2,3))
two=(l*l+l2)/2
three=(l**3+3*l*l2+2*l3)/6
t=F(4,5)
root=a0*three+a1*two*t+a2*l*(t*t+b*t)/2+a3*(t**3+3*b*t*t+2*b*b*t)/6
spine=a1*two+a2*l*t+a3*(t*t+b*t)/2
if not (root<t and spine<1):
    raise RuntimeError((root,spine))
print("exact root < 4/5 and distinguished branch < 1")
