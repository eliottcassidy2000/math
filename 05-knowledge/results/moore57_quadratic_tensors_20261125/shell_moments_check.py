from fractions import Fraction as Q
from math import ceil,floor

allmom=[2700,0,2700*47,0,2700*47*93]
F=lambda j:47**j+1189*7**j+980*(-8)**j
m=[x-F(j) for j,x in enumerate(allmom)]
assert m==[530,-530,3710,-9890,53150]
def certificate(u,v,q):
    # (u+v*x) q(x)^2 has leading coefficient v.
    low=sum(q[i]*q[j]*(u*m[i+j]+(v*m[i+j+1] if i+j<4 else 0))
            for i in range(3) for j in range(3))
    return Q(-low,v)
upper=certificate(6,-1,[-Q(215472,59837),Q(3079,1129),1])
lower=certificate(7,1,[-Q(268,53),1,1])
assert upper==-Q(6202311130,59837) and lower==-Q(14556970,53)
room=(ceil((F(5)+lower)/10),floor((F(5)+upper)/10))
assert room==(21694124,21711223)
print("Conditional pentagon interval on residual shell:",room)
