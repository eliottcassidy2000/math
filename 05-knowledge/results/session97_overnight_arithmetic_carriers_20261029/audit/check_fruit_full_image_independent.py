"""Independent finite arithmetic supporting full-image conditional theorem."""
from fractions import Fraction as F
import json

def need(b,m):
    if not b:raise ArithmeticError(m)
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;z,w=Q
    if x==z and y==-w:return None
    m=(3*x*x+842*x+43904)/(2*y) if P==Q else (w-y)/(z-x)
    xx=m*m-421-x-z
    return xx,m*(x-xx)-y

def times(n,P):
    if n<0:return times(-n,(P[0],-P[1]))
    out=None
    for bit in bin(n)[2:]:
        out=add(out,out)
        if bit=='1':out=add(out,P)
    return out

def reduce3(P):
    if P is None or P[0].denominator%3==0:return None
    return tuple(v.numerator*pow(v.denominator,-1,3)%3 for v in P)

def main():
    discr=16*43904**2*(421**2-4*43904);counts=[]
    for p in [3,11,17]:
        need(discr%p!=0,'good reduction')
        counts.append(1+sum((y*y-x**3-421*x*x-43904*x)%p==0 for x in range(p) for y in range(p)))
    need(counts==[6,12,18],'torsion reduction counts')
    need(all((x**3+300*x*x+9408*x-2458624)%19!=0 for x in range(19)),'rootless cubic')
    z=-F(784,3);need(z**3+421*z*z+43904*z==-F(15366400,27),'real exclusion')
    coeff=[0]*5
    for i,a in enumerate([784,3]):
        for j,b in enumerate([-2458624,9408,300,1]):coeff[i+j]+=a*b
    need(coeff==[-43904**2,0,6*43904,4*421,3],'division polynomial')
    H=(F(-196),F(196));U=(F(0),F(0))
    orbit=[reduce3(times(n,H)) for n in range(6)]
    need(len(set(orbit))==6 and reduce3(times(6,H)) is None,'H mod3 generator')
    need(orbit[3]==(0,0),'3H mod3=U')
    tested=0
    for n in range(-40,41):
        for e in [0,1]:
            Q=add(times(n,H),U if e else None)
            small=reduce3(Q)
            need((small in [None,(0,0)])==(n%3==0),'global-to-reduction criterion on named bank')
            if Q is not None:
                q=-Q[0]/196-1
                guard=q.denominator%3==0 or (q.numerator+q.denominator)%3==0
                need(guard==(n%3==0),'cheap q-coordinate guard')
            tested+=1
    print(json.dumps({'good_prime_counts':dict(zip([3,11,17],counts)),
        'mod3_H_orbit':orbit,'checked_named_targets':tested,
        'global_input':'CITED E(Q)=ZG+C6, not independently saturated'},sort_keys=True))
if __name__=='__main__':main()
