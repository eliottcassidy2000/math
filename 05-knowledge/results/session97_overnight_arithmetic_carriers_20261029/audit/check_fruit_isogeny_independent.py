"""Exact independently implemented binary multiples and local guard audit."""
from fractions import Fraction as R
import json

def need(b,m):
    if not b:raise ArithmeticError(m)

def plus(P,Q,A,B):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u:
        if y+v==0:return None
        slope=(3*x*x+2*A*x+B)/(y+v)
    else:slope=(y-v)/(x-u)
    X=slope*slope-A-x-u
    return X,-y+slope*(x-X)

def times(n,P,A,B):
    out=None
    for bit in bin(n)[2:]:
        out=plus(out,out,A,B)
        if bit=='1':out=plus(out,P,A,B)
    return out

def phi(P):
    x,y=P
    return x*(x-56)**2/(x-4)**2,(x-56)*(x*x+44*x+224)*y/(x-4)**3

def main():
    G=(R(-4),R(28));H=(R(-196),R(196))
    need(times(3,H,421,43904)==phi(G)==(R(-225),R(210)),'3H identity')
    g2=times(2,G,109,224)
    need(g2[0]==R(676,49),'nonintegral doubled point')
    # Integral short model: X=9x+3A, Y=27y, a=81B-27A²,
    # b=54A³-243AB. Nonintegral doubled X rules out torsion by Nagell-Lutz.
    X,Y=9*g2[0]+327,27*g2[1]
    aa=81*224-27*109**2;bb=54*109**3-243*109*224
    need(Y*Y==X**3+aa*X+bb,'short model transport')
    need(X==R(22107,49) and X.denominator>1,'Nagell-Lutz obstruction')
    for n in range(1,11):
        need(phi(times(n,G,109,224))==times(3*n,H,421,43904),'homomorphism test')
    Z,Y=times(35,H,421,43904)
    need(Y*Y==Z**3+421*Z*Z+43904*Z,'35H curve')
    q=-Z/196-1;W=-Y/196
    need(q>R(1,6) and 1+29*q-196*q*q>0,'strict positivity chamber')
    need(W*W==(q+1)*(1+29*q-196*q*q),'square discriminant')
    need(q.numerator%2==0 and q.denominator%2==1,'mod2 coefficient reduction')
    need(all((t**3+t*t+1)%2 for t in [0,1]),'local root obstruction')
    print(json.dumps({'3H':[-225,210],'doubleG_x':str(g2[0]),
        'short_model_doubleG_X':str(X),'short_model_coefficients':[aa,bb],
        'positive_hostile_index':35,'hostile_q_parities':[q.numerator%2,q.denominator%2],
        'checked_map_multiples':10},sort_keys=True))
if __name__=='__main__':main()
