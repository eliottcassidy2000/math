from fractions import Fraction as Q
from math import isqrt
import sys
from pathlib import Path
# Verify through independent collection of the Farey read.
from target_decode import F

def setup(m):
    q=1<<m
    h0=q+1; p0=4*q+2*(3**m-1); L=2*q-1
    return 5*h0,256*p0,147*p0-10*h0,142*h0+256*L,q,h0,p0,L

def check(m,G,delta):
    h,P,A,B,q,h0,p0,L=setup(m)
    r=G//2+1;s=G-1-r
    rev=lambda n:int(format(n,'0'+str(G.bit_length()-1)+'b')[::-1],2)
    C=lambda n:G-1-rev(n)
    bias=20*P*G*G-4*P*G*(s-r)+2*h-4*h*(C(r)-C(s))-4*G*(A-B)
    # y=2k+1: nearest odd, with integer-only round.
    num=delta+(64*P*G*G-h)+bias
    y=2*(num//(6*h))+1
    x=y+1
    u,w=G*(16*x+5)+r,G*(16*y+4)+s
    def H(u,w):
        lk=p0*u+L
        d=2*q-2+h0*(F(u)-1)+(p0-2*h0)*u+q
        return 91*h0+5*(d+h0*(F(w)-1))+142*(lk+h0*w)+256*lk*w
    d=H(u+4*G,w-4*G)-H(u,w)
    if d!=3*h*y-(64*P*G*G-h)-bias or abs(d-delta)>3*h:
        raise ArithmeticError('rounding or physical count failed')
    exact=Q(H(u,w),G**6)
    lead=Q(448**2*P**3,h*h)
    if abs(exact-lead)*G>100*lead:raise ArithmeticError('asymptotic scale failed')

for m in range(2,60,2):
    for t in (8,17,32,70):
        for d in (-2000,0,4,40000):
            check(m,1<<t,d)
phases=[]
for m in range(2,202,2):
    h,P,*_=setup(m)
    c=Q(448**2*P**3,h*h)
    while c>=64:c/=64
    phases.append(c)
phases.sort()
maxratio=max(y/x for x,y in zip(phases,phases[1:]+[64*phases[0]]))
if maxratio>Q(11,10):raise ArithmeticError('100-phase net too wide')
print('464 exact strong-flank scale controls; 100 phases give adjacent count ratios below 1.10.')
