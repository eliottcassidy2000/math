from itertools import product
from fractions import Fraction

def v2(n):
    n=abs(n)
    if not n: return None
    return (n & -n).bit_length()-1

def coeff(word,p,sign):
    A,B,C=1,1,0
    for a in word:
        A,C,B=p*A,p*C+sign*B,B*2**a
    return A,B,C,A-B

def step(n,p,sign):
    z=p*n+sign
    a=v2(z)
    if a is None: raise ValueError('zero')
    return z//2**a,a

def actual_reps(n,w,p,sign,cap=100):
    t=0
    x=n
    while t<cap:
        for a in w:
            y,k=step(x,p,sign)
            if k!=a: return t
            x=y
        t+=1
    return None

count=0
zero=0
switch=0
for p in (3,5):
 for sign in (1,-1):
  words=[w for r in range(1,5) for w in product(range(1,5),repeat=r)]
  bank=[coeff(w,p,sign) for w in words]
  for w,(A,B,C,D) in zip(words,bank):
   S=sum(w)
   for n in range(1,2002,2):
    e=D*n+C
    predicted=None if e==0 else (v2(e)-1)//S
    actual=actual_reps(n,w,p,sign)
    if actual!=predicted: raise AssertionError((p,sign,w,n,e,predicted,actual))
    count+=1
    if e==0: zero+=1
    if predicted and predicted>0:
     y=(A*n+C)//B
     # deterministic second word selected independently by source index
     av,bv,cv,dv=bank[(n+len(w))%len(bank)]
     delta=D*cv-dv*C
     lhs=D*B*(dv*y+cv)
     rhs=A*dv*e+B*delta
     if lhs!=rhs: raise AssertionError('switch')
     switch+=1
print('repetition_checks',count,'zero_defect_checks',zero,'switch_identity_checks',switch)
for H in range(5,126,6):
 n=(2**(H+3)-13)//9
 if 9*n!=2**(H+3)-13: raise AssertionError('integrality')
 x,k=step(n,3,1); y,l=step(x,3,1)
 if (k,l)!=(1,2) or y!=2**H-1: raise AssertionError('path')
 if v2(n+5)!=5 or v2(y+1)!=H: raise AssertionError('budget')
 if actual_reps(n,(1,2),3,1)!=1: raise AssertionError('old repeats')
 if actual_reps(y,(1,),3,1,cap=H+2)!=H-1: raise AssertionError('new repeats')
print('reset_family_checks',len(range(5,126,6)))
for p,sign,w,n in [(3,1,(1,),-1),(3,-1,(1,2),5),(5,1,(1,1,5),13),(3,1,(2,),1)]:
 A,B,C,D=coeff(w,p,sign)
 print('hostile',p,sign,w,n,'D',D,'C',C,'E',D*n+C)
 if D*n+C: raise AssertionError('control cycle')
for H in (5,11,17):
 n=(2**(H+3)-13)//9
 x,a=step(n,3,1); y,b=step(x,3,1)
 print('reset_example',H,(n,x,y),'old_budget',1,'new_budget',H-1,'hidden_center_numerator',9*n+13)
