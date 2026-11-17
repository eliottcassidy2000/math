"""Physically compile arbitrary (p,q) residues, for moduli coprime to 6."""
from math import gcd

A=((0,1),)
B=((1,1),)
def inv(w):return tuple((s,-d) for s,d in reversed(w))
def power(w,k):return (w if k>=0 else inv(w))*abs(k)
W=inv(A)*2+B
def conj(j):return power(A,-j)+W+power(A,j)
def comm(j):return conj(j)+W+inv(W+conj(j))
V=(comm(1),comm(2),sum((power(conj(j),e) for j,e in enumerate((-6,19,-16,4))),()))

def step(v,m,cmd):
 f,g,j,d,x,y,par,f3,x3=v
 i2=pow(2,-1,m); i4=i2*i2
 if cmd==0:
    z=(2*f,3*g*i2,j*i2,2*d,x,y)
    side=(1-par,2*f3%3,x3)
 else:
    z=(4*f+2,(9*g+3)*i4,(j+1)*i4,4*d,x-f,y-g-j)
    side=(par,(f3+2)%3,(x3-f3)%3)
 return tuple(a%m for a in z)+side

def run(v,m,w,output=None):
 for s,e in w:
    if e<0:
        u=step(v,m,s);period=1
        while u!=v:u=step(u,m,s);period+=1
        k=(period-1)*(-e)
    else:k=e
    for _ in range(k):v=step(v,m,s)
    if output is not None:output.append(('0' if s==0 else '10')*k)
 return v

def direct(word,m):
 # Separate unnormalized word recurrence for literal h and P2.
 n=0;D=1;f=0;G=J=0;h=1;P=2;last=0
 for b in map(int,word):
    if not b:h=(h+f+last)%m;P=(2*P+2*G+2*J)%m
    else:P=2*P%m
    f=(2*f+b)%m;G=(3*G+b*D)%m;J=(J+b*D)%m
    D=2*D%m;last=b;n+=1
 return h,P

def pair(v,m):
 f,g,j,d,x,y=v[:6]
 return (x+f)%m,d*(y+2*g-2*j)%m

def solve(v,w,t,m):
 if m==1:return (0,0)
 a,b=v;c,d=w;e,f=t
 I=pow((a*d-b*c)%m,-1,m)
 return I*(d*e-c*f)%m,I*(-b*e+a*f)%m

def initial(word,m):
 n=len(word)
 d=pow(2,n,m)
 inv=pow(d,-1,m)
 f=int(word,2)
 g=sum(int(b)*2**i*3**(n-i-1) for i,b in enumerate(word))*inv % m
 j=sum(int(b)*2**i for i,b in enumerate(word))*inv % m
 h,p=direct(word,m)
 h3,_=direct(word,3)
 return (f%m,g,j,d,(h-f)%m,(p*inv-2*g+2*j)%m,n%2,f%3,(h3-f)%3)

def compile_profile(m,h,P,p_mod3):
 if m<2 or gcd(m,6)!=1:raise ValueError('modulus must be >1 and coprime to 6')
 if p_mod3 not in (1,2):raise ValueError(p_mod3)
 base='10100' if p_mod3==1 else '100'
 start=initial(base,m); phases=tuple(power(w,3) for w in V)
 changes=[]
 for word in phases:
    a=run(start,m,word)
    if a[:4]+a[6:]!=start[:4]+start[6:]:raise ValueError('register shift')
    changes.append(tuple((x-y)%m for x,y in zip(a[4:6],start[4:6])))
 if changes[:2]!=[(-9*pow(4,-1,m)%m,33*pow(9,-1,m)%m),(-45*pow(8,-1,m)%m,255*pow(27,-1,m)%m)]:raise ValueError(changes)
 if m%5==0 and (changes[2][1]+2*changes[2][0])%5!=1:raise ValueError(changes)
 H0,P0=pair(start,m)
 target=((h-H0)%m,(P-P0)*pow(start[3],-1,m)%m)
 m1=1;m0=m
 while m0%5==0:m0//=5;m1*=5
 a,b=solve(*changes[:2],target,m0)
 c,d=solve(changes[0],changes[2],target,m1)
 z0=(a,b,0);z1=(c,0,d)
 powers=[(x+m0*((y-x)*pow(m0,-1,m1)%m1))%m if m1>1 else x for x,y in zip(z0,z1)]
 words=[base];state=start
 for op,n in zip(phases,powers):state=run(state,m,power(op,n),words)
 word=''.join(words)
 if not word.startswith('1') or not word.endswith('0') or not len(word)%2:raise ValueError('not strong even order')
 if direct(word,m)!=(h%m,P%m) or pair(state,m)!=(h%m,P%m):raise ValueError((m,h,P,state))
 if direct(word,3)!=(p_mod3,p_mod3):raise ValueError('bad mod3 phase')
 return word,powers

if __name__=='__main__':
 for m,p,q,c in ((5,4,1,1),(7,0,0,2),(11,2,5,2),(25,9,21,1),(35,14,17,2)):
    h,P=(71*p-44*q)%m,(21*p-13*q)%m
    word,coeff=compile_profile(m,h,P,c)
    h1,P1=direct(word,m)
    if ((44*P1-13*h1)%m,(71*P1-21*h1)%m)!=(p%m,q%m):raise ValueError('wrong pair')
    print(m,p%m,q%m,c,h1,P1,'word length',len(word),'powers',coeff)
