"""Exact arithmetic certificates; no assertions disabled by Python -O."""
from fractions import Fraction as F
from math import gcd, isqrt
import json

def check(condition, label):
    if not condition:
        raise RuntimeError(label)

def sf(n):
    out, p = 1, 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e ^= 1
        if e: out *= p
        p += 1
    return out*n

def roots_at(N,d):
    D=N*d*d
    out=[]
    for u in range(1,isqrt(D)+1):
        if D%u: continue
        v=D//u
        if (u+v)%2: continue
        p,q=(v-u)//2,(u+v)//2
        r=isqrt(q*q+D)
        if p>0 and r*r==q*q+D:
            out.append([p,q,r])
    return out

s16={x*x%16 for x in range(16)}
steps16=sorted({h for a in s16 for h in range(16)
                if (a+h)%16 in s16 and (a+2*h)%16 in s16})
s3={x*x%3 for x in range(3)}
steps3=sorted({h for a in s3 for h in range(3)
               if (a+h)%3 in s3 and (a+2*h)%3 in s3})
check(steps16==[0,8], 'mod16 classification')
check(steps3==[0], 'mod3 classification')
check(roots_at(5,12)==[[31,41,49]], 'least denominator uniqueness')
check(all(not roots_at(5,d) for d in range(1,12)), 'independent small denominator search')
p,q,r,d=31,41,49,12
check(q*q-p*p==r*r-q*q==5*d*d,'root gaps')
a,b,c=F(r-p,d),F(r+p,d),F(2*q,d)
check(a*a+b*b==c*c and a*b/2==5,'right triangle')
x=5*(a+c)/b
y=10*x/b
m=(3*x*x-25)/(2*y)
x2=m*m-2*x
y2=m*(x-x2)-y
check((x,y)==(F(25,4),F(75,8)),'triangle point')
check((x2,-y2)==(F(q*q,d*d),F(p*q*r,d**3)),'signed doubling')
X=F(q*q,d*d); N=5
check(X>2*N, 'direct AP positivity')
check(not ((X*X+N*N)**2>8*N*X*(X*X-N*N)), 'triangle chart positivity failure')
labels=((X-2*N)/2,X/2,(X+2*N)/2)
check(tuple(576*z for z in labels)==(482,3362,6242), 'direct chart scaled labels')
translate_x=[x,-F(N*N,x),N*(x+N)/(x-N),-N*(x-N)/(x+N)]
check(translate_x==[F(25,4),F(-4),F(45),F(-5,9)], 'two-torsion translates')
def rational_square(z):
    return z>=0 and isqrt(z.numerator)**2==z.numerator and isqrt(z.denominator)**2==z.denominator
check(all(not all(rational_square(z+h) for h in (-N,0,N)) for z in translate_x), 'Q not a rational quadruple')
spine=[(3,4,5),(5,12,13),(7,24,25),(9,40,41)]
classes=[sf(A*B//2) for A,B,C in spine]
check(classes==[6,30,21,5],'Berggren classes')
check([C*C>4*A*B for A,B,C in spine]==[False,False,False,True],'positivity exit')
# Independent parameter census: tests the signed-doubling identity on all
# positive ordered triangles with smaller leg designated a, m <= 40.
count=0
for M in range(2,41):
    for n in range(1,M):
        if gcd(M,n)!=1 or (M-n)%2==0: continue
        A,B=sorted((M*M-n*n,2*M*n)); C=M*M+n*n
        N=F(A*B,2)
        X=N*F(A+C,B); Y=2*N*X/B
        slope=(3*X*X-N*N)/(2*Y)
        X2=slope*slope-2*X
        Y2=slope*(X-X2)-Y
        check(X2==F(C*C,4), 'general doubled x')
        check(Y2==F((A*A-B*B)*C,8), 'general doubled y')
        count+=1
print(json.dumps({'residue_steps_mod16':steps16,'residue_steps_mod3':steps3,
 'gap5_least_denominator':12,'unique_positive_roots_at_denominator12':roots_at(5,12),
 'triangle':list(map(str,(a,b,c))),'P':list(map(str,(x,y))),
 'Q_minus_2P':list(map(str,(x2,-y2))), 'Berggren_spine_area_classes':classes,
 'Q_in_2E_not_4E':True,'distinct_positivity_chambers_checked':True,
 'primitive_triangle_identity_checks':count},indent=2))
