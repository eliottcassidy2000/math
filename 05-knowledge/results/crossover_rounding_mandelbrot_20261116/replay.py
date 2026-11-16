from fractions import Fraction
from itertools import combinations, permutations
from math import isqrt, factorial

def F(z):
    out=1
    while z:
        a=z//2
        if not z%2:
            out+=a+(a%2)
        z=a
    return out

def run(z,l):
    """Upper analytic envelope from a reverse address and l carry bits per run."""
    i=0
    W=Fraction()
    D=Fraction(1)
    R=0
    while i<z.bit_length()-1:
        if z>>i&1:
            i+=1
            continue
        a=i
        while not z>>i&1:
            i+=1
        # run [a,i); alpha=(z mod 2^a)/2^a
        w=1-Fraction(1,2**(i-a))
        W+=Fraction(w,2**a)
        alpha_low=Fraction(((z%(2**a))*(2**l))//(2**a),2**l)
        D+=1-alpha_low*w
        R+=1
    return z*W+D, R

def local_residue(z,m):
    mod=2**m
    y=1
    for i in range(z.bit_length()):
        if not z>>i&1:
            u=z>>(i+1)
            y+=(u&(mod-1))+(u&1)
    return y%mod

def recover(z,A,g,J,l,m):
    M=2**m
    fv,R=run(z,l)
    if not g*R <= 2**(m+l):
        raise ValueError("precision budget too small")
    U=A*(z+1)+g*fv+J
    residue=(A*(z+1)+g*local_residue(z,m)+J)%M
    return residue+M*((U-residue)//M)

def components(n):
    arr=[0]*(n+1)
    for j in range(1,n+1):
        arr[j]=2**(j-1)-sum(arr[d] for d in range(1,j) if j%d==0)
    return arr[-1]

def prime(p):
    return p>1 and all(p%d for d in range(2,isqrt(p)+1))

def one_apex_h(A):
    return 1+sum(2**(t-1) for t in A)-sum(2**max(t-s-1,1) for s in A for t in A if s<t)

def class_counts(n):
    edges=list(combinations(range(n),2))
    total=[0,0]
    for p in permutations(range(n)):
        seen=set()
        orbit_lengths=[]
        xor_parities=[]
        for edge in edges:
            if edge in seen:
                continue
            t=edge
            k=x=0
            while t not in seen:
                seen.add(t)
                u,v=p[t[0]],p[t[1]]
                t=(min(u,v),max(u,v))
                x^=int(u>v)
                k+=1
            xor_parities.append(x)
            orbit_lengths.append(k)
        mass=2**len(orbit_lengths)
        for sign in range(2):
            if all(x==(sign*k)%2 for x,k in zip(xor_parities,orbit_lengths)):
                total[sign]+=mass
    return tuple(t//factorial(n) for t in total)

def test(condition,*args):
    if not condition:
        raise ArithmeticError(args)

if __name__=='__main__':
    for z in range(512):
        for m in range(1,7):
            for l in range(5):
                v,R=run(z,l)
                test(F(z)<=v and (v==F(z) if R==0 else v-F(z)<Fraction(R,2**l)),z,l,m)
                test(local_residue(z,m)==F(z)%2**m,z,l,m)
                if 11*R<=2**(l+m):
                    test(recover(z,52,11,0,l,m)==52*(z+1)+11*F(z),z,l,m)
    for P,A,g,J,z in [(30011,18,3,-4,1561),(30013,12,5,0,2363),
                      (1000000007,18,3,-4,47827936),(1000000009,32,5,-4,30268243)]:
        _,R=run(z,0)
        l=0
        m=(g*R-1).bit_length() if R else 1
        test(recover(z,A,g,J,l,m)==P,P)
        l=m-1
        U=A*(z+1)+g*run(z,l)[0]+J
        test(recover(z,A,g,J,l,1)==P,P)
        print(P,'runs',R,'carry bits',l,'upper overshoot',U-P)
    for r in [n for n in range(3,60) if prime(n)]:
        k=r-1
        a=k//2;b=k+1-a
        source=3 if r==3 else one_apex_h((a,b,k+1))
        test(components(r)==source,r)
        test(components(2*r)==(2**r+1)*source,r)
    for p in [2,3,5,7]:
        for q in [11,13,17]:
            x=2**q
            y=sum(2**j*x**(p-1-j) for j in range(p))-1
            test(components(p*q)==components(q)*y,p,q)
    table=[class_counts(n) for n in range(1,8)]
    test(table==[(1,1),(1,1),(2,2),(4,2),(12,8),(56,12),(456,88)],table)
    print("Converse census",table)
    print("Mandelbrot 5, 6, 7, 10, 14", [components(n) for n in (5,6,7,10,14)])
    print("PASS")
