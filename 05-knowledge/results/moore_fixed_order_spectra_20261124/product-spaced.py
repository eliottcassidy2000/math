from math import isqrt, prod
charts = [(29,640,191,158,111),(651,17616,6768,2691,2313)]
def q(n):
    w=0
    while n:
        if n&1:w+=(n//2+1)//2
        n//=2
    return w

def T(h,c,X,Y):
    return 1 << (2*c*Y+h*(X+Y)).bit_length()

for h,a,b,c,k in charts:
    for m in range(3,8):
        Y=1<<m
        t=T(h,c,Y,Y)
        shift=2*h*Y
        products=set()
        for x in range(Y,2*Y):
            v=t*x
            assert q(v)==q(x)
            for u in range(Y,2*Y):
                g=a*v*u+b*v+c*u-h*(q(v)+q(u))
                p=x*(a*u+b)
                assert (g+shift)//t==p
                products.add(p)
        if m==7:
            print('balanced box',h,Y,'distinct product quotients',len(products))

# A primitive-root/order certificate for the only large prime used.
p=61
MP=(1<<p)-1
fac={2:1,3:2,5:2,7:1,11:1,13:1,31:1,41:1,61:1,151:1,331:1,1321:1}
assert prod(r**e for r,e in fac.items())==MP-1
for r in fac:assert all(r%i for i in range(2,isqrt(r)+1))
assert pow(37,MP-1,MP)==1 and all(pow(37,(MP-1)//r,MP)!=1 for r in fac)
E=(1<<(p-1))*MP
for h,a,B,C,k in charts:
    for b,c in [(B,C),(C,B)]:
        for delta in ([-3,-1,1] if h==29 else [-1]):
            G=(E+delta-k)//2
            tested=0
            for i in range(p):
                X=1<<i
                for j in range(p):
                    Y=1<<j
                    shift=h*(X+Y)
                    t=T(h,c,X,Y)
                    if not t*X*(a*Y+b)<=G+shift<t*((2*X-1)*(a*(2*Y-1)+b)+1):continue
                    tested+=1
                    assert t<=1<<(p-2)
                    assert (G+shift)//t==(1<<(p-2))*MP//t
                    ymax=a*(2*Y-1)+b
                    assert 2*X-1<MP and (ymax<MP or (a%2==b%2==0 and ymax<2*MP))
            print('E61 necessary product boxes',h,c,delta,tested)
