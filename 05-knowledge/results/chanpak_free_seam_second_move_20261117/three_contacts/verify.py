from random import Random

def F(n):
    u,z=0,1
    for e in map(int,bin(n)[2:]):
        if not e:z+=u+(u%2)
        u=2*u+e
    return z

def H(u,w):
    return 81+142*(u+w)+256*u*w+5*(F(u)+F(w))+2*(w%2)

def rc(r,g):
    return g-1-int(format(r,'0'+str(g.bit_length()-1)+'b')[::-1],2)

def check(g,r,d):
    s=g-1-r
    A=16384*g*g-5
    B=5120*g*g-1024*g*(s-r)+10-20*(rc(r,g)-rc(s,g))
    z=(-(d+B)*pow(A,-1,15))%15 or 15
    y=(d+A*z+B)//15
    x=y+z
    u,w=g*(16*x+5)+r,g*(16*y+4)+s
    v,t=u+4*g,w-4*g
    if (x%2 or min(x,y)<=0 or (u^v).bit_count()+(w^t).bit_count()!=3
        or (u+w+1)%g or v+t!=u+w or H(v,t)-H(u,w)!=d
        or H(w,u)-H(u,w)!=2 or H(t,v)-H(v,t)!=2):
        raise ArithmeticError((g,r,d,x,y))
    return H(u,w),H(v,t)

if check(4,1,4)!=(22047013898432061,22047013898432065):raise ArithmeticError('fixture')
R=Random(16384)
for gg in range(2,65):
    g=1<<gg
    for _ in range(40):
        r=1+2*R.randrange(g//2)
        d=2*R.randrange(-g*g//2,g*g//2+1)
        check(g,r,d)
print('Exact example and 2520 complementary-suffix controls pass.')
