from random import Random

B,Q,T,S=2855,5760,4034,3914
A=1<<36
U,V=46767159475,48234217403
def need(c):
    if not c: raise AssertionError("unit transfer")
def p(n):
    ans=0
    while n:
        n,b=divmod(n,2)
        if not b:ans+=n+(n&1)
    return ans
def rev(n,d):
    u=0
    for _ in range(d):
        n,b=divmod(n,2)
        u=2*u+b
    return u
def F(x,y):
    return B+Q*x*y+T*x+S*y+45*(p(x)+p(y))
def coeff(a,b,d,t,s,delta):
    l=1<<d
    need(l//2<=a<l and l//2<=b<l)
    need(Q*b+t==s*l+45*(l-1-rev(a,d)))
    need(Q*a+t==s*l+45*(l-1-rev(b,d)))
    need(delta==s*(b-a)+45*(p(b)-p(a)))
coeff(182684216,184270779,28,T,S,2)
coeff(179,179,8,S,T,0)
L=1<<16
D=(S*L-T)//45+L-1
shifts=[]
loops=[]
for a in range(L//2,L):
    n=D-rev(a,16)
    if n%128: continue
    b=n//128
    if not L//2<=b<L or 128*a+rev(b,16)!=D:continue
    d=S*(b-a)+45*(p(b)-p(a))
    coeff(a,b,16,T,S,d)
    shifts.append(d)
    if a==b:loops.append(a)
need(loops==[44603,44729,44858,44984])
need(len(shifts)==16 and set(shifts)=={-82,-50,-32,-18,0,18,32,50,82})
rng=Random(396730)
for k in range(900):
    x,z=rng.getrandbits(k%270),rng.getrandbits((k*47)%290)
    need(F(x,A*z+V)==F(A*x+U,z)+2)
for k in range(25):
    x,z=rng.getrandbits(k),rng.getrandbits(k+7)
    y=A**k*z+V*((A**k-1)//(A-1))
    top=F(x,y)
    for j in range(k+1):
        need(F(x,y)==top-2*j)
        if j<k:x,y=A*x+U,(y-V)//A
    need(y==z)
need(F(0,V)==189077072195367 and F(U,0)==189077072195365)
print('16 exact width-16 slope matches, 4 self matches')
print('900 full-width unit identities and 325 repeated-transport controls')
def global_matches(d):
    a=1<<d
    DT=(a-1)*(T+45)//45
    DS=(a-1)*(S+45)//45
    ans=[]
    for u in range((DS-a+128)//128,DS//128+1):
        v0=DT-rev(u,d)
        if v0%128:continue
        v=v0//128
        if not a//2<=v<a or 128*u+rev(v,d)!=DS:continue
        ans.append((u,v,S*v-T*u+45*(p(v)-p(u))))
    return ans
small=global_matches(12)
need(small==[(2787,2875,2),(2795,2873,-48),(2803,2874,34),(2811,2872,-16)])
need(U%4096 not in [a for a,_,_ in small] and V%4096 not in [b for _,b,_ in small])
need(128*V+rev(U,36)==(A-1)*(T+45)//45)
need(128*U+rev(V,36)==(A-1)*(S+45)//45)
print("complete width-12 list and primitive width-36 certificate")
