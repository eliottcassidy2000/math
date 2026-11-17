from check import need, apply, steer
from itertools import product

def F(x):
    return apply(bin(x)[2:])[1]

def H(u,w):
    return 2*u*w+u+(u%2)+F(u)+F(w)-1

def literal(u,w):
    bits=bin(u)[2:]+'0'+''.join(str(1-int(c)) for c in bin(w)[:1:-1])
    n=len(bits)
    a=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)] + [sum(int(c)<<i for i,c in enumerate(bits))]
    for i,c in enumerate(bits):
        if c=='0':
            a[i]|=1<<n
    # Independent literal Hamilton-path counting via required consecutive arcs.
    layer={(1<<i,i):1 for i in range(n+1)}
    for _ in range(n):
        ext={}
        for (mask,j),v in layer.items():
            pending=a[j]&~mask
            while pending:
                b=pending&-pending; pending-=b
                k=(mask|b,b.bit_length()-1)
                ext[k]=ext.get(k,0)+v
        layer=ext
    return sum(layer.values())

for u,w in product(range(1,14),repeat=2):
    need(H(u,w)==literal(u,w))
for D in (2,4,6,8,14,18,24,36,68,100):
    t=2+(D-1).bit_length()
    r=2**(t-2)-D
    s=r+D
    R=2**t-1-int(format(r,'0'+str(t)+'b')[::-1],2)
    S=2**t-1-int(format(s,'0'+str(t)+'b')[::-1],2)
    def J(x):
        u,w=(x<<t)+r,(x<<t)+s
        return 2*u*w+u+(R+S)*x+F(r)+F(s)-3
    for x in range(1,50):
        u,w=(x<<t)+r,(x<<t)+s
        need(H(w,u)-H(u,w)==D and H(u,w)==2*F(x)+J(x)+2*(x%2))
    for M,N in product((8,15,24,35,72,144),(7,21,30011)):
        eps=((N-J(0)-2)//2)%2
        q=next(z for z in range(1,2*M,2) if (2*z+J(eps)+2*eps-N)%M==0)
        word=steer(M,'1011',str(eps),eps,q)
        x,z,_=apply(word,M=M)
        need((2*z+J(x)+2*eps-N)%M==0)
print('169 literal arc-DPs; exact gap identities and 180 modular witnesses pass.')
