from random import Random
from verify import H,ck
rng=Random(518)
trials=0
for _ in range(2000):
    n=rng.randrange(3,60)
    r=rng.randrange(1<<n)|1
    s=rng.randrange(1<<n)&~1
    split=(r^s)&~1
    A=split&rng.randrange(1<<n)
    C=split&~A&rng.randrange(1<<n)
    if not A or not C: continue
    trials+=1
    e=[(r>>i&1)-(s>>i&1) for i in range(n)]
    a=sum(e[i]<<i for i in range(n) if A>>i&1)
    b=sum(e[i]<<i for i in range(n) if C>>i&1)
    j=sum(e[i]*e[k]*((1<<(abs(i-k)-1))+(abs(i-k)==1))
          for i in range(n) for k in range(n) if A>>i&1 and C>>k&1)
    x,y=r+(7<<n),s+(5<<n)
    d=H(x^(A|C),y^(A|C))-H(x^A,y^A)-H(x^C,y^C)+H(x,y)
    ck(d==-2*(256*a*b+5*j))
    ck(abs(j)<=2*abs(a*b) and abs(d)>=492*abs(a*b))
print(f'Mixed physical digit differences: {trials} signed controls pass.')
