from itertools import product

def check(x):
    if not x:raise AssertionError('height envelope')
def H(w):
    p=s=0;h=1
    for v in w:
        if not v:h+=p+s
        p,s=2*p+v,v
    return h
for n in range(2,15):
    full=[H((1,)+w+(0,)) for w in product((0,1),repeat=n-2)]
    q=n-2
    aa=[2**j+(j==0) for j in range(q)]
    lo=3 if q==0 else 1+2**q+sum(min(i,j) for i,j in zip(aa,aa[::-1]))
    hi=2**n-2**(n//2)-2**((n+1)//2)+3
    check((min(full),max(full))==(lo,hi))
    if q>=2:
        check(lo==1+2**q+(2**(q//2+1) if q%2==0 else 3*2**(q//2)))
print('Exact active-width envelopes through length fourteen.')
