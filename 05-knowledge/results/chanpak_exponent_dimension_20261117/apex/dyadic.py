import argparse
from gap import F, count, pair, physical

def lower(a):
    return count(*pair(a,a+1)[0])

def decode(k, N):
    if k<2 or N%2==0: raise ValueError('k>=2 and odd target')
    K=max(6,k)
    r=1 if N%4==3 else 9
    def address(c):
        x=int("1"*K + "".join(
            "0"*(j+((c>>(j-1))&1)) + "1"*(K+1-((c>>(j-1))&1))
            for j in range(1,k-1)),2)
        return (x<<K)+r
    a0=address(0)
    if (lower(a0)-N)%4: raise ArithmeticError('wrong phase')
    c=((lower(a0)-N)//4)%(1<<(k-2))
    a=address(c)
    P=pair(a,a+1)
    L,U=(count(*p) for p in P)
    if (L-N)%(1<<k) or U!=L+2: raise ArithmeticError('wrong read')
    return a,L,U,P

if __name__=='__main__':
    tests=0
    for k in range(2,12):
        orders=set()
        for N in range(1,1<<k,2):
            a,L,U,P=decode(k,N)
            orders.add(len(physical(*P[0])))
            tests+=1
        if len(orders)!=1: raise ArithmeticError('nonuniform order')
    print(f'{tests} exact dyadic phases with a physical +2 mate (k=2..11)')
