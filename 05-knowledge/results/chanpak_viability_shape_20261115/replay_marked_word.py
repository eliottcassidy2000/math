from itertools import product
from collections import Counter

def check(ok,*args):
    if not ok: raise AssertionError(args)

def p(x):
    z=0
    while x:
        x,b=divmod(x,2)
        if not b:z+=x+(x&1)
    return z

def literal(w):
    n=len(w)+1
    dp=[[0]*n for _ in range(1<<n)]
    out=[((1<<(n-1))-1)^((1<<(i+1))-1) for i in range(n-1)]+[0]
    for i,b in enumerate(w):
        if b=='1':out[-1]|=1<<i
        else:out[i]|=1<<(n-1)
    for i in range(n):dp[1<<i][i]=1
    full=(1<<n)-1
    for mask in range(1,full):
        for i,v in enumerate(dp[mask]):
            if not v:continue
            todo=out[i]&~mask
            while todo:
                bit=todo&-todo
                dp[mask|bit][bit.bit_length()-1]+=v
                todo-=bit
    return sum(dp[-1])

marks=Counter(); addresses=Counter(); words=0
for n in range(2,9):
    for ll in product('01',repeat=n-2):
        w='1'+''.join(ll)+'0'
        H=literal(w)
        opposite=''.join('1' if b=='0' else '0' for b in reversed(w))
        check(H==literal(opposite) and H>=1+(1<<(n-2)),w)
        words+=1
        for i in range(n-1):
            if w[i:i+2]!='10':continue
            x=int(w[:i] or '0',2)
            y=sum((b=='0')<<j for j,b in enumerate(w[i+2:]))
            check(H==3+4*x*y+2*x+2*y+p(x)+p(y),w,i)
            marks[H]+=1
    for lx in range(n-1):
        ly=n-2-lx
        xs=range(1<<(lx-1),1<<lx) if lx else (0,)
        ys=range(1<<(ly-1),1<<ly) if ly else (0,)
        for x in xs:
            for y in ys:addresses[3+4*x*y+2*x+2*y+p(x)+p(y)]+=1
check(marks==addresses)
print('literal one-apex words:',words,'exact marked address presentations:',sum(marks.values()))
