"""Independent DP checks for BALANCE.md."""
from functools import lru_cache
from itertools import product

def K(d):
    return (1 << (d-1)) + (d == 1)

def formula(N,V):
    return 1+2**N+sum(2**i-2**(N-1-i) for i in V)-sum(K(i-j) for i in V for j in V if i>j)

def literal(N,V):
    rows=[sum(1<<j for j in range(i+1,N+1)) + ((i in V or i==N) << (N+1)) for i in range(N+1)]
    rows.append(sum(1<<i for i in range(N) if i not in V))
    @lru_cache(None)
    def visit(i,left):
        if not left:
            return 1
        nxt=left&rows[i]
        out=0
        while nxt:
            b=nxt&-nxt
            nxt-=b
            out+=visit(b.bit_length()-1,left-b)
        return out
    all_v=(1<<len(rows))-1
    return sum(visit(i,all_v-(1<<i)) for i in range(len(rows)))

def need(b,msg):
    if not b:
        raise AssertionError(msg)

count=0
for a,b in product((1,2), repeat=2):
    for n in range(1,4):
        N=a+b+2*n+1
        m=a+n
        for choices in product((0,1),repeat=2*n+1):
            V={a+i for i,c in enumerate(choices) if c}
            need(formula(N,V)==literal(N,V),(a,b,n,V))
            count+=1
        A={2} if n>=2 else set()
        A=A|{0}
        B={1}|({3} if n>=3 else set())
        for orient in product((-1,1),repeat=len(B)):
            for on in product((0,1),repeat=len(A)):
                U={m+sgn*r for r,sgn in zip(sorted(B),orient)}
                Up={m+r for r in B}
                V={m+sgn*r for r,occ in zip(sorted(A),on) if occ for sgn in (-1,1)}
                need(formula(N,U|V)-formula(N,U)==formula(N,Up|V)-formula(N,Up),(a,b,n,U,V))
for n in range(2,80):
    for r in range(1,min(3,n)):
        c=n+2
        N=2*n+4
        for sign in (-1,1):
            full=set(range(2,2*n+3))
            V=full-{c+sign*r}
            exponent_diff=2**(n+r)-2**(n-r)
            H0=9*2**(2*n+1)+1
            need(formula(N,V)==H0-sign*3*exponent_diff,(n,r,sign))
            need(formula(N,V-{c})==H0-sign*3*exponent_diff-2,(n,r,sign))
print(f'{count} literal head/tape graphs; exact rectangles and saturation scaffolds pass')

for R in range(2,12):
    N=2*R+4
    i=R+2
    count_pos=count_neg=0
    last=None
    for y in range(2**(R-1)):
        V={i-R-1}|{i+(d if (y>>(d-2))&1 else -d) for d in range(2,R+1)}
        h=formula(N,V)
        need(formula(N,V|{i})==h+2,(R,y,V,h))
        if last is not None:
            need(h-last>2**R,(R,y,h,last))
        last=h
        for u,v in product((0,1),repeat=2):
            for sign in (-1,1):
                x=2**R-sign-y-u-v
                if not 0<=x<2**R:
                    continue
                W={i-d for d in range(2,R+2) if (x>>(d-2))&1}
                W|={i+d for d in range(2,R+1) if (y>>(d-2))&1}
                W|=({i-1} if u else set())|({i+1} if v else set())
                need(formula(N,W|{i})-formula(N,W)==2*sign,(R,x,y,u,v,W))
                count_pos+=sign==1
                count_neg+=sign==-1
    need(count_pos==2**(R+1) and count_neg==count_pos-4,(R,count_pos,count_neg))
print('pivot-unit exact mask multiplicities and ordered absolute anchors pass')
