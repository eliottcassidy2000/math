from itertools import permutations,product
from math import gcd, factorial
def need(t,m):
    if not t: raise AssertionError(m)
def spine(w):
    n=len(w); a=[0]*(n+1)
    for i in range(n):
        a[i]=((1<<n)-1)^((1<<(i+1))-1)
        if w[i]=='1':a[n]|=1<<i
        else:a[i]|=1<<n
    return a
def brute(a):
    return sum(all(a[x]>>y&1 for x,y in zip(p,p[1:])) for p in permutations(range(len(a))))
def clock(w):
    h,x,b=1,0,0
    for e in map(int,w):
        h+=(1-e)*(x+b); x,b=2*x+e,e
    return h
def mult(a,mod):
    return len(set(pow(2,k,mod) for k in range(a)))
for n in range(2,8):
 for w in map(''.join,product('01',repeat=n)):
    need(clock(w)==brute(spine(w)),('clock',w))
for a,b,c in product(range(1,7),repeat=3):
    need(clock('1'+'0'*a+'1'*b+'0'*c)==(1+2**b*(2**c-1))*2**a+2+(2**b-1)*(2**c-1),'width')
for j in range(1,10):
    P=2*3**(j-1)
    need(mult(P,3**j)==P,'order')
    need(sum((13*pow(2,a,3**(j+1))+11)%3**j==0 and (13*pow(2,a,3**(j+1))+11)%3**(j+1)!=0 for a in range(3*P))==2,'depth')
def words_at(t):
    out=[]; n=2
    while 1+2**(n-2)<=t:
        for bits in product('01',repeat=n-2):
            w='1'+''.join(bits)+'0'
            if clock(w)==t: out.append(w)
        n+=1
    return out
for j in range(9):
    ws=words_at(7*3**j)
    need(len(ws)==(2 if j in (2,8) else 0),('ray',j))
    print(j,7*3**j,ws)
def s(m,k):
    from math import comb
    return sum((-1)**(k-i)*comb(k,i)*i**m for i in range(1,k+1))
for m in range(1,8):
    # substitute TT_m for a vertex of C3. a1=3,a2=1,a>=3=0
    need(clock('1'*m+'0')==3*s(m,1)+s(m,2)==1+2**m,'single module')
print('PASS independent permutation counts and arithmetic controls')
def dp(a):
    n=len(a); d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for m,row in enumerate(d):
        for i,x in enumerate(row):
            if not x: continue
            s=a[i]&~m
            while s:
                b=s&-s;s-=b;d[m|b][b.bit_length()-1]+=x
    return sum(d[-1])
def inflate(blocks):
    rows=[]; o=0; n=sum(map(len,blocks))
    for b in blocks:
        z=((1<<n)-1)^((1<<(o+len(b)))-1)
        rows.extend(z|(r<<o) for r in b); o+=len(b)
    # apex has not yet been included: bits for blocks 0 and 2 are ones
    h=[i for i in range(len(blocks)) if i in (0,2)]
    beats=0;o=0
    for i,b in enumerate(blocks):
        if i in h: beats|=((1<<len(b))-1)<<o
        else:
            for j in range(len(b)): rows[o+j]|=1<<n
        o+=len(b)
    return rows+[beats]
T2=[2,0];C3=[2,4,1]
for s in product(range(2),repeat=3):
    core=inflate([[0],*[(T2,C3)[i] for i in s]])
    need(dp(core)==63*3**sum(s),('subset upgrade',s))
print('PASS 8 apex-block subset upgrades via subset DP')
