from itertools import product

def require(ok, detail):
    if not ok: raise AssertionError(detail)

def F(u):
    s = 1
    while u:
        q = u//2
        if not u%2: s += q + q%2
        u = q
    return s

def one_apex(w):
    return (1 if not w else
        int(w[0]=='1')+int(w[-1]=='0')+
        w.count('01')+
        sum(1<<(j-i-1) for j,b in enumerate(w) for i in range(j)
            if w[i]=='1' and b=='0'))

def ham(a):
    n = len(a)
    require(all((a[i]>>j&1)+(a[j]>>i&1)==1 for i in range(n) for j in range(i)),a)
    d = [[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for s,r in enumerate(d):
        for i,v in enumerate(r):
            if v:
                q=a[i]&~s
                while q:
                    b=q&-q
                    q -= b
                    d[s|b][b.bit_length()-1]+=v
    return sum(d[-1])

def blow(a,mark,m):
    lab=[i for i in range(len(a)) for _ in range(m if i==mark else 1)]
    return [sum(1<<j for j,y in enumerate(lab) if (i<j if x==y else a[x]>>y&1))
            for i,x in enumerate(lab)]

def encode(w):
    n=len(w)
    return [sum(1<<j for j in range(i+1,n))+
            ((1<<n) if b=='0' else 0) for i,b in enumerate(w)]+[
            sum(1<<i for i,b in enumerate(w) if b=='1')]

for k in range(1,8):
    for w in map(''.join,product('01',repeat=k)):
        require(F(int(w,2))==one_apex(w),w)
        require(ham(encode(w))==one_apex(w),w)
        for i in range(k):
            if w[i]!='0' or '1' not in w[:i]:continue
            a=int(w[:i],2)
            b=int(''.join('1' if t=='0' else '0' for t in w[i+1:][::-1]) or '0',2)
            c=a*(b+1)
            d=F(a)+F(b)-a-(a%2==0)
            for m in (1,2,4):
                require(one_apex(w[:i]+'0'*m+w[i+1:])==c*2**m+d,(w,i,m))
        for i in range(k):
            for j in range(i+1,k):
                x,y=int(w[i]),int(w[j]); inner=w[i+1:j]
                L=2**len(inner)
                A=int(w[:i] or '0',2)
                B=int(''.join('1' if t=='0' else '0' for t in w[j+1:][::-1]) or '0',2)
                d=int(inner or '0',2)
                z=int(''.join('1' if t=='0' else '0' for t in inner[::-1]) or '0',2)
                cross=L*(A+x)*(B+1-y)
                # For both 00 and 11 take the separately derived coefficients.
                if (x,y)==(0,0): U,V=A*(1+z-L),(B+1)*d
                if (x,y)==(1,1): U,V=(A+1)*z,B*(d+1-L)
                if (x,y)==(1,0): U,V=(A+1)*(z-L),(B+1)*(d-L)
                if (x,y)==(0,1): U,V=A*(1+z),B*(d+1)
                c0=one_apex(w)-4*cross-2*U-2*V
                for s,m in [(1,2),(2,1),(2,2),(3,4)]:
                    W=w[:i]+w[i]*s+inner+w[j]*m+w[j+1:]
                    require(one_apex(W)==c0+U*2**s+V*2**m+cross*2**(s+m),
                            (w,i,j,s,m,cross,U,V))

for a,r,h,c in [([64,5,49,39,75,403,174,31,223],8,737,368),
                ([254,124,120,48,224,192,8,78],6,59,29),
                ([254,204,248,112,98,194,128,24],6,75,37)]:
    require(min(a[r].bit_count(),len(a)-1-a[r].bit_count())==1,a)
    require(ham(a)==h and ham(blow(a,r,2))==h+2*c,a)
def rooted_word(w):
    out = [2,0]
    for b in map(int,w):
        n=len(out)
        out=[u|((1<<n) if i or not b else 0)
             for i,u in enumerate(out)]+[b]
    return out

for j, edges in enumerate([[(0,2),(0,5)],[(0,6),(2,4)],
                            [(0,3),(5,7)],[(0,1),(3,7),(5,9)]]):
    a=rooted_word('1000'+'1'*(j+1))
    require(ham(a)==(46<<j)+1-20,(j,a))
    for u,v in edges:
        a[u]^=1<<v
        a[v]^=1<<u
    require(ham(a)==(46<<j)+1,(j,a))
print("PASS: one- and two-block formulas, three exact rays, and four finite patches")
