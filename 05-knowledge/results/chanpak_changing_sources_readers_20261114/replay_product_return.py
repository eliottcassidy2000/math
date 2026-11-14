from functools import cache

def check(x):
    if not x:
        raise AssertionError("product return control failed")

readers = [
 (8,5,9,82,39,15,183,63),
 (2,8,3,53,39,71,159,63),
 (2,4,9,51,39,71,159,63),
 (2,12,9,81,7,31,183,63)]

def tt(k):
    return tuple((1<<k)-(2<<i) for i in range(k))

def join(a,b):
    return tuple(x|(((1<<len(b))-1)<<len(a)) for x in a)+tuple(x<<len(a) for x in b)

def cyclic(a,b):
    k=len(a)
    return ((1<<(k+1))-2,)+tuple(x<<1|(((1<<len(b))-1)<<(k+1)) for x in a)+tuple(
           (x<<(k+1))|1 for x in b)

def seam(a,r,b,w):
    out=list(join(a,b))
    u=len(a)+w
    out[r]^=1<<u
    out[u]^=1<<r
    return tuple(out)

def profile(rows):
    n=len(rows)
    full=(1<<n)-1
    back=[sum(1<<i for i in range(n) if (rows[i]>>j)&1) for j in range(n)]
    check(all((rows[i]&back[i])==0 and (rows[i]|back[i])==(full^(1<<i)) for i in range(n)))

    @cache
    def path(s,v,rev):
        z=s^(1<<v)
        if z==0: return 1
        dest=(back if rev else rows)[v]&z
        q=0
        while dest:
            w=(dest&-dest).bit_length()-1
            dest^=1<<w
            q+=path(z,w,rev)
        return q

    @cache
    def h(s):
        return sum(path(s,v,False) for v in range(n) if (s>>v)&1) if s else 1

    def marked(v,rev):
        e=o=0
        for s in range(1,full+1):
            if (s>>v)&1:
                q=path(s,v,rev)*h(full^s)
                if s.bit_count()%2: o+=q
                else: e+=q
        return e,o

    return h(full), marked, (lambda v,rev: path(full,v,rev)), (lambda: sum(
        h(s)*h(full^s) for s in range(full+1)))

for i,R in enumerate(readers):
    a,marks,ends,Q=profile(R)
    qe,qo=marks(6,True)
    check((a,qe,qo)==(73,27+i,46-i))
    for k,m in [(1,1),(2,2),(3,1)]:
        T=seam(cyclic(tt(k),tt(m)),0,R,6)
        h,mt,ep,_=profile(T)
        e,o=mt(0,False)
        check(h==73*(2**k-1)*(2**(m+1)-1)+200+2*i)
        check((ep(0,False),ep(0,True))==(a,ends(6,True)))
        check(e+o==(1+(2**k-1)*2**m)*Q()-sum(marks(6,False)))
        if i==0 and k+m==4:
            expanded=seam(T,0,tt(2),1)
            check(profile(expanded)[0]==h+2*(e+o))
print("Four literal two-axis windows; endpoint equality and rooted budget reanchoring; physical continuation checked.")
