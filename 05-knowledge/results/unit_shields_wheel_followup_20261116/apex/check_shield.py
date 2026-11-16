from itertools import product

def gap(w):
    a,b=zip(*(divmod(z,2) for z in w))
    def F(j,i):return (1-a[i])*b[j]-(1-b[i])*a[j]
    return (b[0]-a[0]+b[-1]-a[-1]
            +sum((1<<(i-j-1))*F(j,i) for i in range(len(w)) for j in range(i))
            +sum(F(i+1,i) for i in range(len(w)-1)))

def rows(w,ab):
    n=len(w)
    A=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[0,0]
    A[n if ab else n+1]|=1<<(n+1 if ab else n)
    for i,z in enumerate(w):
        a,b=divmod(z,2)
        for mark,j in ((a,n),(b,n+1)):
            if mark:A[j]|=1<<i
            else:A[i]|=1<<j
    return A

def literal(w,ab):
    A=rows(w,ab)
    from functools import cache
    @cache
    def f(S,j):
        if S==1<<j:return 1
        rem=S^(1<<j)
        q=rem&A[j]
        ans=0
        while q:
            b=q&-q
            ans+=f(rem,b.bit_length()-1)
            q^=b
        return ans
    return sum(f((1<<len(A))-1,j) for j in range(len(A)))

def fast(w):
    def transition(old,a,b):
        v=[0]*12
        for f,l,m in product(range(2),range(2),range(3)):
            x=old[6*f+3*l+m]
            v[6*a+3*l+m]+=x
            if a or m:v[6*f+3*l+1+b]+=x
            if b or l:v[6*f+3+m]+=x
        return v
    x=y=[1]+[0]*11
    for z in w:
        a,b=divmod(z,2)
        x=transition(x,a,b)
        y=transition(y,b,a)
    return x[0]+x[3]+x[1]+x[4]+y[1]+y[4], y[0]+y[3]+y[1]+y[4]+x[1]+x[4]

def check():
    for n in range(1,7):
        for w in product(range(4),repeat=n):
            A,B=fast(w)
            assert A-B==gap(w)
    print('Checked direct contraction gap against three-run transfer on all 5,460 words through spine 6.')
    for t in range(1,19):
        w=(1,2)+(0,)*t+(2,)
        A,B=fast(w)
        assert A==4*3**t+9*2**t-3 and B-A==2 and gap(w)==-2
        if t<=5:assert (literal(w,True),literal(w,False))==(A,B)
    print('Confirmed single-arc factor-one seam and baseline formula through height 18, five literal DP controls.')
    print('Mirror width, individual anchors, distinct larger endpoints, min, max')
    for k in range(2,14):
        out=set()
        for middle in product((0,3),repeat=k-2):
            u=(3,)+middle+(3,)
            w=u+(0,1,2,0)+tuple(3-z for z in u[::-1])
            A,B=fast(w)
            assert A-B==2 and gap(w)==2
            assert B>=1+2**(2*k+2)
            if k<=4:assert (literal(w,True),literal(w,False))==(A,B)
            out.add(A)
        print(k,2**(k-2),len(out),min(out),max(out))

    pp=tuple(map(int,'133220032220'))
    qq=tuple(3-z for z in pp[::-1])
    assert qq==tuple(map(int,'311103311002'))
    assert (literal(pp,True),literal(pp,False))==(186399,186397)
    assert (literal(qq,True),literal(qq,False))==(186397,186399)
    for n in range(1,9):
        for bits in product((0,1),repeat=n):
            w=sum((pp if b else qq for b in bits),())
            g=2*(2*sum(bits)-n)
            A,B=fast(w)
            assert gap(w)==g and A-B==g
    print('The balanced P/Q atoms were checked by literal subset DP; every P/Q block sequence through 8 atoms obeys the signed gap law.')

if __name__=='__main__':check()
