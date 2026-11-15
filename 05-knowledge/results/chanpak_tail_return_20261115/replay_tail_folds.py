def p(a):
    z=0
    while a:
        q=a//2
        if a%2==0:z+=q+(q%2)
        a=q
    return z

def append(r,m,w):
    r=list(r)
    for c in w:
        n=len(r)
        for k in range(n):r[k]|=1<<n
        r.append(0)
        if c=='1':m|=1<<n
    return tuple(r),m

def apex(r,m):
    n=len(r)
    return tuple(row|(1<<n if not m>>i&1 else 0) for i,row in enumerate(r))+(m,)

def hp(r):
    n=len(r); full=(1<<n)-1
    e=[[0]*n for _ in range(full+1)]
    for s in range(1,full+1):
        a=s
        while a:
            i0=a&-a;i=i0.bit_length()-1;t=s^i0
            b=t&r[i];x=0 if t else 1
            while b:
                k=b&-b;x+=e[t][k.bit_length()-1];b^=k
            e[s][i]=x;a^=i0
    return sum(e[-1])

def main():
    for a in range(4096):
        q=a//2
        if a%2==0: diff=-q-q%2
        else:
            v=a+1;r=(v&-v).bit_length()-1
            w=((v>>r)-1)//2
            diff=(2**r-1)*(2*w+1)-w+1-w%2
            if diff<2:raise AssertionError(a)
        if p(a+1)-p(a)!=diff:raise AssertionError((a,diff))
    print('4096 exact adjacent-address fold checks')
    u0=0
    for k in range(16):
        j=4*k+3;S=2*(17*2**j-15);m=(S-2)//15;a=2*m
        if m%2 or S-15*m!=2 or S+15*(p(a+1)-p(a))!=2:raise AssertionError((j,a))
        if k:
            u=(136*16**(k-1)-1)//15
            if a!=64*u or hex(u)[2:]!='9'+'1'*(k-1):raise AssertionError((k,a,u))
            if p(a)!=p(u)+63*u+1:raise AssertionError((a,u))
            if k>1 and p(u)!=p(u0)+7*u0+1:raise AssertionError(u)
            u0=u
    print('16 preparation heights; all periodic address/correction identities exact')
    for c,j0,expected in ((4,4,36151),(8,5,150401)):
        for k in range(16):
            j=j0+4*k;S=2*(17*2**j-15);m=(S-c)//15;a=2*m
            if m%2 or S-15*m!=c or S+15*(p(a+1)-p(a))!=c:raise AssertionError((j,a))
            hx=('4'+'8'*k+'4') if c==4 else ('8c' if not k else '9'+'1'*(k-1)+'0c')
            if hex(a)[2:]!=hx:raise AssertionError((a,hx))
            base=27+17*2**j+S*a+15*p(a)
            if not k and base!=expected:raise AssertionError((c,base))
            smaller={base+2*i for i in range(4)}
            larger={base+c+2*i for i in range(4)}
            if smaller|larger!=set(range(base,base+c+7,2)):raise AssertionError(c)
        print('16 preparation heights with exact fold step',c,'and first window baseline',expected)
    R=(48,29,25,1,8,30); T=(18,20,17,7,8)
    values=[]
    for i,(r,m,w) in enumerate(((R,25,'1'),(R,16,'11'),(T,1,'11'),(R,1,'11'))):
        for a in (32,33):
            v=''.join('0' if c=='1' else '1' for c in bin(a)[2:][::-1])
            q,mask=append(r,m,w+'0'+v)
            got=hp(apex(q,mask));expected=163+2*i+242*a+15*p(a)
            if got!=expected:raise AssertionError((i,a,got,expected))
            values.append(got)
    if sorted(set(values))!=list(range(8387,8396,2)):raise AssertionError(values)
    print('Eight literal glues cover',sorted(set(values)))
if __name__=='__main__':main()
