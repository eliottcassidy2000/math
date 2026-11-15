from itertools import product

def profile(rows):
    n=len(rows)
    H=[1]+[0]*((1<<n)-1)
    f=[[0]*n for _ in H]
    for S in range(1,1<<n):
        it=S
        while it:
            b=it&-it
            i=b.bit_length()-1
            T=S^b
            if not T:f[S][i]=1
            else:
                js=rows[i]&T
                while js:
                    c=js&-js
                    f[S][i]+=f[T][c.bit_length()-1]
                    js^=c
            H[S]+=f[S][i]
            it^=b
    return H[-1],sum(H[S]*H[-1^S] for S in range(1<<n))

def trans(n):
    return [(1<<n)-(1<<(i+1)) for i in range(n)]

def glue(w,W,U):
    k=len(w)
    R=trans(k)+[0]
    for i,b in enumerate(w):
        if b=='1':R[k]|=1<<i
        else:R[i]|=1<<k
    j=k+1
    l=j+len(W)
    wmask=((1<<len(W))-1)<<j
    umask=((1<<len(U))-1)<<l
    return [r|wmask|(umask if i<k else 0) for i,r in enumerate(R)]+[(r<<j)|umask for r in W]+[(r<<l)|(1<<k) for r in U]

def value(w):
    return int(w,2) if w else 0

def excess(w):
    return sum((1<<(j-i-1))+(j==i+1) for i,x in enumerate(w) if x=='1' for j in range(i+1,len(w)) if w[j]=='0')

def comp(w):
    return ''.join('0' if c=='1' else '1' for c in w)

def cond(d,t,d2):
    if d<d2:return None
    k=(d//d2).bit_length()-1
    r,rem=divmod(d-(d2<<k),t)
    if rem or r%2 or r>=1<<k:return None
    return '' if not k else comp(f'{r:0{k}b}'[::-1])

def signature(U,W,v):
    h,Q=profile(U);g,R=profile(W)
    d,t,a=(Q-h)*R,h*g,(Q-h)*g
    k=len(v)
    D=d*(1<<k)+t*((1<<k)-1-value(v[::-1]))
    G=d-a+2*t+d*value(v)+t*excess(v)
    return D,t,G,int(v.startswith('0'))

def main():
    cases=[
        (([10,12,9,0],[2,4,1],''),([2,0],[14,12,8,16,7],''),36,252,9),
        ((trans(4),[62,124,56,112,32,64,21],''),([0],[76,109,8,32,47,5,60],'1101'),-758,9660,31),
    ]
    for arm1,arm2,delta,D,t in cases:
        a=signature(*arm1);b=signature(*arm2)
        if a[:2]!=b[:2] or a[:2]!=(D,t) or b[2]-a[2]!=delta or a[3]!=b[3]:raise AssertionError((a,b))
        if cond(a[0],a[1],signature(*arm2[:2],'')[0])!=arm2[-1]:raise AssertionError(arm2)
        actual=[]
        for p in ('1','10','11'):
            got=[]
            for U,W,v in (arm1,arm2):
                H,_=profile(glue(p+v,W,U))
                got.append(H)
            if got[1]-got[0]!=delta or got[0]!=a[2]+D*value(p)+t*excess(p):raise AssertionError((p,got))
            actual.append(tuple(got))
        print('common prefix response',D,t,'signed gap',delta,'literal test pairs',actual)
    num=0
    for k in range(1,12):
        for bits in product('01',repeat=k-1):
            z='1'+''.join(bits)
            r=value(comp(z)[::-1])
            E=sum((1<<(j-i-1))+(j==i+1) for i,x in enumerate(z) if x=='0' for j in range(i+1,k) if z[j]=='0')
            if r-excess(z)!=E:raise AssertionError(z)
            for t,d2 in ((1,2),(3,10),(9,40),(25,200)):
                d=(d2<<k)+t*r
                if cond(d,t,d2)!=z:raise AssertionError((d,t,d2,z))
                num+=1
    print('width reconstruction and deficit checks',num)

if __name__=='__main__':main()
