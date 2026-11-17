from slot_profiles import profile,append,ordinal
A=(0,41,11,1,15,29);B=(28,17,18,22,0,31)
def tile(e,r):
    a,q=((A,12),(B,43))[e]
    for _ in range(r):
        q|=1<<len(a);a=ordinal(a,(2,0))
    return profile(a,q)
for e in (0,1):
    for r in range(4):
        h,p,w,z,ee,s,k=tile(e,r)
        if (h,p,z,ee)!=(5,128*4**r,(128*4**r-5)//3,0) or w-3*(h-k)!=10+10*r:
            raise ValueError((e,r,h,p,w,z,ee,s,k))
def compose(seq):
    j,W,K,E=5,162,-45,0
    for e,r in seq:
        o=3*W+2*j
        h,p,w,z,ee,s,k=tile(e,0)
        for _ in range(r):
            k-=w;w=4*w+10;z+=p;p*=4;s=5
        j0=j
        j,W,K,E=append((j,W,K,E),(5,p,w,z,0,s,k))
        if 3*W+2*j!=4**r*(128*o+72*(3*e-2)*j0):raise ValueError("recurrence")
    return 4*W+K
def decode(N):
    out=[]
    if N<594:return out
    d=(N//594).bit_length()-1
    if N>=625*2**d:return out
    for m in range(1,d//7+1):
        if (d-7*m)%2:continue
        s=(d-7*m)//2
        J=5**(m+1)
        W0,rem=divmod(3*N-5**m*(27+10*m+10*s),11)
        if rem:continue
        V=3*W0+2*J;left=s;total=J;word=[]
        for _ in range(m):
            if V<=0:break
            val=(V&-V).bit_length()-1
            e=val%2
            r=(val-(3 if e else 4))//2
            if not 0<=r<=left:break
            total//=5;left-=r
            V,z=divmod((V>>2*r)-72*(3*e-2)*total,128)
            if z:break
            word.append((e,r))
        if len(word)==m and left==0 and V==496:out.append(tuple(reversed(word)))
    if N==603:out.append(())
    return out
for seq in (((0,1),),((1,1),(0,2)),((0,0),(1,5),(1,1)),((1,9),)*3):
    x=compose(seq)
    if seq not in decode(x):raise ValueError((x,seq,decode(x)))
    print(seq,x)
