from replay_rich_admission import profile, old_profiles


def starts(rows):
    n=len(rows); f=(1<<n)-1
    d=[[0]*n for _ in range(f+1)]
    h=[1]+[0]*f
    for m in range(1,f+1):
        for i,r in enumerate(rows):
            if not m>>i&1:continue
            a=m^(1<<i)
            d[m][i]=sum(d[a][j] for j in range(n) if (a&r)>>j&1) if a else 1
        h[m]=sum(d[m])
    weights=[sum(d[m][i]*h[f^m] for m in range(f+1)) for i in range(n)]
    return h[f],weights


def literal(rows, mask, s):
    n=len(rows)
    out=[r|(((1<<s)-1) << n)|((not (mask>>i&1)) << (n+s)) for i,r in enumerate(rows)]
    for i in range(s):
        out.append((((1<<s)-(1<<(i+1)))<<n) | (1<<(n+s+1))) # 1's point at later 1's and final 0, not apex
    out.append(mask|(((1<<s)-1)<<n)) # apex -> all 1's
    if s:
        # core must point at last 0; adjust to reserve it after apex.
        for i in range(n):
            out[i]|=1<<(n+s+1)
        out.append(1<<(n+s))
    return out


def main():
    rows=[
        ((24,1,3,6,14),28,75,41),
        ((18,8,3,5,14),21,75,43),
        ((60,45,8,48,6,20),12,75,45),
        ((16,1,3,7,14),27,68,23),
        ((16,1,3,7,14),18,34,25),
        ((16,1,3,7,14),5,34,27),
        ((2,40,43,49,39,1),1,34,29),
        ((16,9,3,5,14),2,34,31),
        ((16,1,3,7,14),8,34,33),
        ((2,40,43,49,39,1),39,136,35),
        ((16,9,3,5,14),7,68,37),
        ((114,8,99,5,110,74,10),58,510,235),
        ((16,93,9,33,12,23,61),83,510,237),
        ((48,9,3,53,70,22,47),32,255,239),
        ((8,61,89,80,1,93,19),96,255,241),
        ((106,88,115,36,41,2,56),16,255,243),
        ((40,57,51,52,97,64,15),8,255,245),
        ((10,16,35,6,45,75,31),29,510,247),
        ((68,9,2,53,71,23,46),8,255,249),
        ((30,64,42,2,14,91,29),32,255,251),
    ]
    for r,m,d,b in rows:
        h,w=starts(r)
        S=sum(x for i,x in enumerate(w) if m>>i&1)
        (hp,q),strong,rich=profile(r,2)
        if (h+S,hp,strong)!=(d,h,True):raise AssertionError((r,h,S,hp,strong))
        if profile(literal(r,m,0),1)[0][0]!=b:raise AssertionError(r)
        for s in (1,2,3,4):
            (v,),stronger,_=profile(literal(r,m,s),1)
            if v!=d*2**s+b or not stronger:raise AssertionError((r,s,v,stronger))
        print('marked h,Q,S,N ; literal checks at s=1..4:',h,q,S,b,'two triangles?',rich)
    if (23,272) in old_profiles(272):raise AssertionError('old clone')
    print('No old ordinal one-apex profile (23,272); bounded factor bank exhausted.')

if __name__=='__main__':main()
