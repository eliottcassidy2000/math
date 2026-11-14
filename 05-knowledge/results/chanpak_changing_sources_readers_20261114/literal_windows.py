BANKS = {(7, 135, -1): [(-1,
                 (7, 0, 41, 58, 40, 78, 40, 2, 4, 15, 29, 62),
                 (2, 1, 4, 0, 1),
                 1),
                (1,
                 (6, 5, 41, 56, 42, 48, 33, 3, 7, 14, 28),
                 (2, 1, 4, 0, 1),
                 1),
                (3, (3, 0, 1, 4, 4, 0, 1, 3), (4, 5, 32, 8, 1, 3, 6), 0),
                (5,
                 (7, 5, 37, 64, 42, 100, 1, 2, 7, 15, 30, 62),
                 (2, 1, 4, 0, 1),
                 1),
                (7,
                 (6, 3, 37, 62, 44, 50, 36, 1, 7, 14, 28),
                 (2, 1, 4, 0, 1),
                 1),
                (9,
                 (7, 0, 57, 66, 60, 106, 92, 1, 4, 13, 30, 60),
                 (3, 3, 12, 4, 1, 2),
                 3),
                (11,
                 (7, 5, 35, 64, 46, 102, 8, 2, 5, 15, 30, 62),
                 (2, 1, 4, 0, 1),
                 1)],
 (9, 27, -5): [(-5,
                (7, 6, 133, 182, 128, 96, 81, 3, 7, 13, 30, 60),
                (2, 1, 4, 0, 1),
                5),
               (-3, (5, 4, 9, 12, 6, 20, 1, 2, 7, 14), (2, 1, 4, 0, 1), 1),
               (-1,
                (7, 6, 131, 184, 130, 112, 69, 17, 7, 10, 30, 60),
                (2, 1, 4, 0, 1),
                5),
               (1,
                (7, 6, 65, 92, 66, 108, 65, 2, 6, 15, 30, 60),
                (2, 1, 4, 0, 1),
                4),
               (3,
                (7, 6, 131, 180, 134, 100, 65, 18, 7, 11, 30, 60),
                (2, 1, 4, 0, 1),
                5),
               (5,
                (7, 0, 67, 82, 72, 100, 17, 2, 7, 13, 30, 62),
                (2, 1, 4, 0, 1),
                4),
               (7, (4, 0, 1, 8, 8, 0, 1, 3, 7), (2, 1, 4, 0, 1), 0),
               (9,
                (7, 0, 129, 180, 138, 86, 68, 40, 3, 14, 27, 60),
                (2, 1, 4, 0, 1),
                5),
               (11,
                (7, 5, 61, 94, 72, 92, 37, 16, 6, 10, 29, 62),
                (2, 1, 4, 0, 1),
                4)],
 (10, 13, 189): [(189,
                  (7, 6, 91, 136, 98, 104, 77, 17, 4, 11, 30, 60),
                  (1, 1, 2, 0),
                  4),
                 (191,
                  (7, 1, 91, 134, 100, 104, 77, 17, 4, 11, 30, 60),
                  (1, 1, 2, 0),
                  4),
                 (193,
                  (7, 6, 95, 128, 98, 102, 96, 2, 7, 15, 28, 60),
                  (1, 1, 2, 0),
                  4),
                 (195,
                  (7, 5, 81, 140, 114, 100, 41, 2, 5, 15, 28, 62),
                  (1, 1, 2, 0),
                  4),
                 (197,
                  (7, 0, 89, 130, 108, 104, 89, 3, 4, 13, 30, 60),
                  (1, 1, 2, 0),
                  4),
                 (199,
                  (7, 3, 75, 142, 124, 102, 48, 2, 7, 13, 28, 62),
                  (1, 1, 2, 0),
                  4),
                 (201,
                  (7, 1, 75, 140, 126, 88, 37, 1, 6, 14, 29, 62),
                  (1, 1, 2, 0),
                  4),
                 (203,
                  (7, 3, 65, 148, 138, 98, 20, 1, 7, 13, 30, 62),
                  (1, 1, 2, 0),
                  4),
                 (205,
                  (7, 3, 47, 164, 158, 120, 1, 3, 6, 14, 30, 62),
                  (1, 1, 2, 0),
                  4),
                 (207,
                  (7, 1, 35, 174, 172, 102, 8, 2, 5, 15, 30, 62),
                  (1, 1, 2, 0),
                  4)]}


def tt(n):
    return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))


def profiles(r):
    n=len(r)
    F=(1<<n)-1
    E=[[0]*n for _ in range(F+1)]
    U=[1]+[0]*F
    for i in range(n):
        E[1<<i][i]=1
    for m in range(1,F+1):
        U[m]=sum(E[m])
        for v,t in enumerate(E[m]):
            if t:
                avail=r[v]&~m
                while avail:
                    b=avail&-avail
                    avail^=b
                    E[m|b][b.bit_length()-1]+=t
    return U,E


def hq(r):
    U,_=profiles(r)
    F=len(U)-1
    return U[F],sum(h*U[F^m] for m,h in enumerate(U))


def reader(r,q):
    U,E=profiles(r)
    F=len(U)-1
    e=o=0
    for m in range(1,len(U)):
        p=E[m][q]*U[F^m]
        if m.bit_count()%2:o+=p
        else:e+=p
    return U[F],2*o,2*e


def cycle(a,b):
    M,N=len(a),len(b)
    BM=((1<<N)-1)<<(M+1)
    return (((1<<M)-1)<<1,)+tuple((row<<1)|BM for row in a)+tuple((row<<(M+1))|1 for row in b)


def glue(t,r,q):
    m,n=len(t),len(r)
    R=((1<<n)-1)<<m
    s=[row|R for row in t]+[row<<m for row in r]
    s[0]^=1<<(m+q)
    s[m+q]|=1
    return tuple(s)


def check():
    ngraphs=0
    for (length,f,first),records in BANKS.items():
        assert [o for o,*_ in records]==list(range(first,first+2*length,2))
        for c,r,a,v in records:
            size,q,h,oo,ee,*rr=r
            asize,x,Q,*aa=a
            assert len(rr)==size and len(aa)==asize
            assert hq(aa)==(x,Q)
            assert reader(rr,q)==(h,oo,ee)
            assert profiles(converse(rr))[1][-1][q]==(oo-ee)//2
            D=Q-x
            assert (h+(oo+ee)//2)*D==f*(1<<v)
            assert (2*h+ee)*x-h*D==c
            k0=1+max(z[-1] for z in records)
            for k in (k0,k0+1):
                m=k-v
                G=glue(cycle(aa,tt(m)),rr,q)
                assert profiles(G)[0][-1]==f*(1<<k)+c
                ngraphs+=1
    print(f'Checked {ngraphs} finished graphs by independent literal Hamiltonian DP.')



FOUR = (
    (6, (8,5,9,82,39,15,183,63)),
    (6, (2,8,3,53,39,71,159,63)),
    (6, (2,4,9,51,39,71,159,63)),
    (6, (2,12,9,81,7,31,183,63)),
)

def converse(r):
    F=(1<<len(r))-1
    return tuple(F^row^(1<<i) for i,row in enumerate(r))

def ordinal(a,b):
    m=len(a)
    R=((1<<len(b))-1)<<m
    return tuple(r|R for r in a)+tuple(r<<m for r in b)

def check_extra():
    ngraphs=0
    for e,(q,R) in enumerate(FOUR,27):
        h,a,b=reader(R,q)
        assert (h,h+(a+b)//2,b//2)==(73,146,e)
        assert profiles(converse(R))[1][-1][q]==(a-b)//2
        for k,m in ((1,1),(2,1),(1,2),(2,2)):
            assert profiles(glue(cycle(tt(k),tt(m)),R,q))[0][-1]==(2**k-1)*(146*2**m-73)+2*(73+e)
            ngraphs+=1
    C3=(2,4,1)
    for k,m in ((2,1),(2,2),(3,1),(3,2)):
        hi=glue(cycle(tt(k),tt(m)),C3,0)
        lo=glue(cycle(ordinal(tt(m-1),C3),tt(k-1)),tt(1),0)
        assert tuple(profiles(G)[0][-1] for G in (lo,hi))==tuple(3*(2**k-1)*(2**(m+1)-1)+c for c in (6,8))
        ngraphs+=2
    print(f'Checked {ngraphs} finished 2-axis graphs, including literal source/reader swap.')


if __name__=='__main__':
    check()
    check_extra()
