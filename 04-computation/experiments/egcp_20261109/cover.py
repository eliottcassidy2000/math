# coverage closure under operations with proved evenness; bitsets as python ints (bit k = value k)
"""Independent closure check for THM-5040: closes the exact census sets (n <= 9) under padding, cones, E_2/K_1/E_3
multi-joins, two-hub graphs, the top band, and the leaf stable recursion, and compares E(n) with O(n).
Two-hub graph TH(G, t): non-adjacent hubs a, b both joined to all of V(G) (f = |G| vertices), t pendant leaves at a;
d = 2^(n-1) + 2^f - 2^t + d(G) with n = f + t + 2; even iff G is even, except (f, t) = (1, 1) where TH = P_4 is odd
(never used here, n >= 10).  Corroboration only: THM-5040 is proved without it (audit Note 8).
usage: python cover.py NMAX"""
import json, sys
_js=json.load(open('base_sets.json'))
E_b={int(n):set(v) for n,v in _js['E'].items()}; E0_b={int(n):set(v) for n,v in _js['E0'].items()}
C_b={int(n):set(v) for n,v in _js['C'].items()}
NB=9; NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 20
USE_TOP = '--notop' not in sys.argv
def bits(s):
    x=0
    for k in s: x|=1<<k
    return x
E={n:bits(E_b[n]) for n in range(1,NB+1)}
E0={n:bits(E0_b[n]) for n in range(1,NB+1)}
C={n:bits(C_b[n]) for n in range(1,NB+1)}
E0[1]=0; C[1]=0   # K1 is its own universal vertex
Lb={int(n):set(v) for n,v in _js['L'].items()}
L={n:bits(Lb[n]) for n in range(1,NB+1)}
USE_LEAF='--noleaf' not in sys.argv
def Mx(n): return 2**n-n if n%2 else 2**n-n-1
def M0x(n): return 2**n-n-1 if n%2==0 else 2**n-n-4
def oddmask(top):
    # bits at all odd k with 1 <= k <= top
    m=(top+1)//2          # number of odd values <= top
    return ((4**m-1)//3)<<1
def target_E(n):
    s=oddmask(Mx(n))
    for k in (3,7,13): s&=~(1<<k)
    if n%2 and n>=5: s&=~(1<<(2**n-n-2))
    return s
def target_E0(n):
    s=oddmask(M0x(n))
    for k in (3,7,13,M0x(n)-2): s&=~(1<<k)
    return s
def topband(n):
    N=n-1; s=0
    J=list(range(1,N,2)) if N%2==0 else list(range(1,N-1,2))+[N+2]
    for j in J: s|=1<<(2**N-j)
    return s
def listbits(x,limit=20):
    out=[];k=0
    while x and len(out)<limit:
        low=(x&-x).bit_length()-1; out.append(low); x&=x-1
    return out
for n in range(NB+1,NMAX+1):
    e=E[n-1]; e0=E[n-1]; c=E[n-1]
    e|=E0[n-1]<<(2**(n-1))                       # cone
    for g in range(1,n):
        r=n-g
        if r%2==0:                                # E2^t v G
            t=r//2; base=2**n-2**g-2*t
            e|=E[g]<<base; e0|=E0[g]<<base
        else:                                     # E2^t v K1 v G
            t=(r-1)//2; base=2**n-2**g-2*t
            e|=E0[g]<<base
        if r>=3 and (r-3)%2==0:                   # E2^t v E3 v G, G co-connected, not E3
            t=(r-3)//2; base=2**n-2**g-2*t-6
            cg=C[g]
            if g==3: cg&=~(1<<1)
            e0|=cg<<base; e|=cg<<base
    for t in range(1,n-2):                        # two-hub
        f=n-2-t; base=2**(n-1)+2**f-2**t
        e0|=E[f]<<base; c|=E[f]<<base
    if USE_TOP:
        tb=topband(n); e|=tb; e0|=tb; c|=tb
    l=L[n-1]
    for t in range(1,n-2):
        f=n-2-t; base=2**(n-1)+2**f-2**t
        l|=E[f]<<base
    if USE_TOP: l|=topband(n)
    if USE_LEAF and n-2>=1: l|=L[n-2]<<(3*2**(n-3))
    L[n]=l; c|=l; e0|=l
    E[n]=e|e0|c; E0[n]=e0|c; C[n]=c
    missE=target_E(n)&~E[n]; missE0=target_E0(n)&~E0[n]
    extraE=E[n]&~target_E(n)
    print(n,"missing E:",bin(missE).count('1'),listbits(missE,12),"| missing E0:",bin(missE0).count('1'),listbits(missE0,12),"| E beyond target:",listbits(extraE,5),flush=True)

def show(n):
    missE=target_E(n)&~E[n]; missE0=target_E0(n)&~E0[n]
    def fmt(v):
        # express v relative to powers of two
        p=v.bit_length()-1; r=v-2**p
        q=r.bit_length()-1 if r>0 else -1
        return "%d=2^%d+%d"%(v,p,r) if r<64 or q<0 else "%d=2^%d+2^%d+%d"%(v,p,q,r-2**q)
    print("n=%d missing E:"%n, [fmt(v) for v in listbits(missE,200)])
    print("n=%d missing E0 (not in missing E):"%n, [fmt(v) for v in listbits(missE0&~missE,200)])
if '--show' in sys.argv:
    show(NMAX)
if '--check' in sys.argv:
    for v in [61439, 2**16-2**13+2**12-3, 2**15-2**12+2**11-1, 2**14 - 2**11 + 2**10 - 1]:
        n=v.bit_length()
        print(v, "n=",n, "in E[n]:", bool(E[n]>>v&1), "in E0[n]:", bool(E0[n]>>v&1))
    print("missing E at n=16:", [ (v, bin(v)) for v in listbits(target_E(16)&~E[16],200)][-30:])
if '--leafoffsets' in sys.argv:
    for n in range(8,NMAX+1):
        x=L[n]>>(2**(n-1))
        offs=listbits(x,10**6)
        offs=[c for c in offs if c< 2**(n-1)]
        # missing odd offsets below 2^(n-2)+64
        lim=2**(n-2)
        miss=[c for c in range(1,lim,2) if not (x>>c)&1]
        print(n,"L upper offsets: max",max(offs) if offs else None,"missing odd c <",lim,":",len(miss),miss[:15], "...", miss[-5:])
