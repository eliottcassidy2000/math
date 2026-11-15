"""Independent subset-DP receipts for the cross-edit formula and exact decoder."""
from itertools import product, combinations
from random import Random

def tournament(x, y, z, shape="crossing"):
    a,b,c,d = 0,x+1,x+y+2,x+y+z+3
    out = [sum(1<<j for j in range(i+1,d+1)) for i in range(d+1)]
    pairs = {"crossing": [(a,c),(b,d)], "nested": [(a,d),(b,c)],
             "disjoint": [(a,b),(c,d)]}
    if shape in ("same head","chain"):
        d=c
        out=out[:c+1]
        out=[p & ((1<<(c+1))-1) for p in out]
        pairs.update({"same head":[(a,b),(a,c)], "chain":[(a,b),(b,c)]})
    for low, high in pairs[shape]:
        out[low] ^= 1<<high
        out[high] ^= 1<<low
    return out

def path_table(out, root=None):
    n = len(out)
    dp = [[0]*n for i in range(1<<n)]
    for i in range(n):
        dp[1<<i][i] = int(root is None or i==root)
    for s in range(1,1<<n):
        for i in range(n):
            p=s ^ (1<<i)
            if not s & (1<<i) or not p:
                continue
            dp[s][i] = sum(dp[p][j] for j in range(n) if out[j] & (1<<i))
    h=[sum(row) for row in dp]
    h[0]=int(root is None)
    return h

def hamilton(out):
    return path_table(out)[-1]

def rooted(out,u):
    A=path_table(out)
    P=path_table(out,u)
    val=[0,0]
    allv=(1<<len(out))-1
    for s in range(1,allv+1):
        val[s.bit_count()%2]+=P[s]*A[allv^s]
    return A[-1],val[0],val[1]

def ordinal(*parts):
    n=sum(map(len,parts))
    out=[]
    pos=0
    for part in parts:
        end=pos+len(part)
        suffix=((1<<n)-1)^((1<<end)-1)
        out.extend((row<<pos)|suffix for row in part)
        pos=end
    return out

def edit_pairs(out,pairs):
    out=out[:]
    for u,v in pairs:
        out[u]^=1<<v
        out[v]^=1<<u
    return out

def singleton_count(k):
    return ordinal(*([0] for i in range(k)))

def checks_joint_rooted():
    rng=Random(20261115)
    tri=[2,4,1]
    a=b=0
    for ix in range(90):
        P,C,S=(rng.choice(([],[0],tri,singleton_count(2),[2,12,9,1])) for k in range(3))
        # The last arbitrary tournament is used only as prefix or suffix.
        if C not in ([0],tri):
            C=rng.choice(([0],tri))
        W=ordinal(P,C,S)
        u=rng.randrange(len(C))
        h,E,O=rooted(W,len(P)+u)
        hC,eC,oC=rooted(C,u)
        hP=hamilton(P)
        pred=h+2*hP*eC
        if S:
            ht=path_table(S)
            q=sum(ht[A]*ht[(len(ht)-1)^A] for A in range(len(ht)))
            pred=h+hP*(eC+oC)*q
        wH=hamilton(edit_pairs(ordinal(W,[0]),[(len(P)+u,len(W))]))
        if not (wH==pred==h+2*E):
            raise ValueError(('root split',P,C,S,u,wH,pred,h,E))
        a+=1
    for ix in range(160):
        lengths=[rng.randrange(3) for j in range(3)]
        parts=[[0],singleton_count(lengths[0]),[0],singleton_count(lengths[1]),
               rng.choice(([0],tri)),singleton_count(lengths[2]),[0]]
        starts=[]
        starts_quot=[]
        cur=qu=0
        for i,part in enumerate(parts):
            starts.append(cur)
            starts_quot.append(qu)
            cur+=len(part)
            qu+=len(part) if i%2==0 else bool(part)
        anchor=[starts[i]+j for i in (0,2,4,6) for j in range(len(parts[i]))]
        q_anchor=[starts_quot[i]+j for i in (0,2,4,6) for j in range(len(parts[i]))]
        edits=rng.sample(list(combinations(range(len(anchor)),2)),rng.randrange(5))
        out=edit_pairs(ordinal(*parts),[(anchor[i],anchor[j]) for i,j in edits])
        quot=edit_pairs(ordinal(*(x if i%2==0 else [0]*bool(x) for i,x in enumerate(parts))),
                        [(q_anchor[i],q_anchor[j]) for i,j in edits])
        reach=[p|(1<<i) for i,p in enumerate(quot)]
        for k in range(len(quot)):
            for i in range(len(quot)):
                if reach[i]>>k&1:
                    reach[i]|=reach[k]
        prod=1
        for i,t in zip((1,3,5),lengths):
            if not t:
                continue
            v=starts_quot[i]
            is_active=any(w!=v and (reach[v]>>w&1) and (reach[w]>>v&1)
                          for w in range(len(quot)))
            if is_active:
                prod*=2**t-1
        if hamilton(out)<prod:
            raise ValueError(('joint bound',parts,edits,prod,out))
        b+=1
    return a,b

def formula(x,y,z):
    return (1 + 2**(x+y+1) + 2**(y+z+1) + 2**(x+z+1)*3**y
            - (x==0)*2**z - (z==0)*2**x)

def f(gap):
    return 2**gap if gap else 0

def formulas(x,y,z):
    return {
        "crossing":formula(x,y,z),
        "nested":1+f(y)+2**(x+z+1)*(2**(y+1)+3**y-1),
        "disjoint":(1+f(x))*(1+f(z)),
        "same head":1+f(x)+2**(x+y+1)-2**y-(y==0),
        "chain":1+f(x)+f(y)+2**(x+y)+(x+y==0),
    }

def decode_interior(N):
    """All ordered x<=z positive solutions, using exact modular factor filters."""
    if N<1 or N%2==0:
        return []
    candidates = []
    y = 0
    while 1 + 2**(y+3) + 8*3**y <= N:
        S,T=3**y,2**y
        M=T*T+S*((N-1)//2)
        x=1
        while (F:=S*2**x+T)**2 <= M:
            if M%F==0:
                G=(M//F)-T
                if G%S==0:
                    g=G//S
                    if g>0 and g&(g-1)==0:
                        candidates.append((x,y,g.bit_length()-1))
            x+=1
        y+=1
    return candidates

def decode_chamber(N):
    if N<=1 or N%2==0:
        return None
    W=N-1
    v=(W & -W).bit_length()-1
    Q=W>>v
    y=0
    p=1
    while p*3<=Q:
        p*=3
        y+=1
    R=Q-p
    bits=[i for i in range(R.bit_length()) if R>>i&1]
    if len(bits)!=2:
        return None
    x,z=y-bits[1],y-bits[0]
    if 1<=x<z<y and v==x+z+1:
        return x,y,z
    return None

def main():
    checked=0
    for x,y,z in product(range(7),repeat=3):
        if x+y+z>8:
            continue
        N=formula(x,y,z)
        for shape,prediction in formulas(x,y,z).items():
            found=hamilton(tournament(x,y,z,shape))
            if found!=prediction:
                raise ValueError((shape,x,y,z,prediction,found))
            checked+=1
        if 1<=x<=z and (x,y,z) not in decode_interior(N):
            raise ValueError(('factor inverse',x,y,z,N))
    seen=set()
    for x,z,y in combinations(range(1,37),3):
        N=formula(x,y,z)
        if N in seen or decode_chamber(N)!=(x,y,z) or (x,y,z) not in decode_interior(N):
            raise ValueError(('inverse collision',x,y,z,N))
        seen.add(N)
    print(f'{checked} independent tournament DP cases (including zero-gap boundaries)')
    print(f'{len(seen)} distinct exact chamber counts and both target decoders')
    print('first strict-chamber example: gaps (1,3,2), H =',formula(1,3,2))
    a,b=checks_joint_rooted()
    print(f'{a} rooted SCC split checks; {b} joint inflation checks')

if __name__=='__main__':
    main()
