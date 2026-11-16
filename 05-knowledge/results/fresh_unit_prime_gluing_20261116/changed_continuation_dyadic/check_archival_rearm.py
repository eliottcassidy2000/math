from itertools import product

CHILDREN = {
    2:(118,372,504,115,224,192,128,11,249),
    4:(118,124,504,369,224,192,128,11,243),
}
def restrict(a,S):
    return tuple(sum((a[i]>>j&1)<<k for k,j in enumerate(S)) for i in S)

def reachable(a,seed,S):
    x=1<<seed
    while True:
        y=x
        q=x
        while q:
            c=q&-q;q^=c
            y|=a[c.bit_length()-1]&S
        if x==y:return x
        x=y

def terminal(a,skip):
    S=((1<<len(a))-1)^(1<<skip); t=S
    for i in range(len(a)):
        if S>>i&1:t&=reachable(a,i,S)
    return t

def append(a,root):
    return tuple(x+((i!=root)<<len(a)) for i,x in enumerate(a))+(0 if root<0 else 1<<root,)

def profiles(a):
    n=len(a); f=(1<<n)-1
    starts=[[0]*n for _ in range(f+1)]
    H=[1]+[0]*f
    for S in range(1,f+1):
        for i in range(n):
            if S>>i&1:
                M=S^(1<<i)
                starts[S][i]=sum(starts[M][j] for j in range(n) if a[i]>>j&1) if M else 1
        H[S]=sum(starts[S])
    out=[]
    for i in range(n):
        E=sum(starts[S][i]*H[f^S] for S in range(1,f+1) if not S.bit_count()%2)
        O=sum(starts[S][i]*H[f^S] for S in range(1,f+1) if S.bit_count()%2)
        g=H[f^(1<<i)]
        out.append((H[-1],E+O,O-E,g,H[-1]-2*O+3*g))
    return out

for d,a in CHILDREN.items():
    for empty in (False,True):
        p=tuple(i for i in range(9) if not(empty and i==6))
        C=restrict(a,p); root=p.index(7)
        for i in range(len(C)):
            t=terminal(C,i)
            if i==root:assert not (C[i]&t)
            else:assert t.bit_count()==len(C)-1
        assert profiles(C)[root][-1]==d
        for n in range(4):
            for word in product((0,1),repeat=n):
                T=C
                for e in word:
                    T=append(T,root if e==0 else -1)
                for H,B,z,g,Delta in profiles(T):
                    assert z!=g or Delta!=0
        T=C; last=root; tower=[profiles(T)[last]]
        for k in range(1,5):
            prior=last;last=len(T); T=append(T,prior)
            P=profiles(T);tower.append(P[last])
            normalized=[i for i,(_,_,z,g,_) in enumerate(P) if z==g]
            positions=[root]+list(range(len(C),len(T)))
            assert set(normalized)==set(positions[max(0,k-2):k])
            for i in normalized:
                assert P[i][-1]!=0
            if k>=2:assert P[positions[k-1]][-1]==2*tower[k-2][0]
            if k>=3:assert P[positions[k-2]][-1]==2*tower[k-3][0]-tower[k-3][2]-tower[k-3][1]
print('Checked both nonzero digits, empty/nonempty slot, binary reroots, and four moving tower heights.')
