"""Verify exact carry families and width-seven fixture, independently."""
from itertools import combinations


def table(n, V):
    E=sum(1<<v for v in V)
    L=sum(1<<(n-1-v) for v in V)
    J=sum((1<<(v-u-1))+(u+1==v) for u,v in combinations(sorted(V),2))
    return E,L,J


def direct(n, V):
    # x first: one increasing tape path, then t.
    # t before x: divide the tape into slots ending at t / after x.
    # A mark before x: first slot must end in V, suffix must start unmarked.
    full=(1<<n)-1
    mv=sum(1<<v for v in V)
    assert not (mv&1 or mv>>(n-1))
    count=1
    for A in range(1<<n):
        B=full^A
        legal=not (B & -B & mv)
        if legal:
            count+=1
            if A and mv>> (A.bit_length()-1)&1:
                count+=1
    return count


def row(n,V, expected):
    E,L,J=table(n,V)
    assert 1+(1<<n)+E-L-J==expected
    assert direct(n,V)==expected
    return J


def main():
    for m in range(1,5):
        N=4*m+2
        base=(1<<N)-(1<<(2*m))+1
        assert row(N,{2*m},base)==0
        assert row(N,set(range(m+1,3*m+1)),base+2)==(1<<(2*m))-2
    for r in range(2,7,2):
        N=3*r+3; c=(N-1)//2
        deficit=(1<<(r+1))*((1<<r)-1)
        a={r+1}
        b=set(range(c-r,c))|set(range(c+1,c+r+1))
        base=(1<<N)-deficit+1
        assert row(N,a,base)==0
        assert row(N,b,base+2)==deficit-2
    words=['11101111010011','11100010000011','11000100010101','01100111101010',
           '01011001010010','11001000110101','10010011111110']
    Js=[]
    for i,w in enumerate(words):
        V={3+j for j,c in enumerate(reversed(w)) if c=='1'}
        Js.append(row(20,V,1054527+2*i))
    assert Js==[12194,10912,8270,2268,1658,8648,4702]
    assert all(v-u!=-2 for u,v in zip(Js,Js[1:]))
    for i in range(2,6):
        n=3*i+3
        U={i-1,i+1,2*i+2}; V={i-1,2*i+1}
        assert len(U^V)==3
        base=(1<<n)-6*(1<<(2*i))-7*(1<<(i-1))-1
        row(n,U,base)
        row(n,V,base+2)
        e1,l1,j1=table(n,U); e2,l2,j2=table(n,V)
        assert e2-l2-e1+l1==-3*(1<<i)
        assert j2-j1==-3*(1<<i)-2
    for n in range(8,15):
        sites=range(3,n-3)
        G=[]
        for b in range(1<<len(sites)):
            M={v for j,v in enumerate(sites) if b>>j&1}
            E,L,J=table(n,M); G.append((E-L,J))
        for b,(s0,j0) in enumerate(G):
            for c,(s1,j1) in enumerate(G):
                if (b^c).bit_count()<=2 and s0!=s1:
                    assert abs((s1-j1)-(s0-j0))>2
    print('29 literal counts and exhaustive small edit tests passed')


if __name__=='__main__': main()
