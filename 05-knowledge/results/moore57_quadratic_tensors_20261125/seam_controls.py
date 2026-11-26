"""Independent checks for nearest covers and the cyclic Thue–Morse filter."""

from pathlib import Path
import sys

from slot_interface_check import BRANCHES, read_board


def nearest(board_path, nearest_path):
    ref, raw = read_board(board_path)
    words = [[],[],[]]
    for line in Path(nearest_path).read_text().splitlines():
        g,d,*ys=map(int,line.split())
        w=(7,*ys)
        assert d==[7,5,6][g] and len(w)==53
        assert sum(y!=z for y,z in zip(w,ref) if z>=0)==d
        T=[D[y] for D,y in zip(raw,w)]
        assert all(len({t[k] for t in T})==53 for k in range(3) if k!=g)
        words[g].append(w)
    assert list(map(len,words))==[23,1,1] and all(len(set(V))==len(V) for V in words)
    u,v=words[1][0],words[2][0]
    mp={y:i for i,y in enumerate(u)}
    D={i for i,(a,b) in enumerate(zip(u,v)) if a!=b}
    assert len(D)==13
    cycle=[]
    i=min(D)
    while i not in cycle:
        cycle.append(i)
        i=mp[v[i]]
    assert i==cycle[0] and set(cycle)==D
    assert [BRANCHES[h] for h in cycle]==[5,39,25,53,40,33,20,36,15,41,13,27,52] or \
           [BRANCHES[h] for h in cycle]==[5,52,27,13,41,15,36,20,33,40,53,25,39]
    for g,w in [(1,u),(2,v)]:
        assert 53-len({D[y][g] for D,y in zip(raw,w)})==7
    print("Nearest projected repairs = 23, 1, 1; unique pivot seam cycle:",[BRANCHES[h] for h in cycle])

    pool=[set() for _ in range(53)]
    for V in words:
        for w in V:
            for h,y in enumerate(w):pool[h].add(y)
    assert sum(map(len,pool[1:]))==215
    B={h:[(y,*raw[h][y]) for y in pool[h]] for h in range(1,53)}
    # First rule out >=51 changes entirely in the 215-bank, then prove
    # let the optimized C++ replay check all outside contacts.
    def inside_only(A,skips):
        if not A:return True
        h=min(A,key=lambda j:len(A[j]))
        B={j:V for j,V in A.items() if j!=h}
        if skips and inside_only(B,skips-1):return True
        return any(inside_only({j:[f for f in V if all(x!=y for x,y in zip(f[1:],e[1:]))]
                                for j,V in B.items()},skips) for e in A[h])
    assert not inside_only(B,1)
    # The exhaustive import bound has a separate optimized C++ replay.
    print("The 215-contact bank cannot even cover 51 physical branches.")


def rank_mod(M,p):
    A=[[x%p for x in row] for row in M]
    n=len(A);r=0
    for j in range(n):
        pivot=next((i for i in range(r,n) if A[i][j]),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r]
        d=pow(A[r][j],-1,p)
        for i in range(r+1,n):
            s=d*A[i][j]
            A[i]=[(x-s*y)%p for x,y in zip(A[i],A[r])]
        r+=1
    return r


def cyclic():
    c=[1]+[0]*13
    ranks=[]
    for m in range(1,10):
        s=2**(m-1)%14
        c=[c[j]-c[(j-s)%14] for j in range(14)]
        M=[[c[(i-j)%14] for j in range(14)] for i in range(14)]
        assert rank_mod(M,101)==(13 if m==1 else 12)
        r=rank_mod(M,7)
        assert r==max(7-m,0)+max(8-m,0)
        ranks.append(r)
    print("Modulo-seven ranks for depths 1..9:",ranks)


if __name__=='__main__':
    nearest(*sys.argv[1:])
    cyclic()
