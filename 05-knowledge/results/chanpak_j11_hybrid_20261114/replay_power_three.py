"""Literal verifier for a strong tournament with H = 2^k+3, every k>=6."""
import sys
from functools import lru_cache
from collections import Counter

CORES={13:(130,420,265,403,391,477,415,260,1),19:(130,388,465,263,427,399,443,264,1),21:(130,388,409,275,419,463,415,264,1),1:(2,4,1),3:(2,12,17,37,75,23,47),5:(2,12,17,37,43,71,31),7:(2,12,17,69,43,15,55),9:(2,4,25,19,35,79,31)}

def ends(R):
    n=len(R)
    B=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        B[1<<i][i]=1
    for U in range(1,1<<n):
        for t,c in enumerate(B[U]):
            if c:
                V=R[t]&~U
                while V:
                    v=V&-V
                    V^=v
                    B[U|v][v.bit_length()-1]+=c
    return B

def validate(R):
    n=len(R)
    assert all(0<=r<1<<n and not(r>>i&1) for i,r in enumerate(R))
    assert all(((R[i]>>j)^(R[j]>>i))&1 for i in range(n) for j in range(i+1,n))

def graph(k,delta=3):
    core=CORES[delta]
    v=8 if delta in (13,19,21) else 0
    shift={1:0,3:5,5:5,7:5,9:4,13:7,19:7,21:7}[delta]
    n=k-shift
    if n<1:
        raise ValueError("insufficient module length")
    projection=list(range(v))+[v]*n+list(range(v+1,len(core)))
    R=tuple(sum(1<<j for j,b in enumerate(projection)
                if (a==b==v and i<j) or ((core[a]>>b)&1))
            for i,a in enumerate(projection))
    return R

def double_graph(n,m):
    core=CORES[7]
    projection=[0]*n+[1,2,3,4]+[5]*m+[6]
    return tuple(sum(1<<j for j,b in enumerate(projection)
                if (a==b and i<j) or ((core[a]>>b)&1))
            for i,a in enumerate(projection))

def skeleton_counts(R, first, second):
    F=(1<<len(R))-1
    exterior=F^((1<<first)|(1<<second))

    @lru_cache(None)
    def add(last, used, a, b):
        C=Counter()
        if used==exterior and a and b:
            C[a,b]+=1
        options=(F if last<0 else R[last])&~used
        for j in range(len(R)):
            if not(options>>j&1):continue
            if j==first: C.update(add(j,used,a+1,b))
            elif j==second: C.update(add(j,used,a,b+1))
            else: C.update(add(j,used|(1<<j),a,b))
        return C

    return add(-1,0,0,0)

def test_21():
    assert skeleton_counts(CORES[7],0,5)=={
        (1,1):71,(1,2):51,(1,3):10,
        (2,1):32,(2,2):22,(2,3):4,
    }
    for n in (1,2,3,6,9):
        T=double_graph(n,2)
        validate(T)
        assert sum(ends(T)[-1])==76*(1<<n)+21
    print("PASS two-module coefficients and 5 literal counts, including 38933.")

def test():
    for delta,core in CORES.items():
        expected1={1:3,3:67,5:69,7:71,9:41,13:269,19:275,21:277}[delta]
        expected2={1:1,3:32,5:32,7:32,9:16,13:128,19:128,21:128}[delta]
        F=(1<<len(core))-1
        validate(core)
        A=ends(core)
        assert sum(A[F])==expected1
        H=[sum(a) for a in A]
        H[0]=1
        rev=tuple((F^r)^(1<<i) for i,r in enumerate(core))
        B=ends(rev)
        v,w=(8,0) if delta in (13,19,21) else (0,1)
        D=F^(1<<v)^(1<<w)
        by_size=[0]*len(core)
        for Q in range(1,D+1):
            if Q & ~D: continue
            by_size[Q.bit_count()]+=H[D^Q]*B[Q|(1<<w)][w]
        assert sum(by_size)==expected2
        for k in range(8 if v==8 else 6,13):
            R=graph(k,delta)
            validate(R)
            assert sum(ends(R)[-1])==(1<<k)+delta
        print("Offset",delta,"two-visit counts by suffix size",by_size)
    print("PASS exact skeleton counts and 50 literal blow-ups.")

if __name__=='__main__':
    test()
    test_21()
    if len(sys.argv)>1:
        N=int(sys.argv[1])
        found=False
        for delta in CORES:
            q=N-delta
            if q<(256 if delta in (13,19,21) else 64) or q & (q-1):
                continue
            k=q.bit_length()-1
            print("WITNESS", N, graph(k,delta))
            found=True
        D=N-21
        if D>0 and D%76==0:
            q=D//76
            if q>=2 and not(q&(q-1)):
                print("WITNESS",N,double_graph(q.bit_length()-1,2))
                found=True
        if not found:
            print("Not in the certified ray family.")
