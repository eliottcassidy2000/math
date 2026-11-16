from functools import lru_cache

def profile(A,root):
    N=len(A);F=(1<<N)-1
    @lru_cache(None)
    def paths(S,j):
        if S==1<<j:return 1
        return sum(paths(S^(1<<j),i) for i in range(N) if i!=j and S>>i&1 and A[i]>>j&1)
    @lru_cache(None)
    def total(S):
        return sum(paths(S,j) for j in range(N) if S>>j&1) if S else 1
    @lru_cache(None)
    def starts(S,j):
        if S==1<<j:return int(j==root)
        return sum(starts(S^(1<<j),i) for i in range(N) if i!=j and S>>i&1 and A[i]>>j&1)
    r=[0,0]
    for S in range(1,F+1):
        if S>>root&1:
            r[S.bit_count()%2]+=total(F^S)*sum(starts(S,j) for j in range(N) if S>>j&1)
    return total(F),total(F^(1<<root)),*r

def append(A,root):
    n=len(A)
    return tuple(z+((i!=root)<<n) for i,z in enumerate(A))+(1<<root,)

if __name__=='__main__':
    A=(6,12,40,33,15,83,159,63)
    state=(107,37,69,89)
    assert profile(A,0)==state
    for n in range(1,7):
        h,g,r,s=state
        state=(h+2*r,g,r+s-g,r+s)
        A=append(A,0)
        assert profile(A,0)==state
        assert state[0]==3+121*2**n
    assert state[0]==7747==4*(11*2**2)**2+3
    print("Verified rooted seed and literal first six appends; first reopen at 7747.")

    # Relabel the common spacer to the same physical index on both sides.
    from itertools import combinations
    base=(6,12,40,33,15,83,159,63)
    pad=tuple(row|(7<<8) for row in base)+(1536,1024,0)
    old_roles=(6,7,4,3,0,8,9,10,5)
    old_marks={1,3,4,5,7}
    old_before={0,1,2,3,8}
    old_rows=tuple(
        sum(1<<j for j in range(i+1,8))+((i in old_marks)<<8)
        for i in range(8)
    )+(sum(1<<j for j in range(8) if j not in old_marks),)
    flips={(i,j) for i,j in combinations(range(9),2)
           if ((old_rows[i]>>j)^(pad[old_roles[i]]>>old_roles[j]))&1}
    assert flips=={(4,8),(5,8),(7,8)}
    assert old_before=={i for i,r in enumerate(old_roles) if 1<=r<=7}
    assert set(range(8))-set(old_roles)=={1,2}
    print("Verified fixed boundary assignment, 3 old flips, 2 new vertices, unchanged spacer.")

    B=(2,36,65,23,135,93,155,47)
    assert profile(B,1)==(247,31,121,152)
    new_roles=(1,4,2,0,3,6,7,5)
    assert new_roles[0]==1
    assert {(i,j) for i,j in combinations(range(8),2)
            if ((base[i]>>j)^(B[new_roles[i]]>>new_roles[j]))&1}=={
            (0,1),(0,7),(1,6),(3,5)}
    for _ in range(5):
        B=append(B,1)
    assert profile(B,1)[0]==7749
    print("Verified 7749 and the five-flip matched second reopen.")
