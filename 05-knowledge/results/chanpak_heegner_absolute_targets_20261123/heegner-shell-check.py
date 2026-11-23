"""Literal replay for the all-exponents signed Mersenne rail."""
def strong(rows):
    for start in range(len(rows)):
        seen=1<<start
        for _ in rows:
            for i,out in enumerate(rows):
                if seen>>i&1: seen |= out
        if seen != (1<<len(rows))-1: return False
    return True

def weighted(word):
    h,u,previous=1,0,0
    for w in word:
        bit=int(w=="O")
        if not bit: h+=u+previous
        u=2*u+bit; previous=bit
    return h


CORE=(50,84,25,19,32,78,29)
I=35

def endpoint(rows):
    n=len(rows)
    assert all(not (row>>i&1) and row>>n==0 for i,row in enumerate(rows))
    assert all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i in range(n) for j in range(i))
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask,row in enumerate(dp):
        for i,v in enumerate(row):
            if v:
                left=rows[i]&~mask
                while left:
                    b=left&-left;left-=b
                    dp[mask|b][b.bit_length()-1]+=v
    return dp

def reads(rows,tags):
    total=(1<<len(rows))-1
    D=endpoint(rows)
    DR=endpoint([r^total^(1<<i) for i,r in enumerate(rows)])
    G=lost=0
    for a in range(1,total+1):
        b=total^a
        end=sum(v for i,v in enumerate(D[a]) if tags>>i&1)
        G+=end*(sum(DR[b]) if b else 1)
        lost+=end*sum(v for i,v in enumerate(DR[b]) if tags>>i&1)
    return sum(D[-1]),G,lost

def shell(s):
    n=len(CORE)
    f=((1<<s)-1)<<(n+1)
    output=[r|((I>>i&1)<<n) for i,r in enumerate(CORE)]
    output+=[(((1<<n)-1)^I)|f]
    output += [((1<<n)-1)|(f&~((1<<(i+1))-1)) for i in range(n+1,n+1+s)]
    return output

if __name__=='__main__':
    assert reads(CORE,I)==(121,512,128)
    assert strong(CORE) and strong(shell(1))
    for s in range(1,7):
        R=shell(s)
        assert strong(R)
        assert sum(endpoint(R)[-1])==(1<<(s+9))-7
        print('signed shell',s+9,len(R),(1<<(s+9))-7)
    R=(102,4,8,19,39,78,30)
    assert strong(R) and sum(endpoint(R)[-1])==57
    for q,w in [(8,"O"+"I"*4+"O"+"I"*3),
                (9,"O"+"I"*2+"O"+"I"*3+"O"+"I"*2)]:
        assert weighted(w)==(1<<q)-7
    print("small even q=6,8 witnessed")
