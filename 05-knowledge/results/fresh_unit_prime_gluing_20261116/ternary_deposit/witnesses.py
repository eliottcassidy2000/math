from pathlib import Path
f=Path(__file__).with_name("results.txt")
for line in f.read_text().splitlines():
    if " example " not in line:
        continue
    pre,data=line.split(" example ")
    q,m,h,_=map(int,pre.split())
    a=list(map(int,data.split()))
    n=len(a)
    for i in range(n):
        assert 0<=a[i]<1<<n and not(a[i]>>i&1)
        for j in range(i):
            assert (a[i]>>j&1)+(a[j]>>i&1)==1
    F=(1<<n)-1
    dp=[[0]*n for _ in range(F+1)]
    counts=[1]+[0]*F
    for s in range(1,F+1):
        for i in range(n):
            if s>>i&1:
                p=s^(1<<i)
                dp[s][i]=sum(dp[p][j] for j in range(n) if p>>j&1 and a[i]>>j&1) if p else 1
        counts[s]=sum(dp[s])
    e=sum(dp[s][-1]*counts[F^s] for s in range(1,F+1) if not s.bit_count()%2)
    o=sum(dp[s][-1]*counts[F^s] for s in range(1,F+1) if s.bit_count()%2)
    assert tuple(m*x for x in (counts[F],e,o,counts[F^(1<<(n-1))]))==(h,*({79:(378,423,45),163:(756,801,45)}[q]))
print("All recorded terminal witnesses verified via full rooted subset DP.")
