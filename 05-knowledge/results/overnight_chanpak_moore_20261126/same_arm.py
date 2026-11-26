from fixed_bank_inverse import count
from pathlib import Path
from sys import argv
bank={tuple(map(int,l.split())) for l in Path(argv[1]).read_text().splitlines() if l.strip()}
hits=[]
for h in (1,3):
    for R in range(h,h+3):
        for e in range(h+1):
            for v in range(32):
                for i,j in [(0,1),(0,2),(1,2)]:
                    H=count(v,h,e,R)
                    a,b=count(v^(1<<i),h,e,R)-H,count(v^(1<<j),h,e,R)-H
                    if {a,b}=={-2,2}:
                        row=h,e,R
                        hits.append((row,v,i,j,H,(a,b),row in bank))
print(*hits,sep='\n')

def T(v,R,n):
    # h=1 profiles in the valid rows have a singleton or a source in a two-chain on the right.
    rn=1 if R==1 else 2
    rows=[sum(1<<k for k in range(j+1,n))+((v>>j&1)<<n) for j in range(n)]
    rows += [((1<<n)-1)^v | 1<<(n+1)]
    rows += [((1<<n)-1)+(1<<(n+2) if rn==2 else 0)]
    if rn==2:rows += [(1<<(n+1))-1]
    return rows

def dp(rows):
    n=len(rows)
    assert all((r>>i&1)==0 and (r>>n)==0 for i,r in enumerate(rows))
    assert all((r>>j&1)+(rows[j]>>i&1)==1 for i,r in enumerate(rows) for j in range(i))
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for S in range(1,1<<n):
        for i,c in enumerate(d[S]):
            rem=rows[i]&~S
            while rem and c:
                b=rem&-rem;rem^=b
                d[S|b][b.bit_length()-1]+=c
    return sum(d[-1])
for row,v,i,j,H,ab,valid in hits:
    if not valid:continue
    n=max(v.bit_length(),i+1,j+1)
    t,x,y=(T(v ^ z,row[2],n) for z in (0,1<<i,1<<j))
    assert sum((a^b).bit_count() for a,b in zip(t,x))==2
    assert sum((a^b).bit_count() for a,b in zip(t,y))==2
    assert tuple(map(dp,(t,x,y)))==(H,H+ab[0],H+ab[1])
print('literal path DP and single-spoke edits checked for all four genuine controls')
