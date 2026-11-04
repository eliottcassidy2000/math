from itertools import permutations

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def valid(rows):return all(not(row>>i&1)and all(((row>>j)&1)+((rows[j]>>i)&1)==1 for j in range(len(rows))if i!=j)for i,row in enumerate(rows))

def strong(rows):
    n=len(rows);reach=1
    while True:
        nxt=reach
        for i in range(n):
            if reach>>i&1:nxt|=rows[i]
        if nxt==reach:break
        reach=nxt
    if reach!=(1<<n)-1:return False
    rev=[sum(1<<j for j in range(n)if rows[j]>>i&1)for i in range(n)];reach=1
    while True:
        nxt=reach
        for i in range(n):
            if reach>>i&1:nxt|=rev[i]
        if nxt==reach:break
        reach=nxt
    return reach==(1<<n)-1

def endpoints(rows):
    n=len(rows);d=[[0]*n for _ in range(1<<n)]
    for mask in range(1,1<<n):
        for j in range(n):
            if mask>>j&1:
                prev=mask^(1<<j);d[mask][j]=1 if not prev else sum(d[prev][i]for i in range(n)if prev>>i&1 and rows[i]>>j&1)
    return d[-1]

def delete(rows,a):
    ids=[i for i in range(len(rows))if i!=a];return tuple(sum(1<<j for j,v in enumerate(ids)if rows[u]>>v&1)for u in ids),ids

def C2(rows,a):
    n=len(rows);out=list(rows)+[1<<a,1<<n]
    for i in range(n):
        if i!=a:out[i]|=1<<n
        out[i]|=1<<(n+1)
    return tuple(out)

def phi(rows,a):
    root=rows[a].bit_length()-1;need(rows[a]==1<<root,'sole exit');out=list(rows);ids=[i for i in range(len(rows))if i!=a]
    for i in ids:
        for j in ids:
            if i!=j:out[i]^=1<<j
    return tuple(out),root
