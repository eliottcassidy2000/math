"""Bounded discovery of uniform-incidence circulants with different H."""
from collections import defaultdict
from itertools import combinations
import json


def circulant(n, mask):
    S = {i if mask >> (i-1) & 1 else n-i for i in range(1,(n+1)//2)}
    return [[int((v-u)%n in S) for v in range(n)] for u in range(n)], sorted(S)


def tag(A):
    n=len(A)
    c3=sum(A[0][i]==A[i][j]==A[j][0] for i,j in combinations(range(1,n),2))
    s4=0
    for other in combinations(range(1,n),3):
        B=(0,)+other
        s4+=sorted(sum(A[u][v] for v in B) for u in B)==[1,1,2,2]
    return sum(A[0]),c3,s4


def paths(A):
    n=len(A); rows=[sum(A[u][v]<<v for v in range(n)) for u in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for u in range(n):dp[1<<u][u]=1
    for mask in range(1,1<<n):
        for u,value in enumerate(dp[mask]):
            if not value:continue
            left=rows[u]&~mask
            while left:
                bit=left&-left;left-=bit;v=bit.bit_length()-1
                dp[mask|bit][v]+=value
    return sum(dp[-1])


if __name__=='__main__':
    for n in (9,11,13):
        banks=defaultdict(list)
        for mask in range(1<<((n-1)//2)):
            A,S=circulant(n,mask); key=tag(A); H=paths(A)
            for old in banks[key]:
                if H!=old['H']:
                    print(json.dumps({'n':n,'tag':key,'first':old,'second':{'mask':mask,'S':S,'H':H}}),flush=True)
                    raise SystemExit
            banks[key].append({'mask':mask,'S':S,'H':H})
        print(json.dumps({'n':n,'cases':sum(map(len,banks.values())),'tags':len(banks),'collision':False}),flush=True)
