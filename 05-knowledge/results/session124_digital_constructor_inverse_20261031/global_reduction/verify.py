"""Check the explicit H=63 witness without using digital/profile identities."""
from itertools import combinations,permutations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def require(value,message):
    if not value:raise RuntimeError(message)

def reconstruct(record):
    n,code,cut=record['seed'];word=record['blocks'];N=n+len(word)+1
    require(all(b in ('0','1') for b in word),'Expected singleton witness')
    A=[[0]*N for _ in range(N)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j]=(code>>k)&1;A[j][i]=1-A[i][j]
    for j in range(n,N-1):
        for i in range(j):A[i][j]=1
    bits=[(cut>>i)&1 for i in range(n)]+list(map(int,word))
    for i,b in enumerate(bits):A[N-1][i]=b;A[i][N-1]=1-b
    return A

def backward_endpoint_dp(A):
    n=len(A);dp=[[0]*n for _ in range(1<<n)]
    for S in range(1,1<<n):
        for v in range(n):
            if S>>v&1:
                rest=S^(1<<v)
                dp[S][v]=1 if not rest else sum(dp[rest][u] for u in range(n) if rest>>u&1 and A[u][v])
    return sum(dp[-1])

record=json.loads((ROOT/'witness63.json').read_text());A=reconstruct(record);n=len(A)
require(n==record['order']==9,'Wrong order')
require([''.join(map(str,row)) for row in A]==record['adjacency_rows'],'Recipe/adjacency discrepancy')
require(all(A[i][i]==0 for i in range(n)) and all(A[i][j]+A[j][i]==1 for i in range(n) for j in range(i+1,n)),'Not a tournament')
count=backward_endpoint_dp(A)
literal=sum(all(A[u][v] for u,v in zip(word,word[1:])) for word in permutations(range(n)))
require(count==literal==record['target']==63,'Actual H differs from 63')
triangle=[[0,1,0],[0,0,1],[1,0,0]]
require(backward_endpoint_dp(triangle)==3,'Triangle H differs from 3')
print(json.dumps({'status':'PASS','order':n,'target':63,'backward_endpoint_dp':count,'literal_permutation_count':literal,'permutations_examined':362880,'triangle_path_count':3,'profile_or_digital_formula_used':False,'scope':'Actual H=63 and H=3 witnesses only; ordinal-product reduction is proved separately; no coverage census or new exclusion proof.'},sort_keys=True))
