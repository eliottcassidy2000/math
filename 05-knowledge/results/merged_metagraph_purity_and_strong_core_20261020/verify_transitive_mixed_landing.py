from math import comb

def candidate(n):
    A=[[0]*n for _ in range(n)]
    for a in range(n):
        for b in range(a):
            if a-b==1:A[a][b]=1
            else:A[b][a]=1
    P=[0]+list(range(n-1,0,-1))
    assert all(A[a][b] for a,b in zip(P,P[1:]))
    # Path order is descending in the fixed-path model.
    old_at_new=list(reversed(P))
    B=[[A[old_at_new[i]][old_at_new[j]] for j in range(n)] for i in range(n)]
    tiles=[(x,y) for x in range(n) for y in range(x-1)]
    bit={(x,y):B[x][y] for x,y in tiles}
    support={xy for xy,v in bit.items() if v}
    predicted={(n-1,y) for y in range(1,n-2)}
    assert support==predicted
    fixed=all(bit[x,y]==bit[n-1-y,n-1-x] for x,y in tiles)
    assert fixed==(n==3)
    scores=list(map(sum,A)); c3=comb(n,3)-sum(comb(s,2) for s in scores)
    assert scores[0]==scores[1]==n-2 and c3==n-2
    assert scores[P[0]]+scores[P[-1]]==2*n-4
    return dict(n=n,path=P,c3=c3,nonblue_normalized_support=sorted(support),endpoint_score_sum=2*n-4,required_if_blue=n-1)

for n in range(3,101):
    row=candidate(n)
    if n<=10:print(row)
print('ALL CHECKS PASSED for n=3..100; direct matrix/path/normalization/reflection/score and C3 identities, no graph-isomorphism heuristic.')
