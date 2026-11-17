from pathlib import Path
from math import isqrt
import json
from solve import H,profile,literal,count
from find import tournaments,banks
T=Path(__file__).with_name('targets7.jsonl')
check={97309:15,681163:17,8191:11,57337:14,65537:14,458759:17,30011:13,210077:16,30013:13,210091:15}
for line in T.read_text().splitlines():
    z=json.loads(line)
    N=z['target'];D=z['D'];q=z['q'];u=z['u'];w=z['w']
    assert z['order']==check[N]
    assert tuple(z['profile'])==profile(D,q)
    assert count(z['profile'],u,w)==N
    G=literal(D,q,u,w);n=len(G)
    assert G==tuple(z['rows']) and len(G)==z['order']
    assert all(not (r>>i&1) and r>>n==0 for i,r in enumerate(G))
    assert all((G[i]>>j&1)+(G[j]>>i&1)==1 for i in range(n) for j in range(i))
    assert H(G)==N
    print(N,n,z['residue_prune'])
for p in (8191,65537,30011,30013):assert all(p%d for d in range(2,isqrt(p)+1))
assert 97309==31*43*73
assert [len(tournaments(n)) for n in range(1,8)]==[1,1,2,4,12,56,456]
print('PASS all ten literal DP witnesses, prime labels and center enumerator; unique marked profiles',len(banks(7)))
Y=(990,1020,3832,8112,4128,4993,3000,3600,3732,4112,627,5683,1479)
assert all(not (r>>i&1) and r>>len(Y)==0 for i,r in enumerate(Y))
assert all((Y[i]>>j&1)+(Y[j]>>i&1)==1 for i in range(len(Y)) for j in range(i))
assert H(Y)==681163
print('PASS independent order-13 apex witness for 681163')
