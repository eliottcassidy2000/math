from itertools import product
from marked_join import rows,dp,packet,join,joins,close
from math import gcd

A,B=[(rows(6,c),o) for c,o in [(257,51),(8737,29)]]
packets=[packet(*a) for a in (A,B)]
assert packets==[(17,224,135,135,11,11,76),(17,224,135,135,11,11,72)]
assert [dp(close(a,o,[0]))[0][-1] for a,o in (A,B)]==[211,215]
for (a,o),(b,O) in product((A,B),repeat=2):
    c=join(a,b)
    direct=dp(close(c,o|(O<<len(a)),[0]))[0][-1]
    h,q,M,E,s,e,W=joins(packet(a,o),packet(b,O))
    assert direct==h-W+2*M
for k in range(3,25):
    outputs=[]
    for j in range(k+1):
        word=[packets[0]]*(k-j)+[packets[1]]*j
        from functools import reduce
        h,q,M,E,s,e,W=reduce(joins,word)
        for t in range(1,17):
            assert gcd(h-W+M*2**t,17)==1
        outputs.append(h-W+2*M)
    assert all(b-a==4*17**(k-1) for a,b in zip(outputs,outputs[1:]))
print('2 literal primitive hosts, 4 literal two-block closes, 22 repeatable progression lengths')
