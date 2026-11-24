from collections import Counter
from replay import graph, pentagons, components

def moment(k,c):
    # c is a census of exterior cycles of length 5*l, with sheet sum 55.
    p={0:1}
    for _ in range(k):
        q=Counter()
        for i,x in p.items():
            for j in (-1,0,1):q[i+j]+=x
        p=q
    tr=lambda m:m*sum(x for j,x in p.items() if j%m==0)
    center=sum(m*tr(5*l) for l,m in c.items())-tr(5)
    return 52**k+1454*7**k+1245*(-8)**k+(-1)**k*center

if __name__=='__main__':
    for c in ({1:55},{1:53,2:1},{55:1}):
        a=[moment(k,c) for k in range(11)]
        assert a[:5]==[2970,0,154440,0,15907320]
        assert a[5]==10*(36383149-c.get(1,0))
        assert a[6]==2970*52*(5*52**2-6*52+2)+12*(1518993850+5*c.get(1,0))
        assert a[10]==144556853479042730+14530*c.get(1,0)+20*c.get(2,0)
    G=graph(7);seen=0
    for C in pentagons(G):
        X=set.union(*(G[v] for v in C))-set(C)
        Y=set(range(len(G)))-X-set(C)
        assert sorted(map(len,components(G,X)))==[5]*5
        assert sorted(map(len,components(G,Y)))==[5]*4
        seen+=1
    assert seen==1260
    print('exact residual moment constants and 1260 pentagons checked')
