#!/usr/bin/env python3
"""Full C/D owner universe at k=4, without claiming a full mass-layer census."""
from collections import Counter,defaultdict
from itertools import combinations
import json
import extension_check as X
Z=X.Z


def main():
    k=4;rows=[]
    for strong in (False,True):
        recs=[];n=2*k-1;m=n-1 if strong else n+1
        for edges in Z.p3_graphs(k):
            if strong and not Z.indecomp(edges,n):continue
            for ds in Z.compositions(k-1,m,k-1):
                recs.append(('C',edges,tuple(k-d for d in ds),n))
        n=2*k;m=n-1 if strong else n+1
        for edges in Z.matchings(tuple(range(n))):
            if strong and not Z.indecomp(edges,n):continue
            for holes in combinations(range(m),2):
                recs.append(('D',edges,tuple(0 if i in holes else k for i in range(m)),n))
        fibres=defaultdict(list);accepted=Counter();rejected=0
        for kind,edges,gaps,n in recs:
            t,w=Z.model(n,edges,gaps,strong);clones=Z.clones(t)
            if kind=='D':
                actual=(0,)+gaps+(0,) if strong else gaps
                bad=any(b==a+1 and actual[b]==0 for a,b in edges)
                Z.need(bool(clones)==bad,'two-gap clone criterion')
                if bad:rejected+=1;continue
            Z.need(not clones,'retained C/D core')
            if strong:Z.need(Z.strongly_connected(t),'retained strong C/D core')
            key,_=Z.canonical(t,w);fibres[key].append(kind);accepted[kind]+=1
        hist=Counter(tuple(sorted(v)) for v in fibres.values())
        Z.need(set(hist)<=set([('C',),('D',),('C','C','D')]),'only the claimed owner fibre')
        h=60 if strong else 90
        Z.need(hist[('C','C','D')]==h,'all marked span-two owner fibres')
        Z.need(len(fibres)==(14922 if strong else 40605),'C/D consequence class count')
        rows.append(dict(k=k,strong=strong,scope='C/D owner families only; not the complete mass layer',
                         classes=len(fibres),accepted=dict(accepted),clone_rejections=rejected,
                         fibres={','.join(t):v for t,v in hist.items()}))
    print(json.dumps(dict(status='FINITE-EXACT complete k4 C/D owner universe',rows=rows),indent=2,sort_keys=True))


if __name__=='__main__':main()
