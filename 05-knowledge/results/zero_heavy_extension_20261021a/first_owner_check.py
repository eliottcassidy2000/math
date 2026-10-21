#!/usr/bin/env python3
"""Complete small first-owner fibres and independent matching incidence counts."""
from collections import Counter,defaultdict
from itertools import combinations
from math import comb
import json
import extension_check as X
Z=X.Z


def main():
    F=[1];I=[0];H=[0];incidence=[]
    for k in range(1,7):
        F.append(F[-1]*(2*k-1))
        I.append(F[k]-sum(I[j]*F[k-j] for j in range(1,k)))
        H.append((2*k-2)*F[k-1]-sum(H[j]*sum(F[a]*F[k-j-a] for a in range(k-j+1))
                                       for j in range(2,k)))
        if k<2:continue
        allmatch=list(Z.matchings(tuple(range(2*k))))
        strongmatch=[e for e in allmatch if Z.indecomp(e,2*k)]
        ha=sum(sum(b-a==2 for a,b in e) for e in allmatch)
        hs=sum(sum(b-a==2 for a,b in e) for e in strongmatch)
        ka=sum(comb(sum(b-a==1 for a,b in e),2) for e in allmatch)
        ks=sum(comb(sum(b-a==1 for a,b in e),2) for e in strongmatch)
        Z.need(ha==(2*k-2)*F[k-1] and hs==H[k],'span-two edge incidence')
        Z.need(ka==(k-1)*F[k-1],'two marked adjacent edges, all')
        Z.need(ks==(comb(2*k-4,2)*I[k-2] if k>=3 else 0),'two marked adjacent edges, strong')
        incidence.append(dict(k=k,all_H=ha,strong_H=hs,all_K=ka,strong_K=ks))
    layers=[]
    for k in (2,3):
        delta=2*k
        for strong in (False,True):
            fibres=defaultdict(list);accepted=Counter();rejected=Counter()
            recs=list(X.records(k,delta,strong))
            n=2*k;m=n-1 if strong else n+1
            for edges in Z.matchings(tuple(range(n))):
                if strong and not Z.indecomp(edges,n):continue
                for holes in combinations(range(m),2):
                    gaps=tuple(0 if i in holes else k for i in range(m))
                    recs.append(('D',edges,gaps,n))
            for kind,edges,gaps,n in recs:
                t,w=Z.model(n,edges,gaps,strong);clones=Z.clones(t)
                if kind in ('B','D'):
                    actual_gaps=(0,)+gaps+(0,) if strong else gaps
                    bad=any(b==a+1 and actual_gaps[b]==0 for a,b in edges)
                    Z.need(bool(clones)==bad,'complete omitted matched-gap rejection')
                    if bad:rejected[kind]+=1;continue
                Z.need(not clones,'literal retained pair-free class')
                if strong:Z.need(Z.strongly_connected(t),'literal strong first-owner class')
                Z.need(sum(w)==(2*k*k+k if strong else 2*k*k+3*k)-delta,'first-owner total mass')
                key,_=Z.canonical(t,w)
                fibres[key].append(kind);accepted[kind]+=1
            hist=Counter(tuple(sorted(v)) for v in fibres.values())
            Z.need(all(t in [('A',),('B',),('C',),('D',),('C','C','D')] for t in hist),
                   'the only nontrivial owner fibre is exactly C,C,D')
            row=incidence[k-2];h=row['strong_H' if strong else 'all_H'];kk=row['strong_K' if strong else 'all_K']
            Z.need(hist[('C','C','D')]==h,'every span-two matching edge creates exactly one owner fibre')
            E=I[k] if strong else F[k]
            J=(2*k-3)*I[k-1] if strong else F[k]
            P=sum(Z.indecomp(e,2*k-1) for e in Z.p3_graphs(k)) if strong else (k-1)*F[k]
            naive=(E*(comb(m+delta-1,delta)-m*comb(m+delta-k-1,delta-k))
                  +(m*E-J)*comb(m+delta-k-2,delta-k)
                  +P*comb(m+delta-k-3,delta-k-1))
            expected=naive+kk-2*h
            Z.need(len(fibres)==expected,'corrected exact class count')
            layers.append(dict(k=k,strong=strong,classes=len(fibres),raw_records=sum(accepted.values()),
                               accepted=dict(accepted),rejected=dict(rejected),
                               fibre_histogram={','.join(t):v for t,v in hist.items()},
                               corridor_extrapolation=naive,correction=kk-2*h))
    print(json.dumps(dict(status='FINITE-EXACT first-owner fibres and incidence controls',
                         matching_incidence=incidence,layers=layers),indent=2,sort_keys=True))


if __name__=='__main__':main()
