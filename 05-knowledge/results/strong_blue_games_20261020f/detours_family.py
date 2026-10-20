"""Literal controls for the two-block repair and the two-chord family.

R reverses precisely the displayed Hamiltonian path. Certificates for k=2,3
include their final quotient isomorphism explicitly. All checks survive -O.
"""
from itertools import permutations
import json
from pathlib import Path
from detours_probe import G
from detours_square import reverse, opposite

HERE=Path(__file__).resolve().parent


def family(k):
    n=2*k+1
    a=[sum(1<<v for v in range(u+1,n)) for u in range(n)]
    for u,v in ((0,k),(k,2*k)):
        a[u]^=1<<v;a[v]^=1<<u
    return tuple(a)


def repaired_transitive_path(k):
    if k==4:return (0,2,1,3)
    if k==5:return (0,3,2,1,4)
    if k==6:return (0,2,1,4,3,5)
    G.need(k>=7,'repair requires block size >=4')
    m=k-3
    return repaired_transitive_path(m)+(m+1,m,m+2)


def run_r(a,words):
    sums=[]
    for p in words:
        G.need(G.is_path(a,p),'R word must be a Hamiltonian path')
        G.need(not G.blue(a,p),'every certificate move must be black')
        sums.append(a[p[0]].bit_count()+a[p[-1]].bit_count())
        a=reverse(a,p)
    return a,sums


def isomap(a,b):
    for p in permutations(range(len(a))):
        if all(a[u].bit_count()==b[p[u]].bit_count() for u in range(len(a))) and all(
                bool(a[u]>>v&1)==bool(b[p[u]]>>p[v]&1)
                for u in range(len(a)) for v in range(len(a))):
            return p
    raise RuntimeError('no certificate isomorphism')


def main():
    checks=[]
    for k in range(4,101):
        a=family(k);n=len(a);P=tuple(range(n));h=repaired_transitive_path(k)
        tt=tuple(sum(1<<v for v in range(u+1,k)) for u in range(k))
        G.need(G.is_path(reverse(tt,tuple(range(k))),h),'internal endpoint repair')
        L=h;R=tuple(k+1+x for x in h)
        words=[tuple(range(k))+tuple(range(k+1,n))+(k,),L+(k,)+R,(k,)+R[::-1]+L[::-1]]
        final,sums=run_r(a,words)
        G.need(final==reverse(a,P),'exact labelled target')
        G.need(sums==[3*k-1,2*k-1,3*k-2],'symbolic score formula')
        G.need(G.blue(a,P),'the bypassed path must be blue')
        checks.append(dict(k=k,n=n,endpoint_sums=sums))
    small=[]
    for k,left,right,converse in [
            (2,[(1,3,4,2,0)],[(0,4,3,2,1)],False),
            (3,[(0,1,2,4,5,6,3),(0,2,1,3,4,6,5)],[(3,2,1,0,6,5,4)],True)]:
        a=family(k);b=reverse(a,tuple(range(len(a))))
        x,xs=run_r(a,left);y,ys=run_r(b,right)
        if converse:y=opposite(y)
        p=isomap(x,y)
        small.append(dict(k=k,source=a,target=b,source_R_paths=left,
                          target_R_paths=right,converse_target_neighbor=converse,
                          isomorphism_from_source_neighbor=p,
                          source_endpoint_sums=xs,target_endpoint_sums=ys))
    out=dict(status='FINITE-EXACT controls for the all-order proof',
             block_sizes=list(range(4,101)),checks=checks,small_boundaries=small)
    (HERE/'detours_family.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(block_sizes='4..100',literal_family_orders='9..201',small_boundaries=small),indent=2))


if __name__=='__main__':main()
