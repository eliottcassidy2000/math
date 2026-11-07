#!/usr/bin/env python3
"""Independent checks of marked growth, strong cores and two-endpoint support."""
from functools import lru_cache
from pathlib import Path
import hashlib,json
from audit_binary_cycles import binary_matrix,cycles,strong,require
from audit_packet import comb,gatefree


def fixed_end_paths(a,start,end):
    n=len(a)
    @lru_cache(None)
    def rec(v,left):
        if not left:return int(v==end)
        if v==end:return 0
        return sum(rec(w,left^(1<<w)) for w in range(n) if left>>w&1 and a[v][w])
    return rec(start,((1<<n)-1)^(1<<start))


def insert_bit(a,bit):
    old=len(a);z=old-2;x=old-1
    new=[[False]*(old+1) for _ in range(old+1)]
    place=lambda v:v if v<2 else v+1
    for u in range(old):
        for v in range(old):new[place(u)][place(v)]=a[u][v]
    for v in range(old):
        if v==x:out=True
        elif v==z:out=bool(bit)
        else:out=v>=2
        new[2][place(v)]=out;new[place(v)][2]=not out
    return new


def main():
    growth=0
    for m in range(1,7):
        for q in range(2**m):
            for tau in (0,1):
                a=comb(m,q,tau,0)
                for bit in (0,1):
                    new=insert_bit(a,bit)
                    require(new==comb(m+1,2*q+bit,tau,0),'marked induced embedding')
                    require(cycles(new)==2*cycles(a)-tau+bit,'binary growth')
                    require(gatefree(new),'growth gate')
                    growth+=1
    for h in range(1,257):
        a=binary_matrix(h);core=[r[:-1] for r in a[:-1]]
        require(strong(core),'minimal-width core not strong')
        require(fixed_end_paths(a,1,len(a)-1)==h,'both endpoints count')
    data={'status':'PASS','marked_growth_checks':growth,'strong_minimal_width_cores':256,'fixed_start_end_checks':256,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(data,indent=2))
if __name__=='__main__':main()
