#!/usr/bin/env python3
from verify import *
import random

def sfadd(*ps):
    out={}
    for p in ps:
        for m,v in p.items():out[m]=out.get(m,0)+v
    return {m:v for m,v in out.items() if v}
def sfscale(p,c):return {m:v*c for m,v in p.items() if v*c}
def sfmul(p,q):
    out={}
    for m,v in p.items():
        for n,w in q.items():
            if not m&n:out[m|n]=out.get(m|n,0)+v*w
    return {m:v for m,v in out.items() if v}
def sfnegvars(p):return {m:v*(-1)**m.bit_count() for m,v in p.items()}
def literal_kernel(a):
    n=len(a);K=[[{} for _ in a] for _ in a]
    for start in range(n):
        def rec(last,used):
            K[start][last][used]=K[start][last].get(used,0)+1
            for v in range(n):
                if not used>>v&1 and a[last]>>v&1:rec(v,used|1<<v)
        rec(start,1<<start)
    return K

def kernel_checks(a):
    n=len(a);K=literal_kernel(a)
    starts=[sfadd(*row) for row in K]
    ends=[sfadd(*(K[i][j] for i in range(n))) for j in range(n)]
    F=sfadd({0:1},*starts)
    require(sfmul(F,sfnegvars(F))=={0:1},'unrooted squarefree reciprocity')
    h=h_subsets(a)
    require(F=={s:v for s,v in enumerate(h) if v},'literal kernel versus subset DP')
    for i in range(n):
        require(sfmul(F,sfnegvars(ends[i]))==sfscale(starts[i],-1),'full endpoint polynomial identity')
        for j in range(n):
            require(sfmul(F,sfadd(K[i][j],sfnegvars(K[j][i])))==sfmul(starts[i],ends[j]),'full squarefree matrix rank-one relation')
    return n*n

def run():
    rooted=0;matrix_entries=0;witness=None;hist={}
    for n in range(1,6):
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);h=h_subsets(a);K=literal_kernel(a);full=(1<<n)-1
            for r in range(n):
                end=sfadd(*(K[i][r] for i in range(n)))
                d=sum((-1)**used.bit_count()*v*h[full^used] for used,v in end.items())
                start=sum(K[r][j].get(full,0) for j in range(n))
                require(d==-start,'actual alternating endpoint response')
                hist[d]=hist.get(d,0)+1;rooted+=1
                if start==5 and witness is None:witness={'rows':a,'root':r,'H':h[-1],'Hamilton_paths_starting_at_root':start,'alternating_response':d}
            if n<=4:matrix_entries+=kernel_checks(a)
    rng=random.Random(1971147)
    for n,rounds in ((5,32),(6,20)):
        for _ in range(rounds):matrix_entries+=kernel_checks(tournament(n,rng.getrandbits(n*(n-1)//2)))
    require(rooted==5405 and witness is not None,'complete labelled census and nonunit witness')
    print(json.dumps({'status':'PASS','rooted_tournaments':rooted,'full_matrix_entries_checked':matrix_entries,'signed_response_histogram':hist,'nonunit_response_witness':witness},sort_keys=True,indent=2))
if __name__=='__main__':run()
