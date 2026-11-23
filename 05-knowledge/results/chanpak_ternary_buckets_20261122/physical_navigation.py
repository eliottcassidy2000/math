"""Reconstruct the moving ternary forks, including the alignments of the growing blocks."""
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"chanpak_decimal_kernel_20261122"))
from e15_pair import flipped, star_rows, tagged_data
from fermat_receiver import address_height

# Each chart lists the rows of the ordinary tournament, the apex out-mask,
# an odd word, the suffix mode (0 for 2^s*w, 1 for 2^s*w-1), and a shift of s.
EXAMPLES=[
    (((2,4,1),1,7,0,0), ((0,1,11,19,71,31,47),11,1,1,1),
     ((2,3),(4,5)),42,-13),
    (((8,5,9,18,7),11,5,0,0), ((16,9,67,21,38,15,59),68,1,1,1),
     ((2,4),(0,2)),280,7),
    (((10,40,19,148,99,77,143,55),9,1,1,0), ((0,13,41,145,39,75,31,119),50,3,0,0),
     ((3,6),(1,4)),1019,219),
]


def sampled(spec, s):
    rows,mask,w,mode,shift=spec
    y=(w<<(s+shift))-mode
    host=star_rows(rows,mask,0,y)
    n=len(rows)
    # The first s right-hand bits make a uniform transitive module.
    run=range(n,n+s)
    before=[]
    after=[]
    for i in range(len(host)):
        if i in run:
            continue
        incidence={(host[i]>>j)&1 for j in run}
        assert len(incidence)==1
        (before if incidence.pop() else after).append(i)
    return host,before,after


def pad(host,left,right,before,after):
    """Append labelled source and sink vertices so the two outside parts match in size."""
    host=list(host)
    for _ in range(before-len(left)):
        n=len(host)
        host.append((1<<n)-1)
        left.append(n)
    for _ in range(after-len(right)):
        n=len(host)
        host=[r|(1<<n) for r in host]+[0]
        right.append(n)
    return tuple(host)


def align(parent,child,s):
    p,pl,pr=sampled(parent,s)
    c,cl,cr=sampled(child,s)
    before,after=max(len(pl),len(cl)),max(len(pr),len(cr))
    p=pad(p,pl,pr,before,after)
    c=pad(c,cl,cr,before,after)
    L,R=len(pl),len(pr)
    pv,cv=pl+pr,cl+cr
    # Minimize the ordinary orientation differences in this module-preserving gauge.
    best=[len(pv)**2,()]
    def choose(t,used,perm,cost):
        if cost>=best[0]:
            return
        if t==L+R:
            best[:]=cost,tuple(perm)
            return
        domain=range(L) if t<L else range(L,L+R)
        for j in domain:
            if used>>j&1:
                continue
            extra=sum(((p[pv[u]]>>pv[t])^(c[cv[v]]>>cv[j]))&1
                      for u,v in enumerate(perm))
            choose(t+1,used|1<<j,perm+[j],cost+extra)
    choose(0,0,[],0)
    return len(p),best[0],best[1]


def direct(host):
    # Independent Hamilton-path DP for the initial finite instances.
    paths=[[0]*len(host) for _ in range(1<<len(host))]
    for i in range(len(host)):
        paths[1<<i][i]=1
    for mask,ends in enumerate(paths):
        for j,count in enumerate(ends):
            if not count: continue
            free=host[j]&~mask
            while free:
                bit=free&-free
                paths[mask|bit][bit.bit_length()-1]+=count
                free-=bit
    return sum(paths[-1])


def check():
    for parent,child,edges,L,K in EXAMPLES:
        for s in (1,2):
            ph=sampled(parent,s)[0]
            ch=sampled(child,s)[0]
            h=L*(1<<s)+K
            assert direct(ph)==h and direct(ch)==3*h
            assert [direct(flipped(ch,*edge)) for edge in edges]==[3*h-2,3*h+2]
        for spec, edits, mult in ((parent,(),1),(child,edges,3)):
            rows,mask,w,mode,shift=spec
            for offset,edge in [(0,None)]+list(zip((-2,2),edits)):
                trows=flipped(rows,*edge) if edge else rows
                data=tagged_data(trows,mask)
                # Each read is affine in 2^s by the split identity.
                for s in (1,2,17,61):
                    address=(w<<(s+shift))-mode
                    assert address_height(data,0,address)==mult*(L*(1<<s)+K)+offset
        order,dist,perm=align(parent,child,1)
        assert align(parent,child,2)[1:]==(dist,perm)
        print(L,K,"common order s+",order-1,"edits",dist,"outside permutation",perm,
              "at s=1",L*2+K,3*(L*2+K)-2,3*(L*2+K),3*(L*2+K)+2,flush=True)


if __name__=="__main__":
    check()
