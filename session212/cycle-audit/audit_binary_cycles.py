#!/usr/bin/env python3
"""Independent Session212 audit; never imports the researcher's builder/counter.
Directed cycles are modulo rotation, not reversal. No graph-library dependency.
"""
import hashlib, itertools, json
from functools import lru_cache
from pathlib import Path


def require(ok, message="independent audit check failed"):
    if not ok:
        raise RuntimeError(message)


def binary_matrix(h, m=None):
    if h < 1:
        raise ValueError('positive h required')
    if m is None:
        m=(h-1).bit_length()
    if m<0 or h>2**m:
        raise ValueError('insufficient digits')
    n=m+4; z=n-2; x=n-1
    adj=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if j==x:
                winner=x if i==1 else i
            elif j==z:
                winner=z if i==0 or (i>=2 and not ((h-1)>>(i-2)&1)) else i
            else:
                winner=i
            loser=j if winner==i else i
            adj[winner][loser]=True
    require(all(not adj[i][i] for i in range(n)), "loop")
    require(all(adj[i][j]!=adj[j][i] for i in range(n) for j in range(i)), "pair orientation")
    return adj


def cycles(adj, root=0):
    """Top-down suffix/remaining-set DP, distinct from forward author DP."""
    n=len(adj); full=(1<<n)-1
    @lru_cache(None)
    def tail(v, remaining):
        if not remaining:
            return int(adj[v][root])
        out=0
        for w in range(n):
            if remaining>>w&1 and adj[v][w]:
                out+=tail(w, remaining^(1<<w))
        return out
    return tail(root, full^(1<<root))


def brute_cycles(adj,root=0):
    n=len(adj); total=0
    for perm in itertools.permutations(i for i in range(n) if i!=root):
        seq=(root,)+perm+(root,)
        total+=all(adj[a][b] for a,b in zip(seq,seq[1:]))
    return total


def start_paths(adj, root):
    n=len(adj)
    @lru_cache(None)
    def rec(v,remaining):
        if not remaining:return 1
        return sum(rec(w,remaining^(1<<w)) for w in range(n) if remaining>>w&1 and adj[v][w])
    return rec(root,((1<<n)-1)^(1<<root))


def strong(adj):
    n=len(adj)
    for root in range(n):
        seen={root}; stack=[root]
        while stack:
            v=stack.pop()
            for w in range(n):
                if adj[v][w] and w not in seen:
                    seen.add(w);stack.append(w)
        if len(seen)!=n:return False
    return True


def masks(adj):return [sum(1<<j for j,b in enumerate(row) if b) for row in adj]


def main():
    report={'status':'PASS','cycle_convention':'directed Hamiltonian cycles modulo rotation, never reversal', 'checks':{}}
    count=0
    # Exact all-family test, including many leading-zero padded constructions.
    for m in range(9):
        for h in range(1,2**m+1):
            a=binary_matrix(h,m)
            require(cycles(a)==h, (h,m,masks(a)))
            require(strong(a), "strongness")
            require(start_paths([row[:-1] for row in a[:-1]],1)==h, "endpoint count")
            count+=1
    report['checks']['all_digits_m_0_through_8']=count
    # Separately implemented exhaustive rooted permutation oracle.
    count=0
    for m in range(6):
        for h in range(1,2**m+1):
            a=binary_matrix(h,m)
            require(brute_cycles(a)==h, "permutation cycle count")
            count+=1
    report['checks']['permutation_oracle_m_0_through_5']=count
    # Every anchor, to explicitly exclude rooted/unrooted normalization errors.
    anchors=0
    for h in (1,2,3,4,7,21,31,32,127,128,255,256):
        a=binary_matrix(h)
        for r in range(len(a)):
            require(cycles(a,r)==h, "anchor normalization")
            anchors+=1
    report['checks']['all_anchors_selected_targets']=anchors
    # The one-cycle baseline is unchanged when all optional bits are zero.
    for m in range(14):
        a=binary_matrix(1,m); require(cycles(a)==1 and strong(a), "zero-digit padding")
    report['checks']['zero_digit_padding_m_0_through_13']=14
    # Hostile/control tests independent of the family.
    for n in range(3,10):
        tr=[[i<j for j in range(n)] for i in range(n)]
        require(cycles(tr)==0 and not strong(tr), "transitive control")
    report['checks']['transitive_zero_controls_n_3_through_9']=7
    invalid=0
    for h,m in ((0,None),(-1,None),(2,0),(3,1),(1,-1)):
        try:binary_matrix(h,m)
        except ValueError:invalid+=1
        else:raise AssertionError((h,m))
    report['checks']['invalid_inputs_rejected']=invalid
    report['witnesses']=[{'h':h,'n':len(binary_matrix(h)),'outmasks':masks(binary_matrix(h))} for h in (1,2,3,7,21,255,256)]
    report['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text=json.dumps(report,indent=2)+'\n'
    Path(__file__).with_name('binary_audit.json').write_text(text)
    print(text,end='')

if __name__=='__main__':main()
