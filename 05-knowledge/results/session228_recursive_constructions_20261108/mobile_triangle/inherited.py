#!/usr/bin/env python3
"""Exact standard-library controls for the scoped adaptive-diamond report."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
import json

OPS={
'O':((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1)),
'Z':((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1)),
'J':((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0)),
'L':((1,-1,3,-2),(0,0,1,-1),(2,-1,4,-1),(1,-1,1,0)),
'C':((3,-3,9,-6),(1,-1,2,-2),(3,-2,11,-6),(1,-1,3,-2))}
E0=((1,0,0),(0,0,0),(0,1,0),(0,0,1))
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))

def require(p,msg):
    if not p: raise RuntimeError(msg)
def mv(M,x):return tuple(sum(a*b for a,b in zip(row,x)) for row in M)
def mm(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))) for i in range(len(A)))
def word(w,B=I):
    for c in w:B=mm(OPS[c],B)
    return B
def det3(M):
    a,b,c=M
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def tournaments(n):
    pairs=list(combinations(range(n),2))
    for bits in range(1<<len(pairs)):
        rows=[0]*n
        for k,(i,j) in enumerate(pairs):
            a,b=(i,j) if bits>>k&1 else (j,i)
            rows[a]|=1<<b
        yield rows

def profile(rows,r):
    n=len(rows);size=1<<n
    dp=[[0]*n for _ in range(size)];start=[[0]*n for _ in range(size)]
    for i in range(n):dp[1<<i][i]=1
    start[1<<r][r]=1
    for mask in range(1,size):
        for i in range(n):
            if not (dp[mask][i] or start[mask][i]):continue
            todo=rows[i]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;j=bit.bit_length()-1
                dp[mask|bit][j]+=dp[mask][i]
                start[mask|bit][j]+=start[mask][i]
    h=list(map(sum,dp));h[0]=1;s=list(map(sum,start));full=size-1
    return h[full],dp[full][r],sum(s[m]*h[full^m] for m in range(size)),h[full^(1<<r)]

def graft(rows,r,b,br=0):
    n=len(rows);newmask=((1<<len(b))-1)<<n
    out=[v|(newmask if i!=r else 0) for i,v in enumerate(rows)]
    out.extend((v<<n)|(1<<r) for v in b)
    return out,n+br

def graph_word(rows,r,w):
    rows=rows[:]
    for c in w:
        if c=='O':
            bit=1<<len(rows);rows=[v|bit for v in rows]+[0]
        elif c=='Z':
            rows,unused=graft(rows,r,[0])
        else:rows,r=graft(rows,r,{'J':[0],'L':[2,0],'C':[4,1,2]}[c])
    return rows,r

def solve(A,b):
    a=[[Fraction(v) for v in row]+[Fraction(t)] for row,t in zip(A,b)]
    for j in range(len(a)):
        k=next(i for i in range(j,len(a)) if a[i][j]);a[j],a[k]=a[k],a[j]
        v=a[j][j];a[j]=[z/v for z in a[j]]
        for i in range(len(a)):
            if i!=j:
                v=a[i][j];a[i]=[z-v*y for z,y in zip(a[i],a[j])]
    return tuple(row[-1] for row in a)

def main():
    A=word('JOLO',E0);B=word('JJLO',E0)
    require(A==((5,7,0),(0,0,0),(18,20,0),(3,3,0)),'high switch matrix')
    require(B==((5,5,0),(0,0,0),(18,20,0),(3,3,0)),'low switch matrix')
    require(word('JOZJO',E0)==((5,7,0),(0,0,0),(12,14,0),(3,3,0)),'old high switch')
    require(word('JJJJO',E0)==((5,6,0),(0,0,0),(12,14,0),(3,3,0)),'old low switch')
    for D,U,d in product(range(-3,4),repeat=3):
        y=mv(A,(D,U,d));z=mv(B,(D,U,d))
        require(y[1:]==z[1:],'sidecars')
        require(y[0]-z[0]==2*U,'scalar switch')
        require(20*y[3]-3*y[2]==6*D,'predecessor D')
        require(y[2]-6*y[3]==2*U,'predecessor U')
        z2=mv(word('JOZJO',E0),(D,U,d))
        require(z2[2]-4*z2[3]==2*U,'second predecessor U')
        require(14*z2[3]-3*z2[2]==6*D,'second predecessor D')
        require(z2[2]-6*z2[3]==-6*D-4*U,'second sidecar sign identity')
    # All positive-depth mixtures of the four switches through depth six.
    # The sign identifies the final pair; both already verified inverse maps
    # then limit every physical sidecar fiber to two scalar outputs.
    switch_words=('JOLO','JJLO','JOZJO','JJJJO')
    switch_matrices=[word(w) for w in switch_words]
    frontier={(1,0,2,1)};mixed_outputs=set();mixed_fibers=defaultdict(set)
    for depth in range(1,7):
        nxt=set()
        for x in frontier:
            for index,M in enumerate(switch_matrices):
                y=mv(M,x);sidecar_sign=y[2]-6*y[3]
                require((sidecar_sign>0)==(index<2),'disjoint physical sidecar chambers')
                require(sidecar_sign!=0,'no physical sidecar boundary')
                nxt.add(y);mixed_outputs.add(y);mixed_fibers[y[1:]].add(y[0])
        frontier=nxt
    require(max(map(len,mixed_fibers.values()))<=2,'mixed switch fiber bound')
    # Literal controls include every labelled root-at-zero source through order4.
    sources=[(r,0) for n in range(1,5) for r in tournaments(n)]
    literal_moves=literal_switches=0
    for rows,r in sources:
        x=profile(rows,r)
        for c in OPS:
            out,root=graph_word(rows,r,c)
            require(profile(out,root)==mv(OPS[c],x),'literal move '+c)
            literal_moves+=1
        if x[1]==0:
            for w in ('JOLO','JJLO','JOZJO','JJJJO'):
                out,root=graph_word(rows,r,w)
                require(profile(out,root)==mv(word(w),x),'literal diamond '+w)
                literal_switches+=1
    # Recover the first graft column by literal evaluations on a physical basis.
    basis=[([0],0),([2,0],0),([4,1,2],0),([6,4,0],0)]
    columns=[profile(r,k) for r,k in basis]
    weights=solve(tuple(zip(*columns)),(1,0,0,0))
    first_columns=0
    for b,br in sources:
        results=[]
        for a,ar in basis:
            out,root=graft(a,ar,b,br);results.append(profile(out,root))
        recovered=tuple(sum(weights[j]*results[j][i] for j in range(4)) for i in range(4))
        require(recovered==profile(b,br),'physical first graft column')
        first_columns+=1
    # Exact word controls for the sensitivity theorem, on all words throughlength5.
    sensitivity=0
    for k in range(6):
        for chars in product(OPS,repeat=k):
            w=''.join(chars);col=mv(word(w),(1,0,0,0))
            idx=next((i for i,c in enumerate(w) if c in 'JLC'),None)
            if idx is None:require(col==(1,0,0,0),'OZ fixed column')
            else:
                rows={'J':[0],'L':[2,0],'C':[4,1,2]}[w[idx]]
                x=profile(rows,0)
                expected=mv(word(w[idx+1:]),x)
                require(col==expected and col[2]>=1 and col[3]>=1,'graft sensitivity')
            sensitivity+=1
    # Complete bounded operator universe. We deduplicate terminal-O matrices,
    # so counts refer to distinct transfers, not numbers of word witnesses.
    level={E0:''};allmat={E0:''};returns={};levels=[]
    for depth in range(8):
        levels.append(len(level))
        for M,w in level.items():returns.setdefault(mm(OPS['O'],M),w+'O')
        if depth<7:
            nxt={}
            for M,w in level.items():
                for c,X in OPS.items():
                    Y=mm(X,M)
                    if Y not in allmat:nxt.setdefault(Y,w+c)
            allmat.update(nxt);level=nxt
    buckets=defaultdict(list)
    for M,w in returns.items():buckets[(M[2],M[3])].append((M[0],w))
    nontrivial=[]
    for key,values in buckets.items():
        for (a,wa),(b,wb) in combinations(values,2):
            diff=tuple(v-u for u,v in zip(a,b))
            require(any(diff),'distinct returned matrices')
            nontrivial.append((wa,wb,key,a,b))
            require(not (diff[0]==0 and diff[1]==0),'no d-only switch within stated finite universe')
            require(det3((a,*key))==det3((b,*key))==0,'all bounded collision branches rank<=2')
    sample={w:list(mv(word(w),(1,0,2,1))) for w in ('JOLO','JJLO','JOLOZ','JJLOZ')}
    result={'status':'PASS; scoped identities and finite controls only',
      'literal_one_move_checks':literal_moves,'literal_diamond_checks':literal_switches,
      'mixed_switch_depth':6,'mixed_switch_distinct_outputs':len(mixed_outputs),
      'mixed_switch_maximum_scalar_fiber':max(map(len,mixed_fibers.values())),
      'physical_basis_weights':[str(v) for v in weights],
      'general_gadget_first_columns_checked':first_columns,
      'sensitivity_words_through_length5':sensitivity,
      'new_restricted_transfers_by_minimum_prefix_length':levels,
      'distinct_terminal_O_transfers_prefix_length_at_most7':len(returns),
      'sidecar_buckets_with_multiple_scalar_rows':sum(len(v)>1 for v in buckets.values()),
      'unordered_distinct_scalar_switch_pairs':len(nontrivial),
      'd_only_switch_pairs':0,'rank_three_switch_branches':0,
      'sample_seed_outputs':sample,
      'open':'No universal radix synthesis, recursive interval coverage, or Chan-Pak proof.'}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
