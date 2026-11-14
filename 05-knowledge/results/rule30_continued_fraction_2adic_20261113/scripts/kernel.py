# 2-kernel of the Rule 30 center column: are the subsequences c_{2^k n + r} (0<=r<2^k, k<=K) pairwise distinct
# on the computed prefix? (A 2-automatic sequence has a finite 2-kernel.)
import sys
raw=open(sys.argv[1],'rb').read(); N=len(raw); K=int(sys.argv[2])
elems=[]
for k in range(K+1):
    step=1<<k
    for r in range(step):
        seq=raw[r::step]
        elems.append((k,r,seq))
# pairwise distinctness: two kernel elements are certified distinct if they differ at some index < min(len)
import collections
minlen=min(len(e[2]) for e in elems)
groups=collections.defaultdict(list)
for k,r,seq in elems: groups[seq[:minlen]].append((k,r,seq))
undecided=0; total=len(elems)
for key,g in groups.items():
    if len(g)>1:
        # try longer comparisons
        for i in range(len(g)):
            for j in range(i+1,len(g)):
                a=g[i][2]; b=g[j][2]; m=min(len(a),len(b))
                if a[:m]==b[:m]: undecided+=1
print(f"N={N}, K={K}: {total} kernel elements; shortest compared length {minlen}; groups sharing that prefix: {sum(1 for g in groups.values() if len(g)>1)}; pairs not separated on the prefix: {undecided}")
