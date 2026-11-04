import itertools, collections
def ham(A,S):
    S=list(S)
    if not S: return 1
    return sum(1 for p in itertools.permutations(S) if all(A[p[i]][p[i+1]] for i in range(len(S)-1)))
def profile(A,S):
    S=list(S); n=len(S); res=collections.Counter()
    for labels in itertools.product(range(n),repeat=n):
        h=max(labels)+1
        if set(labels)!=set(range(h)): continue
        pr=1
        for b in range(h): pr*=ham(A,[S[i] for i in range(n) if labels[i]==b])
        res[h]+=pr
    return tuple(sorted(res.items()))
def tours(k):
    pairs=list(itertools.combinations(range(k),2))
    for bits in range(1<<len(pairs)):
        B=[[0]*k for _ in range(k)]
        for t,(i,j) in enumerate(pairs):
            if bits>>t&1: B[i][j]=1
            else: B[j][i]=1
        yield B
def ends(A):
    k=len(A); return {v for v in range(k) if sum(A[v])==k-1 or sum(A[v])==0}
for k in (4,5):
    allT=list(tours(k))
    sig={}
    for idx,B in enumerate(allT):
        sig[idx]=tuple(profile(B,S) for r in range(2,k+1) for S in itertools.combinations(range(k),r))
    groups=collections.defaultdict(list)
    for idx,s in sig.items(): groups[s].append(idx)
    stats=collections.Counter()
    for idx,A in enumerate(allT):
        if ham(A,range(k))==1: continue          # transitive: skip
        E=ends(A)
        if not E: stats['no end (strong-ish)']+=1; continue
        for jdx in groups[sig[idx]]:
            Bt=allT[jdx]
            common=E & ends(Bt)
            stats['common end' if common else 'NO common end']+=1
    print('order',k,dict(stats))
