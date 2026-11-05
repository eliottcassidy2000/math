# Check: H(TT_n + x_w) = g(val(w)), B_x = 1 + val(w), e_x, s_x by brute force for n <= 7.
import itertools, sys
def ham_stats(n, out):
    # out[i] = bitmask of out-neighbours; returns h, list e_v, list s_v, list B_v
    full=(1<<n)-1
    # dp over subsets: count paths in subset S ending at v
    end=[[0]*n for _ in range(1<<n)]
    for v in range(n): end[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            c=end[S][v]
            if not c: continue
            nb=out[v]&~S
            while nb:
                w=(nb&-nb).bit_length()-1; nb&=nb-1
                end[S|(1<<w)][w]+=c
    hs=[sum(end[S]) if S else 1 for S in range(1<<n)]
    # start counts: paths in subset S starting at v = paths in converse ending at v; compute directly
    st=[[0]*n for _ in range(1<<n)]
    inn=[0]*n
    for v in range(n):
        for u in range(n):
            if out[u]>>v&1: inn[v]|=1<<u
    for v in range(n): st[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            c=st[S][v]
            if not c: continue
            nb=inn[v]&~S
            while nb:
                w=(nb&-nb).bit_length()-1; nb&=nb-1
                st[S|(1<<w)][w]+=c
    h=hs[full]
    e=[end[full][v] for v in range(n)]
    s=[st[full][v] for v in range(n)]
    B=[]
    for t in range(n):
        tot=0
        for S in range(0,1<<n):
            if S>>t&1: continue
            tot+=hs[S]*st[full&~S][t]
        B.append(tot)
    return h,e,s,B
def g(v):
    if v==0: return 1
    if v&1: return g(v>>1)
    u=v>>1
    return g(u)+u+(u&1)
bad=0
for n in range(1,8):
    for w in itertools.product([0,1],repeat=n):
        # chain 0..n-1, i->j for i<j; x = n; x -> i iff w[i]=1
        out=[0]*(n+1)
        for i in range(n):
            for j in range(i+1,n): out[i]|=1<<j
            if w[i]: out[n]|=1<<i
            else: out[i]|=1<<n
        h,e,s,B=ham_stats(n+1,out)
        val=int(''.join(map(str,w)),2)
        if h!=g(val) or B[n]!=1+val or e[n]!=(1 if (val==0 or w[-1]==0) else 0):
            bad+=1
            if bad<5: print('BAD',w,h,g(val),B[n],1+val,e[n])
print('checked one-apex identities n<=7, failures',bad)
