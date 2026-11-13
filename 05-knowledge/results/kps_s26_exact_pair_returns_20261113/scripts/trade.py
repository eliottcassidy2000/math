# odd-cycle packing census of a tournament: directed odd cycles (up to rotation) by length, and disjoint packings
import itertools, sys
def decode(code,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1: out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def cycles(out):
    n=len(out); cyc=[]
    for L in range(3,n+1,2):
        for S in itertools.combinations(range(n),L):
            s0=S[0]
            for perm in itertools.permutations(S[1:]):
                seq=(s0,)+perm; ok=True
                for a,b in zip(seq,seq[1:]+(s0,)):
                    if not (out[a]>>b&1): ok=False; break
                if ok: cyc.append(frozenset(S) and (L,frozenset(S)))
    return cyc
def packing_profile(out):
    cyc=cycles(out); n=len(out)
    from collections import Counter
    single=Counter(L for L,S in cyc)
    pairs=Counter()
    for (L1,S1),(L2,S2) in itertools.combinations(cyc,2):
        if not (S1&S2): pairs[tuple(sorted((L1,L2)))]+=1
    triples=Counter()
    if n>=9:
        for a,b,c in itertools.combinations(cyc,3):
            if not(a[1]&b[1]) and not(a[1]&c[1]) and not(b[1]&c[1]): triples[tuple(sorted((a[0],b[0],c[0])))]+=1
    return dict(single), dict(pairs), dict(triples)
if __name__=='__main__':
    n=int(sys.argv[1])
    for c in map(int,sys.argv[2:]):
        out=decode(c,n); s,p,t=packing_profile(out)
        print(n,c,'out',out,'cycles',s,'disjoint pairs',p,'triples',t)
