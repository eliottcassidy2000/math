"""
Extend the Klein-four iso-class algebra: complement structure, and the RECURSIVE
score-residual COMPRESSION (score sequence = compressed code; residual = fine structure).

tournament: bit per pair (i<j), bit=1 iff i->j.  n<=6.
"""
from itertools import permutations, combinations

def pairs(n): return list(combinations(range(n),2))

def scores(n, bits):
    pr=pairs(n); s=[0]*n
    for k,(i,j) in enumerate(pr):
        if bits[k]: s[i]+=1
        else: s[j]+=1
    return tuple(sorted(s))

def canon(n, bits):
    pr=pairs(n); idx={p:k for k,p in enumerate(pr)}
    best=None
    for perm in permutations(range(n)):
        out=[]
        for (a,b) in pr:  # new labels a<b
            oi,oj=perm[a],perm[b]
            if oi<oj: out.append(bits[idx[(oi,oj)]])
            else: out.append(1-bits[idx[(oj,oi)]])
        t=tuple(out)
        if best is None or t<best: best=t
    return best

def complement(bits): return tuple(1-b for b in bits)

def all_iso(n):
    pr=pairs(n); m=len(pr); seen={}
    for x in range(1<<m):
        bits=tuple((x>>k)&1 for k in range(m))
        c=canon(n,bits)
        if c not in seen: seen[c]=scores(n,bits)
    return seen   # canon -> score seq

print("n | #iso(A000568) | #score-seqs | residual (iso classes per score seq) | SC | NS-pairs")
for n in [3,4,5,6]:
    iso=all_iso(n)
    by_score={}
    for c,sc in iso.items(): by_score.setdefault(sc,[]).append(c)
    nres=sorted(set(len(v) for v in by_score.values()))
    # complement: SC = canon == canon(complement)
    sc=0; pairset=set()
    for c in iso:
        cc=canon(n,complement(c))
        if cc==c: sc+=1
        else: pairset.add(frozenset([c,cc]))
    print(f"{n} | {len(iso):>3} | {len(by_score):>3} | residual sizes {nres}, max={max(len(v) for v in by_score.values())} | SC={sc} | NSpairs={len(pairset)}")

# n=4: the Klein four + complement automorphism
print("\nn=4 Klein-four + complement automorphism:")
NM={(0,1,2,3):'T',(0,2,2,2):'+',(1,1,1,3):'-',(1,1,2,2):'S'}
iso4=all_iso(4)
for c,sc in iso4.items():
    cc=canon(4,complement(c))
    print(f"  class {NM[sc]} (score {sc}): complement -> {NM[scores(4,cc)]}  {'[SC/fixed]' if cc==c else '[swapped]'}")
print("  => complement = the (+ <-> -) automorphism, fixing {T,S}. The SC subgroup {T,S}=<S> is the")
print("     +1 eigenspace of complement; the NS pair {+,-} is the -1 eigenspace. The Klein-four IS")
print("     the complement-involution's eigenspace decomposition, made into a group (S = + * -).")

# RECURSIVE compression: how small is the 'residual' beyond the score sequence?
print("\nRECURSIVE COMPRESSION view: score sequence as the compressed code")
for n in [4,5,6]:
    iso=all_iso(n)
    by_score={}
    for c,sc in iso.items(): by_score.setdefault(sc,[]).append(c)
    import math
    res_bits=sum(max(0,math.ceil(math.log2(len(v)))) for v in by_score.values())
    score_bits=math.ceil(math.log2(len(by_score)))
    full=len(pairs(n))-(n-1)  # tiling bits C(n-1,2)
    print(f"  n={n}: tiling bits={full}; score-code ~{score_bits} bits selects 1 of {len(by_score)} score seqs;")
    print(f"        residual: {sum(len(v)-1 for v in by_score.values())} classes need a tie-break ("
          f"max {max(len(v) for v in by_score.values())} per score). Score seq alone resolves "
          f"{sum(1 for v in by_score.values() if len(v)==1)}/{len(by_score)} score seqs uniquely.")
