#!/usr/bin/env python3
import hashlib,itertools,json
from pathlib import Path
from verify import require,subset_counts,transitive,tournament

def substitute(Q,blocks):
    offsets=[]; at=0
    for b in blocks:offsets.append(at);at+=len(b)
    rows=[]
    for i,b in enumerate(blocks):
        ext=sum(((1<<len(blocks[j]))-1)<<offsets[j] for j in range(len(Q)) if Q[i]>>j&1)
        rows.extend((r<<offsets[i])|ext for r in b)
    return tuple(rows)

def H(rows):return subset_counts(rows)[-1]
C3=(2,4,1)

def check(Q,sites):
    k=len(sites); index={v:i for i,v in enumerate(sites)}; deck=[]; values=[]
    for S in range(1<<k):
        blocks=[transitive(1 if (i not in index or S>>index[i]&1) else 3) for i in range(len(Q))]
        deck.append(H(substitute(Q,blocks)))
    for E in range(1<<k):
        blocks=[C3 if i in index and E>>index[i]&1 else transitive(3 if i in index else 1) for i in range(len(Q))]
        values.append(H(substitute(Q,blocks)))
        subsets=[s for s in range(1<<k) if s&E==s]
        exact=sum((1<<s.bit_count())*deck[s] for s in subsets)
        require(values[-1]==exact,'switch expansion failure')
        for r in range(1,k+2):
            truncated=sum((1<<s.bit_count())*deck[s] for s in subsets if s.bit_count()<r)
            require((values[-1]-truncated)%(1<<r)==0,'two-adic truncation failure')
    if k:
        types={((deck[1<<i]%4)-1)//2 for i in range(k)}
        offsets={0,2} if types=={0} else {0,6} if types=={1} else {0,2,4,6}
        require({v%8 for v in values}=={(deck[0]+x)%8 for x in offsets},'mod8 criterion failure')
    if k>=2:
        require(max(values)-min(values)>=3**k-1,'range floor failure')
        require(len(set(values))<(max(values)-min(values))//2+1,'interval gap failure')
    return deck,values

def main():
    families=0; outputs=0
    for n in range(1,4):
        for c in range(1<<(n*(n-1)//2)):
            Q=tournament(n,c)
            for mask in range(1<<n):
                sites=[i for i in range(n) if mask>>i&1]
                d,v=check(Q,sites);families+=1;outputs+=len(v)
    # source followed by C3; only source and one cyclic vertex selected.
    Q=(14,4,8,2)
    # all four nodes are TT3 in the example, not only the selected nodes.
    base=[transitive(3)]*4; sites=[0,1]; deck=[];values=[]
    for s in range(4):
        bs=[transitive(1) if i in sites and s>>sites.index(i)&1 else b for i,b in enumerate(base)]
        deck.append(H(substitute(Q,bs)))
        bs=[C3 if i in sites and s>>sites.index(i)&1 else b for i,b in enumerate(base)]
        values.append(H(substitute(Q,bs)))
    require(deck==[2721,2721,51,51],'mixed deck mismatch')
    require(values==[2721,8163,2823,8469],'mixed output mismatch')
    for e in range(4):require(values[e]==sum((1<<s.bit_count())*deck[s] for s in range(4) if s&e==s),'mixed expansion failure')
    hostile_base=H(substitute(C3,[transitive(3),transitive(1),transitive(1)]))
    hostile_changed=H(substitute(C3,[C3,transitive(1),transitive(1)]))
    require((hostile_base,hostile_changed)==(9,15),'hostile mismatch')
    require((hostile_changed-3*hostile_base)%8!=0,'false shortcut not refuted')
    print(json.dumps({'status':'finite controls PASS; proofs in ADDENDUM.md','complete_small_quotient_families':families,'literal_small_outputs':outputs,'mixed_12_vertex_deck':deck,'mixed_12_vertex_outputs':values,'false_mod8_shortcut_hostile':[hostile_base,hostile_changed],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2,sort_keys=True))
if __name__=='__main__':main()
