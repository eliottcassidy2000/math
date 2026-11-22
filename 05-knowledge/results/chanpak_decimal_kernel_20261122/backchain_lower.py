"""Exact exponent-61 backchain lower and one-arc contacts at the last two gaps."""
from argparse import ArgumentParser
from functools import lru_cache
from itertools import permutations
from math import prod
from e15_pair import flipped, strong, subpaths

WORDS = ((2,), (2,), (1,0,0,1),
         (1,0,0,0,1,0,0,4,0,0,0,2),
         (1,0,0,0,0,1,0,0,0,2,0,0,0,1,0,1))
PRIMES = (2,3,5,7,11,13,17,19,23,29,31,41,47,59,71)


def height(word):
    d,u,v=0,0,1
    for width in word:
        if width:
            d,u,v=(1<<width)*(d+2*u+v),0,d+u+v
        else:
            d,u,v=d+2*u,v,d+u
    return d+u+v


def rows_for(word):
    bounds=[0]
    for m in word:
        bounds.append(bounds[-1]+m+1)
    size=bounds[-1]+1
    rows=[(1<<size)-(2<<i) for i in range(size)]
    for i,j in zip(bounds,bounds[1:]):
        rows[i]^=1<<j
        rows[j]^=1<<i
    return tuple(rows)


def sparse_paths(rows, specials, final):
    """Allocate an increasing sequence of ordinary vertices between ordered specials."""
    ordinary=[v for v in range(len(rows)) if v not in specials]
    assert all(rows[u]>>v&1 for i,u in enumerate(ordinary) for v in ordinary[i+1:])
    total=0
    for order in permutations(specials):
        places=len(order)+1
        powers=[3**j for j in range(places)]
        states=[tuple(i//p%3 for p in powers) for i in range(3**places)]
        dp=[0]*len(states)
        dp[0]=1
        @lru_cache(None)
        def moves(entry,leave):
            dest=[]
            for idx,state in enumerate(states):
                dest.append([idx+powers[j]*((1 if leave[j] else 2)-d)
                             for j,d in enumerate(state) if d or entry[j]])
            return dest
        for v in ordinary:
            entry=tuple(j==0 or rows[order[j-1]]>>v&1 for j in range(places))
            leave=tuple((final>>v&1) if j==len(order) else rows[v]>>order[j]&1
                        for j in range(places))
            nxt=[0]*len(states)
            for i,targets in enumerate(moves(entry,leave)):
                if dp[i]:
                    for k in targets: nxt[k]+=dp[i]
            dp=nxt
        for i,state in enumerate(states):
            if 2 in state: continue
            if not state[-1] and not final>>order[-1]&1: continue
            if all(state[j] or rows[order[j-1]]>>order[j]&1 for j in range(1,len(order))):
                total+=dp[i]
    return total


def readings(s,t,contact=None):
    rows=rows_for((s,t))
    if contact is not None:
        rows=flipped(rows,0,contact)
    special=(0,s+1,s+t+2)
    full=(1<<len(rows))-1
    # k uses the physical boundary; t replaces it by a two-way separator.
    k=sparse_paths(rows,special,full)
    nontags=(full^1)&~rows[0]
    separator=(full^1,)+tuple(r|1 for r in rows[1:])
    return k,sparse_paths(separator,special,nontags)


def main():
    for s,t in ((1,2),(2,1),(2,3)):
        for i in (None,1,s+1,s+t+2):
            rows=rows_for((s,t))
            if i is not None: rows=flipped(rows,0,i)
            k,two=readings(s,t,i)
            full=(1<<len(rows))-1
            assert k==subpaths(rows)[-1]
            sep=(full^1,)+tuple(r|1 for r in rows[1:])
            # The separator has arcs both ways. Reverse by transposing its adjacency.
            rev=tuple(sum(((sep[j]>>i)&1)<<j for j in range(len(sep))) for i in range(len(sep)))
            assert two==subpaths(rev,(full^1)&~rows[0])[-1]
    F=prod(height(word) for word in WORDS)
    center=(1<<60)*((1<<61)-1)
    assert F*32769*((1<<61)+1)==center-1
    prefix=sum(WORDS,())
    whole=prefix+(15,61)
    assert height(whole)==center-1 and len(rows_for(whole))==133 and strong(rows_for(whole))
    options={}
    for i,word in enumerate(WORDS):
        for w in (word,word[::-1]):
            G=(F//height(w))*height(w[:-1])
            options[G]=sum((v for j,v in enumerate(WORDS) if j!=i),())+w
    best=(0,None)
    at_three=(0,None)
    nearest=None
    admissible=[]
    for s,t in ((15,61),(61,15)):
        k,two=readings(s,t)
        for i in range(1,s+t+3):
            k1,two1=readings(s,t,i)
            dk,dt=k1-k,two1-two
            n=2-F*dt
            if dk==dt:
                if n==0: admissible.append((s,t,i,"every G"))
            elif n%(dk-dt)==0 and 0<n//(dk-dt)<=F:
                admissible.append((s,t,i,n//(dk-dt)))
            for G,word in options.items():
                change=G*dk+(F-G)*dt
                score=prod(p for p in PRIMES if (change-2)%p==0)
                if i!=s+1 and change>0 and (nearest is None or change<nearest[0]):
                    nearest=(change,(s,t,i,G))
                if i!=s+1 and score%3==0 and score>at_three[0]:
                    at_three=(score,(s,t,i,G,word,change))
                if i!=s+1 and score>best[0]:
                    host=rows_for(word+(s,t))
                    start=len(rows_for(word))-1
                    assert strong(flipped(host,start,start+i))
                    best=(score,(s,t,i,G,word,change))
    assert not admissible
    for score,(s,t,i,G,word,change) in (best,at_three):
        graph=rows_for(word+(s,t))
        start=len(rows_for(word))-1
        assert strong(graph) and strong(flipped(graph,start,start+i))
        assert (change-2)%score==0
    print("integer solutions for 0<G<=F",admissible)
    print("smallest upward jump in this prefix bank",nearest)
    print("strong contact best simultaneous modulus",best,flush=True)
    print("strong contact meeting upper mod three",at_three,flush=True)


def all_tail_contacts():
    target=(1<<60)*((1<<61)-1)-1
    widths=[m for m in range(1,target.bit_length()+1) if target%((1<<m)+1)==0]
    pairs=contacts=0
    eligible=[]
    for s in widths:
        for t in widths:
            F,rem=divmod(target,((1<<s)+1)*((1<<t)+1))
            if rem: continue
            pairs+=1
            k,two=readings(s,t)
            for i in range(1,s+t+3):
                contacts+=1
                k1,two1=readings(s,t,i)
                dt=two1-two
                numerator=2-F*dt
                denominator=k1-k-dt
                if denominator==0:
                    if numerator==0: eligible.append((s,t,i,"all"))
                elif numerator%denominator==0 and 0<numerator//denominator<=F:
                    eligible.append((s,t,i,numerator//denominator))
    print("eligible widths",widths,"ordered pairs",pairs,"contacts",contacts,"registers",eligible)
    assert not eligible


if __name__=="__main__":
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--all-widths",action="store_true")
    args=parser.parse_args()
    if args.all_widths:
        all_tail_contacts()
    else:
        main()
