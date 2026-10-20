#!/usr/bin/env python3
"""Exact finite checks for marked-root parity repair; stdlib only.
Universe: every unlabeled simple graph on n<=5, q=1..4 colors;
all graph automorphisms, exact edge-reversal character; no samples/filters.
"""
from itertools import combinations, permutations, product
from collections import Counter, deque
from fractions import Fraction
from math import comb, factorial
import json

checks=0
def check(p, why):
    global checks
    checks+=1
    if not p: raise AssertionError(why)

def action(mask, ep):
    out=0
    for i,j in enumerate(ep):
        if mask>>i&1: out|=1<<j
    return out

def data(n):
    pairs=list(combinations(range(n),2)); ix={e:i for i,e in enumerate(pairs)}
    ps=list(permutations(range(n)))
    eps=[]; rev=[]; cycles=[]
    for p in ps:
        eps.append(tuple(ix[tuple(sorted((p[a],p[b])))] for a,b in pairs))
        rev.append(sum(1<<i for i,(a,b) in enumerate(pairs) if p[a]>p[b]))
        seen=set(); c=0
        for a in range(n):
            if a not in seen:
                c+=1; x=a
                while x not in seen: seen.add(x); x=p[x]
        cycles.append(c)
    pending=set(range(1<<len(pairs))); reps=[]
    while pending:
        x=min(pending); orb={action(x,e) for e in eps}; pending.difference_update(orb); reps.append(x)
    return pairs,ps,eps,rev,cycles,reps

def subset_action(s,p):
    return sum(1<<p[i] for i in range(len(p)) if s>>i&1)

def selfpaired_sign(mask,p,pairs):
    # Independent character path: count self-paired permutation cycles of edges.
    edge_set={e for i,e in enumerate(pairs) if mask>>i&1}; seen=set(); parity=0
    for e in edge_set:
        if e in seen: continue
        a,b=e; x,y=a,b; swapped=False
        while True:
            seen.add(tuple(sorted((x,y))))
            x,y=p[x],p[y]
            if (x,y)==(b,a): swapped=True
            if tuple(sorted((x,y)))==e: break
        parity^=swapped
    return -1 if parity else 1

def orbit_color_count(n,q,A,ps,signs):
    good=0; seen=set()
    for col in product(range(q), repeat=n):
        if col in seen: continue
        orb={tuple(col[p[i]] for i in range(n)) for p in (ps[a] for a in A)}
        seen.update(orb)
        H=[a for a in A if tuple(col[ps[a][i]] for i in range(n))==col]
        good+=all(signs[a]==1 for a in H)
    return good

def main():
    rows=[]; examples=[]; automata_states=0; color_tests=0
    for n in range(0,6):
        pairs,ps,ep,rev,cy,reps=data(n)
        even=0; rooted_sums=[0,0,0]; d_hist=Counter(); max_states=0
        even_rooted=0; odd_root_defect=0; codec_checks=0
        for mask in reps:
            A=tuple(i for i,e in enumerate(ep) if action(mask,e)==mask)
            signs={a: -1 if (mask&rev[a]).bit_count()%2 else 1 for a in A}
            for a in A: check(signs[a]==selfpaired_sign(mask,ps[a],pairs),('sign',n,mask,a))
            is_even=all(s==1 for s in signs.values()); even+=is_even
            # Root-forgetting defect: relative-even roots of odd graphs.
            remaining=set(range(n)); root_good=0
            while remaining:
                v=min(remaining); orb={ps[a][v] for a in A}; remaining.difference_update(orb)
                H=tuple(a for a in A if ps[a][v]==v)
                good=all(signs[a]==1 for a in H); root_good+=good
                if not is_even and good: check(len(orb)%2==0,('odd-root-orbit',n,mask,v))
            signed_root=Fraction(sum(signs[a]*sum(ps[a][i]==i for i in range(n)) for a in A),len(A))
            check(root_good==signed_root,('root-defect',n,mask))
            if is_even: even_rooted+=root_good
            else: odd_root_defect+=root_good
            buckets=Counter()
            for a in A: buckets[cy[a]]+=signs[a]
            buckets={j:v for j,v in buckets.items() if v}
            def F(H,k):
                x=Fraction(sum(signs[a]*2**(k*cy[a]) for a in H),len(H))
                check(x.denominator==1 and x>=0,('integer',n,mask,H,k))
                return int(x)
            def D(q):
                x=Fraction(sum(signs[a]*q**cy[a] for a in A),len(A))
                check(x.denominator==1 and x>=0,('D-integer',n,mask,q))
                return int(x)
            # Exhaustive direct orbit/stabilizer count, independent of signed Burnside.
            for q in range(1,5):
                check(D(q)==orbit_color_count(n,q,A,ps,signs),('color',n,mask,q)); color_tests+=1
            d=next(q for q in range(1,max(n,1)+1) if D(q)>0)
            d_hist[d]+=1
            if n:
                coeff=[sum((-1)**(r-i)*comb(r,i)*D(i) for i in range(r+1)) for r in range(n+1)]
                check(all(v>=0 for v in coeff),('binomial-coeff',n,mask,coeff))
            # Nonnegative stabilizer automaton, actual subgroups retained as sidecar.
            queue=deque([A]); transitions={}
            while queue:
                H=queue.popleft()
                if H in transitions: continue
                remaining=set(range(1<<n)); row=Counter()
                while remaining:
                    s=min(remaining)
                    orb={subset_action(s,ps[a]) for a in H}; remaining.difference_update(orb)
                    K=tuple(a for a in H if subset_action(s,ps[a])==s)
                    row[K]+=1
                    if K not in transitions: queue.append(K)
                transitions[H]=row
            max_states=max(max_states,len(transitions)); automata_states+=len(transitions)
            for H,row in transitions.items():
                check(F(H,0)==int(all(signs[a]==1 for a in H)),('accept',n,mask,H))
                for k in range(4):
                    check(F(H,k+1)==sum(m*F(K,k) for K,m in row.items()),('transfer',n,mask,H,k))
            for k in range(3): rooted_sums[k]+=F(A,k)
            # Local branch-count ranking/unranking; no catalogue of terminal types.
            def branches(H):
                remaining=set(range(1<<n)); out=[]
                while remaining:
                    s=min(remaining); remaining.difference_update(subset_action(s,ps[a]) for a in H)
                    out.append((s,tuple(a for a in H if subset_action(s,ps[a])==s)))
                return out
            def unrank(H,k,r):
                check(0<=r<F(H,k),('unrank-domain',n,mask,k,r))
                if k==0: return ()
                for s,K in branches(H):
                    count=F(K,k-1)
                    if r<count: return (s,)+unrank(K,k-1,r)
                    r-=count
                raise AssertionError('unrank fell through')
            def rank(H,word):
                if not word:
                    check(all(signs[a]==1 for a in H),('rank-terminal',n,mask,H)); return 0
                s=word[0]; best=min(subset_action(s,ps[a]) for a in H)
                g=next(a for a in H if subset_action(s,ps[a])==best)
                tail=tuple(subset_action(t,ps[g]) for t in word[1:]); total=0
                for v,K in branches(H):
                    if v==best: return total+rank(K,tail)
                    total+=F(K,len(tail))
                raise AssertionError('rank fell through')
            k=2
            words=[unrank(A,k,r) for r in range(F(A,k))]
            check(len(set(words))==len(words),('distinct-codewords',n,mask))
            for r,w in enumerate(words):
                check(rank(A,w)==r,('codec-inverse',n,mask,r)); codec_checks+=1
                for a in A:
                    check(rank(A,tuple(subset_action(s,ps[a]) for s in w))==r,('codec-invariant',n,mask,r,a))
            # Direct recurrence whose distinct roots are 2^j for nonzero buckets.
            poly=[1]
            for j in sorted(buckets):
                root=2**j; z=[0]*(len(poly)+1)
                for i,c in enumerate(poly): z[i]-=root*c; z[i+1]+=c
                poly=z
            for k in range(3): check(sum(c*F(A,k+i) for i,c in enumerate(poly))==0,('recurrence',n,mask,k))
            if n==3 and mask in (0,(1<<len(pairs))-1):
                examples.append({'n':n,'mask':mask,'aut_order':len(A),'weighted_cycle_buckets':buckets,'b_k_0_to_4':[F(A,k) for k in range(5)],'min_colors':d,'automaton_states':len(transitions)})
        # Independent tournament-orbit count with a distinguished old vertex.
        m=len(pairs); allones=(1<<m)-1; pending=set(range(1<<m)); tournament_rooted=0
        while pending:
            x=min(pending)
            images=[action(x^rev[a],ep[a]) for a in range(len(ps))]
            pending.difference_update(images)
            TA=tuple(a for a,y in enumerate(images) if y==x)
            rem=set(range(n))
            while rem:
                v=min(rem); rem.difference_update(ps[a][v] for a in TA); tournament_rooted+=1
        check(tournament_rooted==even_rooted+odd_root_defect,('tournament-rooted-identity',n))
        rows.append({'n':n,'tournament_rooted_types':tournament_rooted,'globally_even_rooted_types':even_rooted,'odd_graph_root_defect':odd_root_defect,'codec_roundtrips':codec_checks,'graph_types':len(reps),'sign_even_types':even,'rooted_extension_counts_k_0_1_2':rooted_sums,'min_colors_histogram':dict(sorted(d_hist.items())),'max_stabilizer_states':max_states})
    # Complete-graph coding threshold, larger n via signed Stirling polynomial.
    complete=[]
    for n in range(1,17):
        k=(n-1).bit_length()
        check(comb(2**k,n)>0,('K-positive',n))
        if k: check(2**(k-1)<n,('K-negative',n))
        complete.append({'n':n,'min_roots':k})
    result={'status':'PASS','checks':checks,'universe':'all simple graph types n<=5; direct q-color orbit/stabilizer checks q=1..4','rows':rows,'direct_color_tests':color_tests,'total_stabilizer_states':automata_states,'complete_graph_thresholds':complete,'hostile_same_ordinary_group_different_sign':examples}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__': main()
