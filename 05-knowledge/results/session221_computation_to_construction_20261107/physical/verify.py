#!/usr/bin/env python3
"""Exact tape-to-tournament exclusion calculus. No third-party imports."""
import itertools, json, random
from pathlib import Path

def require(p, msg):
    if not p: raise RuntimeError(msg)

def weight(d):
    require(d >= 1, 'positive distance')
    return 1 if d == 1 else 1 << (d-2)

def energy(w):
    return sum((1-w[i])*w[j]*weight(j-i) for i in range(len(w)) for j in range(i+1,len(w)))

def H(w): return 1+2*energy(w)

def graph(w):
    # Vertex 0 is apex; chain vertices are 1..len(w).
    out=[0]*(len(w)+1)
    for i,b in enumerate(w,1):
        if b: out[i] |= 1
        else: out[0] |= 1<<i
        for j in range(i+1,len(out)): out[i] |= 1<<j
    return out

def flip(out,u,v):
    out=list(out); out[u]^=1<<v; out[v]^=1<<u
    return out

def relabel(out,perm):
    ans=[0]*len(out)
    for v,row in enumerate(out):
        for u in range(len(out)):
            if row>>u&1: ans[perm[v]] |= 1<<perm[u]
    return ans

def single_flip_check(a,b):
    sites=[j for j in range(len(a)) if a[j]!=b[j]]
    require(len(sites)==2 and sites[1]==sites[0]+1,'two adjacent changed tape cells')
    u,v=sites[0]+1,sites[1]+1
    changed=flip(graph(a),u,v)
    perm=list(range(len(changed))); perm[u],perm[v]=v,u
    require(relabel(changed,perm)==graph(b),'single physical arc flip and endpoint relabeling')
    require(sum((x^y).bit_count() for x,y in zip(graph(a),changed))==2,'exactly one reversed arc')
    return changed,[u,v]

def subset_h(out):
    n=len(out); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,c in enumerate(dp[mask]):
            if c:
                avail=out[i]&~mask
                while avail:
                    bit=avail&-avail; avail-=bit
                    dp[mask|bit][bit.bit_length()-1]+=c
    return sum(dp[-1])

def hp_paths(out):
    n=len(out); result=set()
    def go(path,mask):
        if len(path)==n:
            result.add(tuple(path)); return
        avail=out[path[-1]]&~mask
        while avail:
            bit=avail&-avail; avail-=bit
            go(path+[bit.bit_length()-1],mask|bit)
    for v in range(n): go([v],1<<v)
    return result

def injection_check(w,i):
    old=graph(w); u,v=i+1,i+2; new=flip(old,u,v)
    before=hp_paths(old); after=hp_paths(new); images=set()
    for path in before:
        q=list(path)
        for k in range(len(q)-1):
            if q[k:k+2]==[u,v]:
                q[k:k+2]=[v,u]; break
        images.add(tuple(q))
    require(len(images)==len(before) and images<=after,'literal Hamilton-path injection')
    surplus=after-images
    expected={p for p in after if any(p[k:k+3] in ((0,v,u),(v,u,0)) for k in range(len(p)-2))}
    require(surplus==expected,'exact unmatched surplus path language')
    require(len(surplus)==2*response(w,i),'surplus count')

def partition_h(w):
    n=len(w); ans=0
    for mask in range(1<<n):
        left=[i for i in range(n) if mask>>i&1]
        right=[i for i in range(n) if not(mask>>i&1)]
        ans += (not left or w[left[-1]]==1) and (not right or w[right[0]]==0)
    return ans

def response(w,i):
    require(tuple(w[i:i+2])==(1,0),'legal 10 switch')
    L=sum((1-w[j])*(1<<(i-j-2)) for j in range(i-1))
    R=sum(w[j]*(1<<(j-i-3)) for j in range(i+3,len(w)))
    return 1+L+R

def switched(w,sites):
    z=list(w)
    for i in sites:
        require(tuple(w[i:i+2])==(1,0),'parallel old-state eligibility')
        z[i:i+2]=[0,1]
    return tuple(z)

def rule184(w):
    # Closed boundaries: no particle enters or leaves.
    z=(0,)+tuple(w)+(1,)
    return tuple((z[i-1] and not z[i]) or (z[i] and z[i+1]) for i in range(1,len(z)-1))

def f(n):
    total=0
    while n:
        u,b=divmod(n,2)
        if b: total += u+(u&1)
        n=u
    return total

def diff_pair(N):
    require(N>=1,'positive difference parameter')
    m=N-1
    digits=tuple(map(int,bin(m)[2:])) if m else ()
    prefix=tuple(1-b for b in digits)+(1,)
    return prefix+(1,0),prefix+(0,1)

def strong_pair(N):
    if N == 1: return (0,1,0,1),(0,0,1,1)
    a,b=diff_pair(N)
    return a+(1,),b+(1,)

def is_strong(out):
    def reach(adj):
        seen=1; fringe=1
        while fringe:
            bit=fringe&-fringe; fringe-=bit
            fresh=adj[bit.bit_length()-1]&~seen
            seen |= fresh; fringe |= fresh
        return seen == (1<<len(adj))-1
    rev=[sum(1<<j for j,row in enumerate(out) if row>>i&1) for i in range(len(out))]
    return reach(out) and reach(rev)

def words(n): return itertools.product((0,1),repeat=n)

def main():
    counts={}
    literal=0
    for n in range(10):
        for w in words(n):
            require(H(w)==subset_h(graph(w)),f'literal {w}')
            require(H(w)==partition_h(w),f'partition {w}')
            literal+=1
    counts['literal_all_tapes_length_0_to_9']=literal
    local=parallel=subsets=injections=0; histograms={}; strong_histograms={}
    hostile=None; by_old={}
    for n in range(13):
        for w in words(n):
            sites=[i for i in range(n-1) if w[i:i+2]==(1,0)]
            for i in sites:
                z=switched(w,[i]); delta=H(z)-H(w)
                require(delta==2*response(w,i),f'local {w} {i}')
                if n<=7:
                    injection_check(w,i); injections+=1
                histogram=histograms.setdefault((n,i),{})
                histogram[delta//2]=histogram.get(delta//2,0)+1
                if w[0]==0 and w[-1]==1:
                    hist=strong_histograms.setdefault((n,i),set()); hist.add(delta//2)
                key=(n,H(w),i)
                old=by_old.get(key)
                if old and old[1]!=delta and hostile is None:
                    hostile={'tape_a':old[0],'tape_b':w,'site_zero_based':i,'common_H':H(w),'increments':[old[1],delta]}
                by_old[key]=(w,delta); local+=1
            z=switched(w,sites)
            require(z==rule184(w),f'CA rule {w}')
            predicted=2*(sum(response(w,i) for i in sites)+sum(weight(j-i-1) for i,j in itertools.combinations(sites,2)))
            require(H(z)-H(w)==predicted,f'parallel {w}')
            require((H(z)>H(w))==bool(sites),f'Lyapunov {w}')
            parallel+=1
            if n<=9:
                for mask in range(1<<len(sites)):
                    chosen=[i for k,i in enumerate(sites) if mask>>k&1]
                    increment=2*(sum(response(w,i) for i in chosen)+sum(weight(j-i-1) for i,j in itertools.combinations(chosen,2)))
                    require(H(switched(w,chosen))-H(w)==increment,f'subset {w} {chosen}')
                    subsets+=1
    for (n,i),hist in histograms.items():
        a=max(i-1,0); b=max(n-i-3,0); u=n-2-a-b
        require(set(hist)==set(range(1,(1<<a)+(1<<b))),'fixed-position exact interval')
        for N,mult in hist.items():
            convolution=max(0,min((1<<a)-1,N-1)-max(0,N-(1<<b))+1)
            require(mult==(1<<u)*convolution,'fixed-position exact multiplicity')
    for (n,i),values in strong_histograms.items():
        a=i-1; b=n-i-3
        low=1+((1<<(a-1)) if a else 0)+((1<<(b-1)) if b else 0)
        high=(1<<a)+(1<<b)-1
        require(values==set(range(low,high+1)),'strong fixed-position exact interval')
    counts['literal_path_injections_length_2_to_7']=injections
    counts['strong_fixed_position_intervals']=len(strong_histograms)
    counts['exact_fixed_position_histograms']=len(histograms)
    counts.update(local_switches=local,parallel_all_tapes_length_0_to_12=parallel,arbitrary_switch_subsets_length_0_to_9=subsets)
    trajectories=0
    for n in range(11):
        for w in words(n):
            z=w; t=0
            while any(z[i:i+2]==(1,0) for i in range(n-1)):
                z=tuple(map(int,rule184(z))); t+=1
                require(t<=n*n//4,'termination bound')
            r=sum(w); zeros=n-r
            require(z==(0,)*zeros+(1,)*r,'sorted endpoint')
            expected=1 if not r or not zeros else (1<<n)-(1<<zeros)-(1<<r)+3
            require(H(z)==expected,'endpoint energy')
            trajectories+=1
    counts['full_trajectories_length_0_to_10']=trajectories
    examples=[]
    for N in range(1,4097):
        a,b=diff_pair(N); m=N-1
        A=4*N-3-f(m); B=6*N-3-f(m)
        require(H(a)==A and H(b)==B and B-A==2*N,'difference family')
        require(len(graph(a))==m.bit_length()+4,'vertex cost')
        require(sum((x^y).bit_count() for x,y in zip(graph(a),graph(b)))==4,'exactly two reversed arcs')
        flipped,arc=single_flip_check(a,b)
        if N<=64:
            require(subset_h(graph(a))==A and subset_h(flipped)==B,'literal pair')
        if N in (1,2,3,7,21,47,63,1000,4096):
            examples.append({'N':N,'tape_a':a,'tape_b':b,'H_a':A,'H_b':B,'vertices':len(graph(a)),'out_a':graph(a),'out_b_tape_gauge':graph(b),'out_b_single_flip':flipped,'flipped_chain_arc':arc})
    strong_examples=[]
    for N in range(1,4097):
        a,b=strong_pair(N); m=N-1
        A,B=(9,11) if N==1 else (12*N-9-f(m),14*N-9-f(m))
        require(H(a)==A and H(b)==B and B-A==2*N,'strong difference formula')
        flipped,arc=single_flip_check(a,b)
        require(is_strong(graph(a)) and is_strong(flipped),'strongness')
        require(len(a)+1==5+m.bit_length(),'optimal strong order')
        if N in (1,2,7,21,47,63):
            strong_examples.append({'N':N,'tape_a':a,'tape_b':b,'H_a':A,'H_b':B,'vertices':len(a)+1,'out_a':graph(a),'out_b_single_flip':flipped,'flipped_chain_arc':arc})
        if N<=64:
            require(subset_h(graph(a))==A and subset_h(flipped)==B,'strong literal')
    counts['strong_difference_pairs_formula']=4096; counts['strong_difference_pairs_literal']=64
    counts['difference_pairs_formula']=4096; counts['difference_pairs_literal']=64
    # Large arbitrary target: the theorem, not an extrapolation, supports this test.
    huge=(1<<1000)+(1<<601)+137
    a,b=diff_pair(huge)
    require(H(b)-H(a)==2*huge,'1001-bit difference')
    continuation=0
    for n in range(11):
        for w in words(n):
            comp=sum((1-b)<<(n-1-i) for i,b in enumerate(w))
            for gap in (1,2,5):
                require(H(w+(0,)*gap+(1,))-H(w)==(1<<gap)*(comp+1),'continuation observer')
                continuation+=1
    counts['suffix_observer_checks']=continuation
    # Hostile: isolated switches do not add independently in a parallel CA step.
    w=(1,0,1,0); sites=(0,2)
    require(H(switched(w,sites))-H(w)==sum(2*response(w,i) for i in sites)+2,'interaction hostile')
    result={'status':'PASS','counts':counts,'difference_examples':examples,'strong_single_flip_examples':strong_examples,'same_H_local_response_hostile':hostile,'parallel_interaction_hostile':{'tape':w,'old_H':H(w),'new_H':H(switched(w,sites)),'sum_individual_increments':sum(2*response(w,i) for i in sites),'interaction':2},'large_target':{'N_bit_length':huge.bit_length(),'vertices':len(a)+1,'exact_difference':True},'scope':'Exact physical difference constructor; no absolute target coverage or scalar radix gate.'}
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
