#!/usr/bin/env python3
"""Literal tournament controls, independent of arithmetic check.py."""
import itertools
import json
import math

SEEDS = [(16,1,3,7,14),(20,1,2,7,14)]


def need(x, msg):
    if not x:
        raise RuntimeError(msg)


def path_counts(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    total=[0]*(1<<n)
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1 and mask != 1<<v:
                old=mask^(1<<v)
                dp[mask][v]=sum(dp[old][u] for u in range(n)
                               if old>>u&1 and rows[u]>>v&1)
        total[mask]=sum(dp[mask])
    return total


def profile(rows):
    n=len(rows); h=path_counts(rows)
    dp=[[0]*(n+1) for _ in h]; dp[0][0]=1
    for mask in range(1,1<<n):
        bit=mask&-mask; sub=mask
        while sub:
            if sub&bit:
                for j in range(1,n+1):
                    dp[mask][j] += h[sub]*dp[mask^sub][j-1]
            sub=(sub-1)&mask
    return dp[-1]


def ordered(p):
    return [v*math.factorial(i) for i,v in enumerate(p)]


def ordered_tt(r):
    S=[[0]*(r+1) for _ in range(r+1)]; S[0][0]=1
    for n in range(1,r+1):
        for k in range(1,n+1): S[n][k]=k*S[n-1][k]+S[n-1][k-1]
    return ordered(S[-1])


def cycle_formula(A,B,C):
    def a(i): return A[i] if 0<=i<len(A) else 0
    def b(i): return B[i] if 0<=i<len(B) else 0
    return sum(C[q]*(a(q-1)*b(q-1)+a(q-1)*b(q)+a(q)*b(q-1)
                    +3*a(q)*b(q)+a(q+1)*b(q)+a(q)*b(q+1)+a(q+1)*b(q+1))
               for q in range(1,len(C)))


def wrap(word,r,split=1,query=None):
    blocks=[]
    for i,bit in enumerate(word):
        rows=list(SEEDS[bit])
        if query == i:
            rows=[row | ((1<<5) if j != 0 else 0) for j,row in enumerate(rows)]+[1]
        blocks.append(rows)
    sizes=[len(x) for x in blocks]
    N=sum(sizes); rows=[0]*(N+r); offsets=[]; off=0
    for b in blocks:
        offsets.append(off)
        for v,row in enumerate(b):rows[off+v] |= row<<off
        off += len(b)
    for i in range(len(blocks)):
        for j in range(i+1,len(blocks)):
            for u in range(offsets[i], offsets[i]+sizes[i]):
                for v in range(offsets[j], offsets[j]+sizes[j]):rows[u] |= 1<<v
    cut=sum(sizes[:split])
    for c in range(N,N+r):
        for a in range(cut): rows[c] |= 1<<a
        for b in range(cut,N): rows[b] |= 1<<c
        for d in range(c+1,N+r):rows[c] |= 1<<d
    return rows,N


def strong(rows, removed):
    keep=[v for v in range(len(rows)) if v not in removed]
    for start in keep:
        seen={start}; todo=[start]
        while todo:
            u=todo.pop()
            for v in keep:
                if rows[u]>>v&1 and v not in seen:
                    seen.add(v);todo.append(v)
        if len(seen)!=len(keep):return False
    return True


def extend(rows, mask):
    n=len(rows)
    return [row | (0 if mask>>i&1 else 1<<n) for i,row in enumerate(rows)]+[mask]


def main():
    base=[profile(b) for b in SEEDS]
    need(base[0] == base[1] == [0,9,33,31,10,1], 'base')
    query=[]
    for b in SEEDS:
        q=[row | ((1<<5) if j != 0 else 0) for j,row in enumerate(b)]+[1]
        query.append(profile(q))
    delta=[b-a for a,b in zip(*query)]
    need(delta == [0,2,14,12,2,0,0], ('query delta',delta))
    records=[]; graphs=0; cuts=0
    for r in range(1,5):
        Hs=[]; responses={}
        for word in itertools.product((0,1),repeat=2):
            rows,N=wrap(word,r)
            H=path_counts(rows)[-1];graphs+=1
            expected=cycle_formula(ordered(base[word[0]]),ordered(base[word[1]]),ordered_tt(r))
            need(H==expected, ('H',r,word,H,expected));Hs.append(H)
            for size in range(r+1):
                for rm in itertools.combinations(range(N+r),size):
                    is_cut=(set(rm)==set(range(N,N+r)))
                    need(strong(rows,set(rm)) != is_cut, ('cut',r,word,rm))
                    cuts+=1
            for address in range(2):
                rows2,_=wrap(word,r,query=address)
                h=path_counts(rows2)[-1];graphs+=1
                x,y=[ordered(query[bit] if i==address else base[bit])
                     for i,bit in enumerate(word)]
                need(h==cycle_formula(x,y,ordered_tt(r)), ('read',r,word,address))
                key=(address,word[address])
                need(key not in responses or responses[key]==h, ('other-bit',r,key))
                responses[key]=h
        need(len(set(Hs))==1, ('blind',r))
        for a in range(2):need(responses[a,1]-responses[a,0]>=6*9, ('gap',r,a))
        records.append({'r':r,'order':10+r,'H':Hs[0],
                        'reads':[responses[0,0],responses[0,1]],
                        'read_gap':responses[0,1]-responses[0,0]})
    # Full profiles of the four literal 2-vertex-strong storage graphs.
    full=[profile(wrap(w,2)[0]) for w in itertools.product((0,1),repeat=2)]
    need(all(p==full[0] for p in full), 'strong full profile')
    # Higher-connectivity hostiles to an articulation-based recovery argument.
    larger=[]
    for r in (5,9,10,20):
        split=(r+5)//5; M=2*split
        rows,N=wrap([0]*M,r,split)
        need(strong(rows,set()) and not strong(rows,set(range(N,N+r))), ('large cut',r))
        for v in range(N+r):need(strong(rows,{v}), ('large single',r,v))
        larger.append([r,M,N+r])
    coupled=[]
    for word in [*itertools.product((0,1),repeat=2),*itertools.product((0,1),repeat=3)]:
        rows,_=wrap(word,1)
        mask=sum(1<<(5*i) for i in range(len(word)))
        ext=extend(rows,mask);H=path_counts(ext)[-1]
        need(H % (3**max(len(word)-2,0)) == 0, ('port budget',word,H))
        v=0;q=H
        while q%3==0:v+=1;q//=3
        coupled.append({'word':word,'rows':ext,'H':H,'v3':v})
    need(coupled[0]['H']==30431 and all(30431%d for d in range(2,math.isqrt(30431)+1)),
         'prime global read')
    triple={tuple(x['word']):x['H'] for x in coupled if len(x['word'])==3}
    need(triple[0,1,0]==12573525 and triple[0,1,0]%9==3,'sharp M3 port budget')
    cubic=sum((-1)**(3-sum(w))*H for w,H in triple.items())
    need(cubic==0,('third mixed difference',cubic))
    need(coupled[3]['H']-coupled[1]['H']-coupled[2]['H']+coupled[0]['H']==312,
         'sharp quadratic interaction')
    print(json.dumps({'status':'FINITE-EXACT controls for TM1–TM3 in TOURNAMENT.md',
                      'query_profiles':query,'query_profile_difference':delta,
                      'literal_H_graphs':graphs,'all_small_deletion_sets':cuts,
                      'small_wrappers':records,'r2_common_profile':full[0],
                      'coupled_one_vertex_reads':coupled,
                      'coupled_third_mixed_difference':cubic,
                      'larger_connectivity_controls_r_M_order':larger},indent=2,sort_keys=True))


if __name__=='__main__':main()
