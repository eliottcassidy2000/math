#!/usr/bin/env python3
"""Exact cap-semiring BFS for ribbon H; source transfer, independent graph DP.
All writes are in this packet. No unbounded/count-magnitude efficiency claim.
"""
import argparse,collections,hashlib,importlib.util,itertools,json,pathlib,time
ROOT=pathlib.Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('pinned_ribbon',ROOT/'sources/ribbon.py')
r=importlib.util.module_from_spec(sp);sp.loader.exec_module(r)
def require(p,msg):
    if not p: raise RuntimeError(msg)
def cap_apply(edges,x,C):
    y=[0]*len(x)
    for j,i,m in edges:
        y[j]=min(C,y[j]+m*x[i])
    return tuple(y)
def exact_apply(edges,x):
    y=[0]*len(x)
    for j,i,m in edges:y[j]+=m*x[i]
    return tuple(y)
def independent_graph(word,w):
    n=len(word);a=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if j-i<=w and ((word[j]>>(j-i-1))&1):a[j]|=1<<i
            else:a[i]|=1<<j
    return a
def independent_h(word,w):
    a=independent_graph(word,w);n=len(a)
    if not n:return 1
    # Endpoint-major flat table, no source profile DP imported.
    dp=[0]*((1<<n)*n)
    for i in range(n):dp[(1<<i)*n+i]=1
    for s in range(1,1<<n):
        for v in range(n):
            c=dp[s*n+v]
            if not c:continue
            free=a[v]&~s
            while free:
                bit=free&-free;free-=bit;u=bit.bit_length()-1
                dp[(s|bit)*n+u]+=c
    return sum(dp[((1<<n)-1)*n:])
def word_for(k,parents,letters):
    ans=[]
    while k>=0:
        ans.append(letters[k]);k=parents[k]
    return list(reversed(ans))
def system(w):
    order,ids,mats,final=r.transfer_system(w)
    masks=list(range(1<<w))
    alph=[frozenset(i+1 for i in range(w) if (m>>i)&1) for m in masks]
    edges=[tuple((j,i,n) for (j,i),n in sorted(mats[a].items())) for a in alph]
    require(all(n>=0 for es in edges for _,_,n in es),'negative transfer coefficient')
    return order,edges,final
def solve(w,B,max_states):
    require(w>=1 and B>=1,"span and bound must both be positive")
    require(max_states>=1,"max-states must be positive")
    order,edges,final=system(w);d=len(order);C=B+1
    init=[0]*d;init[0]=1
    x=cap_apply(edges[0],init,C) # a single vertex; every first letter has same effect
    require(all(cap_apply(e,init,C)==x for e in edges),'initial letter dependence')
    states=[x];seen={x:0};parents=[-1];letters=[0];depths=[1];witness={1:0};q=0;transitions=0;pruned=0
    while q<len(states):
        x=states[q]
        for a,e in enumerate(edges):
            y=cap_apply(e,x,C);transitions+=1
            h=min(C,sum(y[j] for j in final))
            if h>B:pruned+=1;continue
            require(h>=1 and h%2==1,'tournament oddness failed')
            if y not in seen:
                require(len(states)<max_states,f'state budget {max_states} exceeded at span={w}, B={B}; no completeness claim')
                k=len(states);seen[y]=k;states.append(y);parents.append(q);letters.append(a);depths.append(depths[q]+1)
                witness.setdefault(h,k)
        q+=1
    words={str(h):word_for(k,parents,letters) for h,k in sorted(witness.items())}
    # Every emitted realiser is independently reconstructed and counted.
    for h,word in words.items():require(independent_h(word,w)==int(h),f'independent H mismatch: {h}, {word}')
    digest=hashlib.sha256()
    for x in states:digest.update((','.join(map(str,x))+'\n').encode())
    result={'span':w,'bound':B,'cap':C,'profile_dimension':d,'complete':True,'states':len(states),'transitions':transitions,'pruned_over_bound':pruned,'maximum_shortest_state_length':max(depths),'maximum_output_witness_length':max(map(len,words.values())),'state_stream_sha256':digest.hexdigest(),'missing_odd':[h for h in range(1,B+1,2) if h not in witness],'witnesses':words,'scope':'all word lengths, nonempty tournaments only'}
    return result

def cap_row_apply(edges,row,C):
    out=[0]*len(row)
    for j,i,m in edges:out[i]=min(C,out[i]+row[j]*m)
    return tuple(out)
def marked_replacement(word,pos,new_letter,edges,final,C):
    d=1+max(max(i,j) for es in edges for j,i,_ in es)
    left=[0]*d;left[0]=1
    for a in word[:pos]:left=cap_apply(edges[a],left,C)
    right=[int(j in final) for j in range(d)]
    for a in reversed(word[pos+1:]):right=cap_row_apply(edges[a],right,C)
    middle=cap_apply(edges[new_letter],left,C)
    return min(C,sum(a*b for a,b in zip(right,middle)))

def controls():
    # Brute-force graph DP against exact/capped source transfer. Short labels with
    # non-existing offsets are allowed and harmless by the declared grammar.
    count=0
    for w in (1,2,3):
        order,edges,final=system(w);d=len(order)
        for n in range(1,5):
            for word in itertools.product(range(1<<w),repeat=n-1):
                word=(0,)+word;x=[0]*d;x[0]=1;z=x
                for a in word:x=exact_apply(edges[a],x);z=cap_apply(edges[a],z,6)
                h=sum(x[i] for i in final)
                require(h==independent_h(word,w),'profile/graph mismatch')
                require(tuple(min(6,v) for v in x)==z,'cap coordinate mismatch')
                require(min(6,h)==min(6,sum(z[i] for i in final)),'cap readout mismatch')
                count+=1
    # Output-only lumping fails on the actual one-apex family.
    def f(q):
        if q<2:return 1
        if q&1:return f(q//2)
        return f(q//2)+q//2+(q//2)%2
    require(f(4)==f(6)==5 and (f(8),f(12))==(9,11),'one-apex hostile')
    # Capping cannot be applied blindly to signed reduced transitions.
    require(min(6,10-9)!=min(6,min(6,10)-min(6,9)),'signed hostile')
    # Nor to native decrement microsteps.
    require(min(5,12)==min(5,18) and min(5,12-8)!=min(5,18-8),'countdown hostile')
    edit_checks=0;hostile=None;fibres={}
    for w in (1,2):
        order,edges,final=system(w)
        for n in range(1,5):
            alph=[range(1<<min(w,j)) for j in range(n)]
            for word in itertools.product(*alph):
                h=independent_h(word,w)
                for pos in range(1,n):
                    for bit in range(min(w,pos)):
                        b=word[pos]^(1<<bit);edited=word[:pos]+(b,)+word[pos+1:]
                        hn=independent_h(edited,w)
                        require(marked_replacement(word,pos,b,edges,final,6)==min(6,hn),'marked replacement sidecar failed')
                        edit_checks+=1
                        key=(w,n,h,pos,word[pos],b)
                        if key in fibres and fibres[key][1]!=hn and hostile is None:
                            hostile={'span':w,'words':[list(fibres[key][0]),list(word)],'old_H':h,'position_0based':pos,'old_letter':word[pos],'new_letter':b,'new_H':[fibres[key][1],hn]}
                        fibres.setdefault(key,(word,hn))
    require(hostile is not None,'expected scalar-edit hostile missing')
    return {'graph_cap_checks':count,'marked_edit_checks':edit_checks,'same_H_edit_hostile':hostile,'one_apex_same_H_next_zero':[5,5,9,11],'signed_cap_hostile':{'input':[10,9],'uncapped_difference':1,'difference_of_caps':0},'countdown_hostile':[12,18,4,5]}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--span',type=int,default=2);p.add_argument('--bound',type=int,default=99);p.add_argument('--max-states',type=int,default=2000000);p.add_argument('--out',default='receipt.json');args=p.parse_args()
    require(args.span>=1 and args.bound>=1 and args.max_states>=1,'span, bound, and max-states must be positive')
    res={'method':'positive profile cap-semiring fixed point, all word lengths','controls':controls(),'run':solve(args.span,args.bound,args.max_states)}
    (ROOT/args.out).write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res['run'].items() if k!='witnesses'},indent=2));print(json.dumps(res['controls']))
