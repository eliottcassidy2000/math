#!/usr/bin/env python3
"""Independent Hamilton-language/core reconstruction audit. No producer imports.
All labels 0..n-1. Dense masks: lex unordered pairs, bit1 i->j.
X_n masks: lex nonconsecutive pairs, bit1 i->j, i->i+1 fixed.
"""
import argparse,itertools,json,random,time,pathlib

def need(condition, context=None):
    if not condition: raise AssertionError(context)

def tournament(n,bits,fixed=False):
    a=[0]*n;j=0
    for u in range(n):
        for v in range(u+1,n):
            f=1 if fixed and v==u+1 else (bits>>j)&1
            if not(fixed and v==u+1):j+=1
            if f:a[u]|=1<<v
            else:a[v]|=1<<u
    return a

def core_dfs(a):
    n=len(a);full=(1<<n)-1;out=[0]*n;count=0;p=[]
    def walk(mask,last):
        nonlocal count
        if mask==full:
            count+=1
            for u,v in zip(p,p[1:]):out[u]|=1<<v
            return
        rest=(full^mask) if last<0 else a[last]&(full^mask)
        while rest:
            bit=rest&-rest;rest-=bit;v=bit.bit_length()-1;p.append(v);walk(mask|bit,v);p.pop()
    walk(0,-1)
    return out,count

def reach(a,remove=-1):
    n=len(a);allow=((1<<n)-1)^((1<<remove) if remove>=0 else 0)
    r=[a[u]&allow if u!=remove else 0 for u in range(n)]
    for k in range(n):
        for u in range(n):
            if r[u]>>k&1:r[u]|=r[k]
    return r

def components(a):
    n=len(a);r=reach(a);todo=set(range(n));out=[]
    while todo:
        u=min(todo);c={v for v in todo if v==u or (r[u]>>v&1 and r[v]>>u&1)}
        todo-=c;out.append(sorted(c))
    # order by how many other components are reachable
    keys={c[0]:-sum(bool(r[c[0]]>>d[0]&1) for d in out if c!=d) for c in out}
    out.sort(key=lambda c:keys[c[0]])
    return out

def decode(U):
    n=len(U);T=U.copy();cs=components(U)
    for i,C in enumerate(cs):
        for D in cs[i+1:]:
            for x in C:
                for y in D:T[x]|=1<<y
        if len(C)==1:continue
        ix={v:i for i,v in enumerate(C)}
        core=[sum(((U[x]>>y)&1)<<ix[y] for y in C) for x in C]
        rr=[reach(core,z) for z in range(len(C))]
        for x,y in itertools.combinations(C,2):
            if (U[x]>>y&1) or (U[y]>>x&1):continue
            directions=set()
            for z,r in enumerate(rr):
                if z in (ix[x],ix[y]):continue
                xy=(r[ix[x]]>>ix[y])&1;yx=(r[ix[y]]>>ix[x])&1
                if xy!=yx:directions.add((x,y) if xy else (y,x))
            need(len(directions)==1,(U,x,y,directions))
            u,v=directions.pop();T[u]|=1<<v
    return T

def check(a):
    U,H=core_dfs(a);need(H>0 and H%2==1)
    need(decode(U)==a,(a,U,decode(U)))
    cs=components(a);need(cs==components(U))
    for C in cs:
        if len(C)<3:continue
        ix={v:i for i,v in enumerate(C)}
        b=[sum(((a[x]>>y)&1)<<ix[y] for y in C) for x in C]
        c=[sum(((U[x]>>y)&1)<<ix[y] for y in C) for x in C]
        need(all(reach(b,z)==reach(c,z) for z in range(len(C))))
    return H,sum(x.bit_count() for x in U)

def relocation_controls():
    def hp(a,p):return all(a[u]>>v&1 for u,v in zip(p,p[1:]))
    def signature(a):
        n=len(a);out=[]
        for i in range(n):
            for j in range(i+1,n):
                for left in (False,True):
                    p=list(range(n))
                    if left:p.insert(i,p.pop(j))
                    else:p.insert(j,p.pop(i))
                    out.append(hp(a,p))
        return tuple(out)
    for n in range(1,4):
        M=(n-1)*(n-2)//2
        need(len({signature(tournament(n,b,True)) for b in range(1<<M)})==1<<M)
    A=tournament(4,0,True);B=tournament(4,5,True)
    need(signature(A)==signature(B))
    p=(1,2,3,0);q=(0,2,3,1);r=tuple(p[t] for t in q)
    need(hp(A,p) and hp(B,p));need(not hp(A,r) and hp(B,r))
    return dict(minimal_collision_order=4,masks=[0,5],rows=[A,B],named_trace=[p,q],composed_permutation=r,first_step_legal=[True,True],second_step_legal=[False,True])

def negative_controls():
    def rejected(f):
        try:f()
        except AssertionError:return True
        return False
    tt=tournament(3,7);core,_=core_dfs(tt)
    n1=rejected(lambda:need(reach(tt,1)==reach(core,1)))
    T=tournament(5,10);U,_=core_dfs(T)
    n2=rejected(lambda:need(U==T))
    # Delete a genuinely traceable arc, rather than a nontraceable absent arc.
    bad=U.copy();need(bad[4]>>3&1);bad[4]^=1<<3
    n3=rejected(lambda:need(decode(bad)==T))
    need(n1 and n2 and n3)
    return dict(nonstrong_TT3_reachability_claim_rejected=n1,naive_core_equals_tournament_rejected=n2,corrupted_core_omitting_4_to_3_rejected=n3)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--max-fixed',type=int,default=7);ap.add_argument('--samples',type=int,default=500);ap.add_argument('--receipt');args=ap.parse_args()
    results=[];rng=random.Random(70473);start=time.time();local=relocation_controls();negative=negative_controls();print(json.dumps(dict(local_guard_controls=local,negative_controls=negative)),flush=True)
    for n in range(1,args.max_fixed+1):
        M=(n-1)*(n-2)//2;total=1<<M;mh=mu=0
        for bits in range(total):
            H,used=check(tournament(n,bits,True));mh=max(mh,H);mu=max(mu,n*(n-1)//2-used)
        row=dict(universe='all X_n fixed-identity-path',n=n,states=total,max_H=mh,max_unused_arcs=mu,decoder_failures=0,reachability_failures=0)
        results.append(row);print(json.dumps(row),flush=True)
    for n in (8,9,10):
        for _ in range(args.samples):check(tournament(n,rng.getrandbits(n*(n-1)//2)))
        # Strong maximal-nontraceability upset controls, plus transitive controls.
        check(tournament(n,0));check(tournament(n,(1<<(n*(n-1)//2))-1))
        a=tournament(n,(1<<(n*(n-1)//2))-1)
        for i in range(0,n-2,2):a[i]^=1<<(i+2);a[i+2]^=1<<i
        check(a)
        row=dict(universe='seed70473 random dense masks plus 3 structured controls',n=n,states=args.samples+3,decoder_failures=0)
        results.append(row);print(json.dumps(row),flush=True)
    receipt=dict(schema="session70-independent-hamilton-core-audit-v1",seed=70473,local_guard_controls=local,negative_controls=negative,results=results)
    if args.receipt:pathlib.Path(args.receipt).write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(dict(results=results,seconds=time.time()-start)),flush=True)
if __name__=='__main__':main()
