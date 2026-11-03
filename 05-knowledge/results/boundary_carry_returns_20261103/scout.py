#!/usr/bin/env python3
"""Exhaust a declared two-complement/common-frame/one-birth language.

Transitive substitutions use fixed clone addresses. Finite polynomial
identities are promoted only with the inherited ordered-run bound.
"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import argparse, json


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def reverse(t):
    f=(1<<len(t))-1
    return tuple(f ^ row ^ (1<<i) for i,row in enumerate(t))


def endpoint(t):
    n=len(t);f=(1<<n)-1
    dp=[[0]*n for _ in range(f+1)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,f+1):
        for i,c in enumerate(dp[s]):
            if not c:continue
            bits=t[i]&~s
            while bits:
                b=bits&-bits;bits-=b
                dp[s|b][b.bit_length()-1]+=c
    return dp


@lru_cache(None)
def hp(t):
    return sum(endpoint(t)[-1])


def paths(t):
    def walk(p,left):
        if not left:
            yield p;return
        bits=left if not p else left&t[p[-1]]
        while bits:
            b=bits&-bits;bits-=b
            yield from walk(p+(b.bit_length()-1,),left^b)
    return walk((),(1<<len(t))-1)


def native(t,p):
    require(len(p)==len(t) and len(set(p))==len(t),'path names')
    require(all(t[a]>>b&1 for a,b in zip(p,p[1:])),'legal native path')
    out=list(reverse(t))
    for a,b in zip(p,p[1:]):
        out[a]^=1<<b;out[b]^=1<<a
    return tuple(out)


def birth(t,mask):
    n=len(t)
    return tuple(row|((1<<n) if mask>>i&1 else 0) for i,row in enumerate(t))+(((1<<n)-1)^mask,)


def strong(t):
    for rows in (t,reverse(t)):
        seen=1
        while True:
            nxt=seen
            for i,row in enumerate(rows):
                if seen>>i&1:nxt|=row
            if nxt==seen:break
            seen=nxt
        if seen!=(1<<len(t))-1:return False
    return True


def clone(t,v,m):
    address=[(i,j) for i in range(len(t)) for j in range(m if i==v else 1)]
    out=[]
    for x,(a,ai) in enumerate(address):
        out.append(sum(1<<y for y,(b,bi) in enumerate(address)
                       if (a==b and ai<bi) or (a!=b and t[a]>>b&1)))
    return tuple(out)


def frames(p):
    return {v:(p[i-1] if i else -1,p[i+1] if i+1<len(p) else -1) for i,v in enumerate(p)}


def run_bound(t,v):
    d=t[v].bit_count()
    return min(d,len(t)-1-d)+1


@lru_cache(None)
def all_births(t):
    n=len(t);f=(1<<n)-1
    ends=endpoint(t);starts=endpoint(reverse(t))
    w=[[0]*n for _ in range(n)]
    for s in range(1,f):
        for i,a in enumerate(ends[s]):
            if a:
                for j,b in enumerate(starts[f^s]):
                    w[i][j]+=a*b
    linear=[ends[f][i]-starts[f][i]+sum(w[i]) for i in range(n)]
    counts=[sum(ends[f])]*(f+1)
    for mask in range(1,f+1):
        bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
        counts[mask]=counts[rest]+linear[i]-sum(w[i][j]+w[j][i] for j in range(n) if rest>>j&1)
    return tuple(counts)


def search(t, ports, excluded=(), first_limit=None):
    n=len(t);h=hp(t)
    srcvals={v:[hp(clone(t,v,m)) for m in range(1,6)] for v in ports}
    programmes={};stats={'source_H':h,'first_paths':0,'path_pairs':0,'common_frames':0,'scalar_hits':0,'uniform_programmes':0}
    for rank,p in enumerate(paths(t)):
        if first_limit is not None and rank>=first_limit:break
        stats['first_paths']+=1;u=native(t,p);fp=frames(p)
        for q in paths(u):
            stats['path_pairs']+=1;fq=frames(q)
            eligible=[v for v in ports if fp[v]==fq[v]]
            if not eligible:continue
            stats['common_frames']+=1;w=native(u,q)
            if hp(w)>3*h+4:continue
            for mask,value in enumerate(all_births(w)):
                c=value-3*h
                if c not in (2,4):continue
                stats['scalar_hits']+=1;target=birth(w,mask)
                if not strong(target):continue
                for v in eligible:
                    bound=max(run_bound(t,v),run_bound(target,v))
                    if all(hp(clone(target,v,m))==3*srcvals[v][m-1]+c for m in range(2,bound+1)):
                        key=(v,c,target)
                        if key in programmes:continue
                        require(hp(clone(target,v,bound+1))==3*srcvals[v][bound]+c,'independent next clone size')
                        programmes[key]={'clone':v,'offset':c,'source':t,'target':target,'first_path':p,'second_path':q,
                                         'mask':mask,'run_bound':bound,'source_H_values':srcvals[v],
                                         'target_H_values':[hp(clone(target,v,m)) for m in range(1,bound+2)],
                                         'frame':fp[v]}
    result=list(programmes.values());stats['uniform_programmes']=len(result)
    stats['clone_offset_counts']={str(v):{str(c):sum(r['clone']==v and r['offset']==c for r in result) for c in (2,4)}for v in ports}
    return {'statistics':stats,'programmes':result}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input');ap.add_argument('--output',default='scout31.json');ap.add_argument('--limit',type=int)
    args=ap.parse_args()
    data=json.loads(Path(args.input).read_text()) if args.input else {'source':[10,36,1,22,7,29],'ports':[0,2,3,4]}
    result=search(tuple(data['source']),data['ports'],first_limit=args.limit)
    Path(__file__).with_name(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['statistics'],indent=2),flush=True)


if __name__=='__main__':main()
