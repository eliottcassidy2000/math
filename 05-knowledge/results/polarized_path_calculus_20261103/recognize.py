"""Recognize and recover a unique rooted Hamilton path in O(n^3) time."""
import json
from explore import tiling, induced, opposite, require, tables, funnel


def ordered_components(T):
    """Kosaraju: O(n^2) adjacency inspections, source-to-sink order."""
    seen=set();finish=[]
    def dfs(v):
        seen.add(v)
        for w in range(len(T)):
            if T[v] >> w & 1 and w not in seen:
                dfs(w)
        finish.append(v)
    for v in range(len(T)):
        if v not in seen:
            dfs(v)
    R=opposite(T);seen=set();blocks=[]
    def collect(v,block):
        seen.add(v);block.append(v)
        for w in range(len(T)):
            if R[v] >> w & 1 and w not in seen:
                collect(w,block)
    for v in reversed(finish):
        if v not in seen:
            block=[];collect(v,block);blocks.append(tuple(block))
    return blocks


def unique_start_path(T,r):
    """Return the unique Hamilton path starting at r, or None if count !=1."""
    if len(T)==1:
        return (r,)
    keep=[v for v in range(len(T)) if v!=r]
    U=induced(T,keep);blocks=ordered_components(U)
    if any(len(block)>1 for block in blocks[1:]):
        return None
    core=blocks[0]
    entries=[j for j,v in enumerate(core) if T[r] >> keep[v] & 1]
    if len(entries)!=1:
        return None
    inner=unique_start_path(induced(U,core),entries[0])
    if inner is None:
        return None
    return (r,)+tuple(keep[core[j]] for j in inner)+tuple(keep[block[0]] for block in blocks[1:])


def class_counts(n):
    u=[0]*(max(n,4)+1)
    u[3]=1;u[4]=3
    for k in range(5,n+1):
        u[k]=3*u[k-1]-u[k-2]
    return u


def unrank(n,index):
    """Decode a strong unique-start rooted class; canonical root is zero."""
    u=class_counts(n)
    require(n>=3 and 0<=index<u[n],('rank range',n,index))
    for m in [1]+list(range(3,n)):
        t=n-m-1;c=1 if t==0 else 1 << (t-1)
        core_count=1 if m==1 else u[m]
        mass=c*core_count
        if index>=mass:
            index-=mass;continue
        word_index,core_index=divmod(index,core_count)
        core=(0,) if m==1 else unrank(m,core_index)
        word=() if t==0 else tuple((word_index >> i)&1 for i in range(t-1))+(1,)
        T,r=funnel(core,{0},word)
        return induced(T,[r]+list(range(r)))
    raise RuntimeError('unrank exhausted')


def rank_class(T,r):
    n=len(T);u=class_counts(n)
    require(n>=3 and len(ordered_components(T))==1,'rank strong input')
    keep=[v for v in range(n) if v!=r]
    U=induced(T,keep);blocks=ordered_components(U)
    require(all(len(b)==1 for b in blocks[1:]),'rank singleton tail')
    core=blocks[0];m=len(core);t=n-m-1
    entries=[j for j,v in enumerate(core) if T[r] >> keep[v] & 1]
    require(len(entries)==1,'rank unique entry')
    core_index=0 if m==1 else rank_class(induced(U,core),entries[0])
    word=tuple((T[keep[b[0]]] >> r)&1 for b in blocks[1:])
    require(not word or word[-1]==1,'rank strong tail')
    word_index=sum(bit << i for i,bit in enumerate(word[:-1]))
    offset=0
    for size in [1]+list(range(3,n)):
        if size==m:
            break
        tail=n-size-1
        offset+=(1 if tail==0 else 1 << (tail-1))*(1 if size==1 else u[size])
    return offset+word_index*(1 if m==1 else u[m])+core_index


def main():
    all_counts=[];strong_counts=[];instances=0
    for n in range(1,8):
        a=s=0
        for bits in range(1 << (max(0,(n-1)*(n-2)//2))):
            T=tiling(n,bits)
            path=unique_start_path(T,0)
            if path is not None:
                require(path==tuple(range(n)),('canonical unique path',T,path))
                a+=1
                s+=int(n>=3 and len(ordered_components(T))==1)
            if n<=5:
                ends=tables(opposite(T))[1][-1]
                for r in range(n):
                    p=unique_start_path(T,r)
                    require((p is not None)==(ends[r]==1),('independent decision',T,r))
                    if p is not None:
                        require(p[0]==r and sorted(p)==list(range(n)) and all(T[u] >> v & 1 for u,v in zip(p,p[1:])),('recovered path',T,r,p))
                    instances+=1
            tables.cache_clear()
        all_counts.append(a);strong_counts.append(s)
    require(all_counts==[1,1,2,5,13,34,89],all_counts)
    require(strong_counts==[0,0,1,3,8,21,55],strong_counts)
    ranks=0
    for n in range(3,11):
        for j in range(class_counts(n)[n]):
            T=unrank(n,j)
            require(rank_class(T,0)==j and unique_start_path(T,0)==tuple(range(n)),('rank round trip',n,j))
            # Relabel to make the physical root nonzero.
            R=induced(T,list(reversed(range(n))))
            require(rank_class(R,n-1)==j,('rank label invariance',n,j))
            ranks+=1
    long=[]
    for n in (50,100):
        total=class_counts(n)[n]
        for j in (0,total//2,total-1):
            T=unrank(n,j)
            require(rank_class(T,0)==j and unique_start_path(T,0)==tuple(range(n)),('long rank',n,j))
        long.append({'n':n,'classes':total,'rank_bits':(total-1).bit_length()})
    print(json.dumps({'all_rooted_classes':all_counts,'strong_rooted_classes':strong_counts,'all_root_DP_decision_controls':instances,'canonical_charts':33868,'small_rank_round_trips':ranks,'long_rank_controls':long},indent=2))


if __name__=='__main__':
    main()
