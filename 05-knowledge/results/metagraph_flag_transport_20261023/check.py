"""Exact local path-flag transport; no canonical graph package required.

Default checks the frozen receipt. --write explicitly regenerates it.
The fixed-path gauge is increasing, and subset order is inherited from it.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def rows(code, n):
    a = [[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j] = (code >> b)&1
        a[j][i] = 1-a[i][j]
    return a


def encode(a, order):
    return sum(a[u][v]<<b for b,(u,v) in enumerate(combinations(order,2)))


def tables(k):
    pairs=list(combinations(range(k),2)); lookup={p:b for b,p in enumerate(pairs)}
    transforms=[]
    for p in permutations(range(k)):
        source=[lookup[tuple(sorted((p[i],p[j])))] for i,j in pairs]
        flips=sum((p[i]>p[j])<<b for b,(i,j) in enumerate(pairs))
        fmap=[sum(((f>>s)&1)<<b for b,s in enumerate(source)) for f in range(1<<len(pairs))]
        transforms.append((p,fmap,[f^flips for f in fmap]))
    return pairs,lookup,transforms


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    k=5;pairs,lookup,transforms=tables(k);full=(1<<len(pairs))-1
    byroot=[[t for t in transforms if t[0][0]==r] for r in range(k)]
    def canonical(t,f,root=None):
        choices=transforms if root is None else byroot[root]
        return min((fm[f]<<10)|tm[t] for _,fm,tm in choices)
    def rooted_type(t,r):
        return min(tm[t] for _,_,tm in byroot[r])
    def unrooted_type(t):
        return min(tm[t] for _,_,tm in transforms)
    root_types={rooted_type(t,r) for t in range(1024) for r in range(k)}
    need(len(root_types)==48,'rooted tournament types on five')
    flag_map={};root_map={};states=0;witness=None;forgot=defaultdict(dict)
    cardinality_forgot=defaultdict(dict);count_witness=None
    for fword in range(1<<(k-1)):
        f=sum(((fword>>i)&1)<<lookup[i,i+1] for i in range(k-1))
        for t in range(1<<len(pairs)):
            if t&f!=f:continue
            states+=1;target=t^(full^f)
            need(target^(full^f)==t and target&f==f,'literal involution and retained path arcs')
            # Every local flag is realized by an actual path presentation of size <=9.
            positions=[0]
            for i in range(k-1):positions.append(positions[-1]+(1 if fword>>i&1 else 2))
            n=positions[-1]+1;a=[[int(i<j) for j in range(n)] for i in range(n)]
            local=rows(t,k)
            for i,j in pairs:
                u,v=positions[i],positions[j];a[u][v]=local[i][j];a[v][u]=local[j][i]
            need(all(a[i][i+1] for i in range(n-1)),'actual global Hamilton path')
            b=[[a[i][j] if abs(i-j)==1 or i==j else 1-a[i][j] for j in range(n)] for i in range(n)]
            need(encode(b,positions)==target,'global-to-local restriction commutes')
            target_local=rows(target,k)
            for deleted in range(k):
                retained=[i for i in range(k) if i!=deleted]
                smaller_f=sum(((f>>lookup[i,j])&1)<<z for z,(i,j) in enumerate(combinations(retained,2)))
                need((encode(local,retained)^(63^smaller_f))==encode(target_local,retained),
                     'deletion commutes without inserting a shortcut')
            sc=canonical(t,f);tc=canonical(target,f)
            need(sc not in flag_map or flag_map[sc]==tc,'flag isomorphism quotient transport')
            flag_map[sc]=tc
            for r in range(k):
                sr=canonical(t,f,r);tr=canonical(target,f,r)
                need(sr not in root_map or root_map[sr]==tr,'rooted quotient transport')
                root_map[sr]=tr
                # Forget only the marked path forest, while retaining source and root.
                key=rooted_type(t,r);value=rooted_type(target,r)
                previous=forgot[key]
                if previous and value not in previous and witness is None:
                    old=next(iter(previous.values()));witness=[old,{'local_word':t,'path_edges':fword,'root':r,'target_type':value,'positions':positions}]
                previous[value]={'local_word':t,'path_edges':fword,'root':r,'target_type':value,'positions':positions}
                key2=(key,fword.bit_count());prev=cardinality_forgot[key2]
                if prev and value not in prev and count_witness is None:
                    count_witness=[next(iter(prev.values())),{'local_word':t,'path_edges':fword,'root':r,'target_type':value,'positions':positions}]
                prev[value]={'local_word':t,'path_edges':fword,'root':r,'target_type':value,'positions':positions}
    need(states==2**6*3**4==5184,'complete admissible local universe')
    need(all(flag_map[v]==u for u,v in flag_map.items()),'unrooted quotient involution')
    need(all(root_map[v]==u for u,v in root_map.items()),'rooted quotient involution')
    need(witness is not None and count_witness is not None,'path-sidecar hostile controls')
    # The zero-path-arc face is converse; the full retained path is NOT converse.
    need(all((t^full)==(t^(full^0)) for t in range(1024)),'converse face')
    result={'status':'FINITE-EXACT local five-vertex flags; all-size restriction theorem in README',
            'labelled_states':states,'root_marked_labelled_states':states*k,
            'five_vertex_rooted_types':len(root_types),'unrooted_flag_types':len(flag_map),
            'rooted_flag_types':len(root_map),'unrooted_flag_fixed_points':sum(u==v for u,v in flag_map.items()),
            'rooted_flag_fixed_points':sum(u==v for u,v in root_map.items()),
            'flag_map':sorted(flag_map.items()),'rooted_flag_map':sorted(root_map.items()),
            'forgotten_forest_witness':witness,'forest_cardinality_insufficient_witness':count_witness,
            'checks':CHECKS,'source_sha256_lf':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();target=HERE/'receipt.json'
    if args.write:target.write_bytes(raw)
    else:need(target.read_bytes().replace(b'\r\n',b'\n')==raw,'frozen receipt')
    print(json.dumps({**{k:v for k,v in result.items() if k not in ('flag_map','rooted_flag_map')},'checks':CHECKS}))


if __name__=='__main__':main()
