"""Independent rooted-five census and source-recursive transitive counts.

This uses only the two published connection sets, imports no producer code,
and enumerates only five-subsets containing a fixed root. TT6 is counted by
unique-source recursion rather than the producer's six-vertex score test.
"""
from itertools import combinations,permutations
from functools import lru_cache
from pathlib import Path
from hashlib import sha256
from math import comb,isqrt
import argparse,json

CHECKS=0
def need(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(message)

SETS=((1,2,3,4,6,10,12,14,15,16,18),(1,2,3,4,7,9,10,12,15,17,18))

def rooted_dictionary():
    pairs=list(combinations(range(5),2));lookup={e:i for i,e in enumerate(pairs)}
    transformations=[]
    for perm in permutations(range(1,5)):
        p=(0,)+perm;indices=[];flip=0
        for k,(i,j)in enumerate(pairs):
            a,b=p[i],p[j]
            indices.append(lookup[tuple(sorted((a,b)))])
            if a>b:flip|=1<<k
        transformations.append((indices,flip))
    canonical=[]
    for code in range(1024):
        canonical.append(min(sum(((code>>old)&1)<<k for k,old in enumerate(ix))^flip for ix,flip in transformations))
    need(len(set(canonical))==48,'48 distinguished-root tournament types on five vertices')
    need(all(canonical[canonical[c]]==canonical[c]for c in range(1024)),'root canonicalizer idempotent')
    return canonical

def analyze(S,canonical):
    n=23;S=set(S)
    need(len(S)==11 and all((a in S)!=(23-a in S)for a in range(1,23)),'actual cyclic tournament')
    out=[sum(1<<j for j in range(n)if i!=j and (j-i)%n in S)for i in range(n)]
    need(all(x.bit_count()==11 for x in out),'all vertices have eleven outneighbors')
    counts={c:0 for c in sorted(set(canonical))};pairs=list(combinations(range(5),2))
    for other in combinations(range(1,n),4):
        vertices=(0,)+other
        code=sum(((out[vertices[i]]>>vertices[j])&1)<<k for k,(i,j)in enumerate(pairs))
        counts[canonical[code]]+=1
    need(sum(counts.values())==7315,'entire rooted-five subset universe')
    @lru_cache(None)
    def transitive(mask,k):
        if k==0:return 1
        if mask.bit_count()<k:return 0
        if k==1:return mask.bit_count()
        total=0;todo=mask
        while todo:
            bit=todo&-todo;todo-=bit;v=bit.bit_length()-1
            total+=transitive(mask&out[v],k-1)
        return total
    tt=[transitive((1<<n)-1,k)for k in range(7)]
    need(tt[:4]==[1,23,253,1265],'small unique-source recursion controls; TT3=sum_v binomial(outdegree(v),2)')
    return {'connection_set':sorted(S),'rooted_five_counts':counts,'transitive_counts_0_through_6':tt,
            'rooted_five_records_sha256':sha256(json.dumps(counts,sort_keys=True).encode()).hexdigest()}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    canonical=rooted_dictionary();a,b=[analyze(S,canonical)for S in SETS]
    need(a['rooted_five_counts']==b['rooted_five_counts'],'all 48 rooted-five incidence records agree exactly')
    need(a['transitive_counts_0_through_6'][:6]==b['transitive_counts_0_through_6'][:6],'all smaller transitive counts agree')
    need([a['transitive_counts_0_through_6'][6],b['transitive_counts_0_through_6'][6]]==[1173,1265],'independent exact TT6 separation')
    need(all(1063%d for d in range(2,isqrt(1063)+1)),'prime order for explicit counting bound')
    need(2**531>1062*(comb(1063,5)+1)**12,'exact sufficient five-census pigeonhole inequality')
    result={'status':'FINITE-EXACT independent fixed-root census and unique-source TT recursion',
            'cases':[a,b],'checks':CHECKS}
    raw=(json.dumps(result,indent=2,sort_keys=True)+'\n').encode();p=Path(__file__).with_suffix('.json')
    if args.write:p.write_bytes(raw)
    else:need(p.read_bytes().replace(b'\r\n',b'\n')==raw,'frozen independent audit receipt')
    print(json.dumps({'checks':CHECKS,'rooted_types':len(a['rooted_five_counts']),'rooted_records_agree':True,'TT6':[1173,1265]}))

if __name__=='__main__':main()
