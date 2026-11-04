"""Independent permutation controls for counts and full canonical labellings."""
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import check as producer

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
checks=0


def require(ok, message):
    global checks
    checks+=1
    if not ok: raise RuntimeError(message)


def literal(rows):
    n=len(rows)
    start=[0]*n;end=[0]*n;c=[[0]*n for _ in rows]
    for word in itertools.permutations(range(n)):
        bad=[i for i in range(n-1) if not (rows[word[i]]>>word[i+1])&1]
        if not bad:
            start[word[0]]+=1;end[word[-1]]+=1
            for i in range(n-1): c[word[i]][word[i+1]]+=1
        elif len(bad)==1:
            i=bad[0];c[word[i]][word[i+1]]+=1
    return tuple(start),tuple(end),tuple(tuple(row) for row in c)


def brute_canon(rows,data):
    s,e,c=data;n=len(rows)
    best=None;codes=set();aut=0
    for order in itertools.permutations(range(n)):
        key=tuple(z for u in order for z in (s[u],e[u]))+tuple(c[u][v] for u in order for v in order)
        code=sum(((rows[order[i]]>>order[j])&1)<<bit
                 for bit,(i,j) in enumerate(itertools.combinations(range(n),2)))
        if best is None or key<best: best,codes,aut=key,{code},1
        elif key==best: codes.add(code);aut+=1
    return best,codes,aut


def main():
    bank=ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    raw=bank.read_bytes();data=json.loads(gzip.decompress(raw))['members']
    orders=[];counted=0;large=[]
    for n in range(1,7):
        seen=set()
        for row in data[str(n)]:
            t=producer.rows_from_code(n,row['code'])
            exact=literal(t);counted+=1
            require(exact==producer.sidecar(t),('literal counts',n,row['code']))
            key,codes,aut=brute_canon(t,exact)
            require(len(codes)==1,('extra sidecar automorphism',n))
            require(aut==row['aut'],('bank automorphism',n))
            _,fast_codes,fast_aut,_=producer.canonical(t,exact)
            require(fast_aut==aut and len(fast_codes)==1,('refined automorphism',n))
            require(key not in seen,('sidecar class collision',n));seen.add(key)
        orders.append([n,len(seen)])
    # Two named deterministic controls at each larger order, no sampling claim.
    for n in (7,8,9):
        for index in (len(data[str(n)])//2,len(data[str(n)])-1):
            row=data[str(n)][index]
            t=producer.rows_from_code(n,row['code'])
            require(literal(t)==producer.sidecar(t),('large literal counts',n,index));counted+=1
            large.append([n,index,row['code']])
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==SOURCE,'source changed')
    print(json.dumps({'status':'PASS','source_sha256':SOURCE,'bank_sha256':hashlib.sha256(raw).hexdigest(),
                      'checks':checks,'literal_counted_classes':counted,'brute_canonical_orders':orders,
                      'larger_literal_controls':large},indent=2,sort_keys=True))


if __name__=='__main__':main()
