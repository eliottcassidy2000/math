"""Exhaustive independent orbit coverage of the small source tournament bank."""
import hashlib,itertools,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def main():
    source=Path(__file__).read_bytes();path=ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives.json'
    raw=path.read_bytes();bank=json.loads(raw)['members'];receipt=[]
    for n in range(1,6):
        pairs=list(itertools.combinations(range(n),2));covered=set()
        for entry in bank[str(n)]:
            code=entry['code'];a=[[0]*n for _ in range(n)]
            for k,(u,v) in enumerate(pairs):a[u][v]=(code>>k)&1;a[v][u]=1-a[u][v]
            images=[sum(a[p[u]][p[v]]<<k for k,(u,v) in enumerate(pairs)) for p in itertools.permutations(range(n))]
            orbit=set(images)
            if covered&orbit:raise RuntimeError('duplicate source class')
            if images.count(code)!=entry['aut']:raise RuntimeError('source automorphism mismatch')
            covered.update(orbit)
        if covered!=set(range(1<<len(pairs))):raise RuntimeError('missing labelled tournament')
        receipt.append({'n':n,'classes':len(bank[str(n)]),'labelled_tournaments':len(covered)})
    if source!=Path(__file__).read_bytes() or raw!=path.read_bytes():raise RuntimeError('source changed')
    print(json.dumps({'status':'FINITE-EXACT PASS','orders':receipt,'source_sha256':sha(source),'bank_sha256':sha(raw)},indent=2,sort_keys=True))
if __name__=='__main__':main()
