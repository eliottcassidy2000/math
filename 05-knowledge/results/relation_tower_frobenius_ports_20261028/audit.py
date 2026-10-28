"""Alternate-path checks and hostile mutations, independent of search decisions."""
import gzip
import hashlib
import importlib.util
import json
from itertools import product, permutations
from pathlib import Path
from math import prod, gcd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checks=0
def need(b,msg):
    global checks
    checks+=1
    if not b: raise RuntimeError(msg)

def main():
    bank=json.loads((HERE/'prime_bank.json').read_text())
    cert={int(k):v for k,v in bank['prime_certificates'].items()}
    verified={2}
    for p in sorted(set(cert)-{2}):
        fs={int(k):v for k,v in cert[p]['factors'].items()}; a=cert[p]['witness']
        need(all(q in verified for q in fs) and prod(q**v for q,v in fs.items())==p-1,'certificate factors')
        need(pow(a,p-1,p)==1 and all(gcd(pow(a,(p-1)//q,p)-1,p)==1 for q in fs),'certificate witness')
        verified.add(p)
    for row in bank['rows']:
        d=prod(p for i,p in enumerate(bank['coordinate_primes']) if row['mask']>>i&1)
        need(row['prime'] in verified and row['prime']+1==d*row['a']**2,'prime-bank column')
    need([r['prime'] for r in bank['rows'][:15]]==[7,2,5,19,89,59,29,251,13,83,41,139,5669,419,839],'inherited bank prefix')
    # Direct integer products of three successor values, for every PG(10,2)
    # line. The square root is supplied exactly by the even exponents.
    from math import isqrt
    ps=[0]+[r['prime']+1 for r in bank['rows']]
    for a in range(1,2048):
        for b in range(a+1,2048):
            c=a^b
            if b<c:
                total=ps[a]*ps[b]*ps[c]
                need(isqrt(total)**2==total,'integer square line')
    need(isqrt(ps[1]*ps[2]*ps[4])**2!=ps[1]*ps[2]*ps[4],'independent-column hostile')
    result=json.loads((HERE/'root_switch.json').read_text())
    native=json.loads((ROOT/'05-knowledge/results/session67_metagraph_enumeration_information_20261026/session67-root/output/census.json').read_text())
    for x,y in zip(result['results'],native):
        need(x['order']==y['n'] and x['classes']==y['ordinary_classes'],'native component counts')
        need(x['tc_histogram']==y['tc_histogram'],'entire native tc histogram')
    # Independent permutation-oracle exhaustions on all labelled tournaments
    # of size<=5, using lexicographic all-permutation encodings.
    for n in range(1,6):
        pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
        keys=set()
        for code in range(1<<len(pairs)):
            a=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(pairs): a[i][j]=(code>>k)&1; a[j][i]=1-a[i][j]
            keys.add(min(sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs)) for p in permutations(range(n))))
        need(len(keys)==result['results'][n-1]['classes'],'labelled orbit oracle')
    # Rank-15 attachment hostile. Verify the whole explicit wreath subgroup;
    # exactness follows from preserved block classes and Aut(T5)=C5.
    t5=tuple(sum(1<<((i+d)%5) for d in (1,2)) for i in range(5))
    t15=tuple((t5[i]<<(5*k))|(31<<(5*((k+1)%3))) for k in range(3) for i in range(5))
    M=tuple(31<<(5*(i//5)) for i in range(15))
    aut=0; kernel=0
    for root_shift in range(3):
        for shifts in product(range(5),repeat=3):
            p=[5*((k+root_shift)%3)+(i+shifts[k])%5 for k in range(3) for i in range(5)]
            need(all((t15[u]>>v&1)==(t15[p[u]]>>p[v]&1) for u in range(15) for v in range(15)),'wreath automorphism')
            aut+=1
            kernel+=all((M[u]>>v&1)==(M[u]>>p[v]&1) for u in range(15) for v in range(15))
    need(kernel==125 and aut==375,'fifteen chart kernel')
    need(all(row.bit_count()==5 for row in M),'equitable row')
    need(all(sum(row>>j&1 for row in M)==5 for j in range(15)),'equitable column')
    raw=gzip.decompress((HERE/'tournament_bank.json.gz').read_bytes())
    tournaments=json.loads(raw)
    need(tournaments['order']==9 and len(tournaments['rows'])==191536,'n9 bank')
    for rows in tournaments['rows']:
        need(len(rows)==9 and all(not (rows[i]>>i&1) for i in range(9)) and all((rows[i]>>j&1)+(rows[j]>>i&1)==1 for i,j in __import__('itertools').combinations(range(9),2)),'bank tournament type')
    out={'status':'PASS','checks':checks,'native_metagraph_crosscheck_orders':[1,8],'full_labelled_permutation_oracle_orders':[1,5],'arithmetic_triples_integer_checked':698027,'rank15_chart_kernel':kernel,'rank15_wreath_subgroup':aut,'bank_order':9,'bank_classes':len(tournaments['rows']),'bank_uncompressed_sha256':hashlib.sha256(raw).hexdigest()}
    (HERE/'audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out))

if __name__=='__main__': main()
