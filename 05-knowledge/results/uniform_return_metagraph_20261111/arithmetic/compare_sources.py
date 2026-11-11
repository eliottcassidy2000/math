"""Compare both enumerators to pinned OEIS numerical tables.

The default replay is offline. --refresh deliberately replaces the primary
numeric source files; no OEIS update or submission is made by this script.
"""
from hashlib import sha256
from pathlib import Path
from urllib.request import Request,urlopen
import argparse
import json
import maxima
import quadratic_primes

HERE=Path(__file__).resolve().parent
SOURCES={
    'A135244':'https://oeis.org/A135244/b135244.txt',
    'A002496':'https://oeis.org/A002496/b002496.txt',
}


def require(ok,msg):
    if not ok:raise RuntimeError(msg)


def read_table(data):
    rows=[]
    for line in data.decode().splitlines():
        fields=line.split()
        if fields and not fields[0].startswith('#'):
            require(len(fields)==2,'unexpected source row')
            rows.append(tuple(map(int,fields)))
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--refresh',action='store_true')
    args=ap.parse_args()
    reports={}
    tables={}
    for name,url in SOURCES.items():
        path=HERE/(name.lower()+'_source.txt')
        if args.refresh:
            data=urlopen(Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=30).read()
            path.write_bytes(data)
        data=path.read_bytes()
        rows=read_table(data)
        tables[name]=rows
        reports[name]={'url':url,'bytes':len(data),'sha256':sha256(data).hexdigest(),
                       'rows':len(rows),'first':rows[0],'last':rows[-1]}
    X=tables['A135244'][-1][0]
    vals,_,_=maxima.enumerate_maxima(X)
    require(all(vals[n]==v for n,v in tables['A135244']),'A135244 published table mismatch')
    last=tables['A002496'][-1][1]
    from math import isqrt
    N=isqrt(last-1)
    require(N*N+1==last,'square-plus-one endpoint')
    least,_=quadratic_primes.quadratic_sieve(N)
    primes=[2]+[r*r+1 for r in range(2,N+1,2) if not least[r]]
    require(list(enumerate(primes,1))==tables['A002496'],'A002496 published table mismatch')
    result={'status':'PASS','retrieved_utc_date':'2026-11-11','sources':reports,
            'scope':'Primary numerical prefixes agree; no new OEIS term or record claim.'}
    encoded=json.dumps(result,indent=2,sort_keys=True)+'\n'
    (HERE/'source_comparison.json').write_text(encoded,encoding='utf-8')
    print(encoded)


if __name__=='__main__':main()
