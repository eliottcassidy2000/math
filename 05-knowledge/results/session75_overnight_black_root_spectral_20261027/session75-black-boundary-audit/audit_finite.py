#!/usr/bin/env python3
"""Second independent complete replay, restoring every converse-built route
back to its original named source before native K conversion. No producer imports.
"""
import csv,gzip,pathlib,json,hashlib,collections
from audit_checks import *
H=pathlib.Path(__file__).resolve().parent

def decode(n,code):
    # Read binary word from left to right (MSB-first upper triangle).
    word=f'{code:0{n*(n-1)//2}b}';need(len(word)==n*(n-1)//2,'extra source bits')
    a=[[False]*n for _ in range(n)];k=0
    for i in range(n):
        for j in range(i+1,n):a[i][j]=word[k]=='1';a[j][i]=not a[i][j];k+=1
    return tuple(tuple(row) for row in a)

def main():
    reports={};steps=0;record_hash=hashlib.sha256();literal_hash=hashlib.sha256()
    with gzip.open(H/'inputs/witnesses.jsonl.gz','rt') as certs:
        for n in range(4,10):
            stats=collections.Counter();dist=collections.Counter();p=H/'inputs'/(f'classes_n{n}.tsv' if n<9 else 'n9_classes.tsv.gz')
            with (gzip.open if n==9 else open)(p,'rt') as f:
                for ix,source in enumerate(csv.DictReader(f,delimiter='\t')):
                    line=next(certs);record_hash.update(line.encode());c=json.loads(line);code=int(source['canon' if n<9 else 'code']);need((c['n'],c['i'],c['c'])==(n,ix,code),'source record misbound')
                    a=decode(n,code);stats['ordinary_classes']+=1;st=strong(a)
                    if not st:
                        need(c['t']=='N','nonstrong filter mismatch');stats['nonstrong']+=1;continue
                    stats['strong']+=1
                    if c['t']=='P':
                        hp=list(paths(a));need(hp and all(anti_failure(a,p) is None for p in hp),'false pure-blue exception');need(n==5 and len(hp)==15,'unexpected pure-blue exception');stats['pure_blue']+=1;continue
                    need(c['t']=='W','strong source omitted');lam,delta,_,_=cut_certificate(a);need(lam==delta and lam>0,'wrong finite lambda/delta equality');need(len(c['p'])==delta,'wrong geodesic length')
                    dual=c['d'];v=c['v'];raw=[tuple(p) for p in c['p']]
                    need(all(p[-1]==v for p in raw),'certificate lost named fixed endpoint')
                    # Do not silently change source to its converse.
                    pp=[p[::-1] if dual else p for p in raw]
                    ds=degrees(a);need((ds[v] if dual else n-1-ds[v])==delta,'port not extremal in original source')
                    terminal,moves=replay(a,pp,require_core=True)
                    need(v in core(terminal),'wrong terminal named singleton')
                    need(degrees(terminal)[v]==(0 if dual else n-1),'wrong terminal source/sink gauge')
                    dist[delta]+=1;stats['strong_black_active']+=1;stats['literal_R_edges']+=len(pp);stats['literal_K_edges']+=len(pp);steps+=len(pp)
                    literal_hash.update(json.dumps({'n':n,'code':code,'dual_construction':dual,'original_named_rows':rows(a),'moves':moves,'terminal_rows':rows(terminal)},sort_keys=True,separators=(',',':')).encode())
            need(stats['ordinary_classes']=={4:4,5:12,6:56,7:456,8:6880,9:191536}[n],'ordinary source count')
            need(stats['strong']=={4:1,5:6,6:35,7:353,8:6008,9:178133}[n],'strong source count')
            reports[str(n)]={**dict(stats),'distance_counts':dict(sorted(dist.items()))}
        need(next(certs,None) is None,'extra certificate')
    c3=from_rows([2,4,1]);p3=list(paths(c3));need(strong(c3) and len(p3)==3 and all(anti_failure(c3,p) is None for p in p3),'C3 exception failure')
    print(json.dumps({'status':'PASS','scope':'every ordinary strong non-pure-blue source n4..9, with C3 boundary control; native K starts at original named source even for converse-built certificates','orders':reports,'source_certificate_uncompressed_sha256':record_hash.hexdigest(),'literal_replay_summary_sha256':literal_hash.hexdigest(),'literal_R_edges':steps,'literal_K_edges':steps,'conclusion':'both target distances equal lambda=delta on every eligible source through n9; no converse-merging'},sort_keys=True,indent=2))
if __name__=='__main__':main()
