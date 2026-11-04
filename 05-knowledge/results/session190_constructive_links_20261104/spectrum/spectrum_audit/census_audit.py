#!/usr/bin/env python3
from pathlib import Path
from hashlib import sha256
from itertools import permutations, product
from collections import defaultdict
import json,subprocess,time
D=Path(__file__).resolve().parent
P=D.parent/'spectrum'
if not P.exists():P=D.parent/'session190-spectrum'
COVER=D/'inherited_cover8.json'

def need(x,msg):
    if not x:raise RuntimeError(msg)
def canon(r):
    groups=defaultdict(list)
    for i,z in enumerate(r):groups[z.bit_count()].append(i)
    keys=sorted(groups)
    best=None
    for blocks in product(*(permutations(groups[k]) for k in keys)):
        p=sum(blocks,())
        # Encode upper triangle as one integer, independent from producer tuple minimum.
        code=0;bit=0
        for i in range(len(r)):
            for j in range(i+1,len(r)):
                code|=((r[p[i]]>>p[j])&1)<<bit;bit+=1
        if best is None or code<best:best=code
    return best

def main():
    start=time.time();bank=json.loads((P/'bank8.json').read_text())
    need(sha256((P/'bank8.json').read_bytes()).hexdigest()=='55409706be4f8c71f5e2b2cd0777c9c0f0f8de97eac4b174d648f50b638dc903','bank8 pin')
    # Recheck induction witnesses from the inherited complete rooted extension cover.
    if COVER.exists():
        need(sha256(COVER.read_bytes()).hexdigest()=='7af3e95ebe2c20d594f8e7623a10fa08682a4963e6a5d1b61b9600d56ea06175','inherited cover pin')
        cover=json.loads(COVER.read_text());levels=cover['levels'];transports=0
        need(levels[0]['representatives']==[[]],'empty base')
        for n in range(1,9):
            prev=levels[n-1]['representatives'];cur=levels[n]['representatives'];records=levels[n]['extension_records']
            need(len(records)==len(prev)*(1<<(n-1)),'all extensions')
            for i,(target,packed) in enumerate(records):
                parent,mask=divmod(i,1<<(n-1));r=prev[parent];g=[a|((1-((mask>>j)&1))<<(n-1)) for j,a in enumerate(r)]+[mask]
                p=[(packed>>(3*j))&7 for j in range(n)]
                need(sorted(p)==list(range(n)) and packed>>(3*n)==0,'transport bijection')
                need(0<=target<len(cur),'target index');t=cur[target]
                need(all(((g[p[a]]>>p[b])&1)==((t[a]>>b)&1) for a in range(n) for b in range(n)),'edgewise isomorphism')
                transports+=1
        need(bank==levels[8]['representatives'],'exact bank output')
    else:raise RuntimeError('Inherited cover required for audit; use pinned prior bank-completeness certificate in publication')
    print('PASS inherited cover: '+str(transports)+' transports',flush=True)
    subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',str(D/'independent_census.cpp'),'-o',str(D/'independent_census')],check=True)
    def run(mode,n=None,rows=None):
        text=None if rows is None else f'{n} {len(rows)}\n'+''.join(' '.join(map(str,r))+'\n' for r in rows)
        result=subprocess.check_output([str(D/'independent_census'),mode],input=text,text=True)
        label='small' if n is None else mode+str(n)
        (D/(label+'.out')).write_text(result)
        print('PASS '+label,flush=True)
        return result
    small=run('small');base=run('bank',8,bank);nine=run('extend',8,bank)
    lowrows=[tuple(map(int,l.split()[1:])) for l in nine.splitlines() if l.startswith('LOW ')]
    canonical={}
    for r in lowrows:canonical.setdefault(canon(r),r)
    need(len(canonical)==545,'545 low nine classes')
    bank9=sorted(canonical.values());(D/'independent_low9_bank.json').write_text(json.dumps(bank9,separators=(',',':'))+'\n')
    # Check that both banks represent exactly the same classes without requiring same gauges.
    producer9=json.loads((P/'low9_bank.json').read_text())
    need({canon(r) for r in producer9}==set(canonical),'independent canonical bank agrees')
    ten=run('extend',9,bank9)
    expected=json.loads((P/'census_receipt.json').read_text())
    need(small.splitlines()==expected['small'],'small outputs match')
    summaries=[next(l for l in s.splitlines() if l.startswith('SUMMARY')) for s in [base,nine,ten]]
    need(summaries==expected['summaries'],'all summary counts match')
    targets=[{'order':n,'rows':[l for l in raw.splitlines() if l.startswith('TARGET')]} for n,raw in [(8,base),(9,nine),(10,ten)]]
    need(targets==expected['targets'],'complete target F sets and multiplicities match')
    result=dict(status='PASS',method='Independent actual Hamilton-word traversal with safe global cutoff at 202, actual partial-word enumeration for every candidate profile, Floyd-Warshall strong test, integer-edge-code degree-preserving canon.',inherited_extension_transports=transports,small=small.splitlines(),summaries=summaries,low9_classes=len(bank9),targets=targets)
    (D/'census_audit_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS complete finite census',flush=True)
if __name__=='__main__':main()
