"""Two minimal forbidden induced cores, with witnesses for all proper subsets."""
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def paths(rows):
    n=len(rows); full=(1<<n)-1
    e=[[0]*n for _ in range(full+1)]; h=[1]+[0]*full
    for mask in range(1,full+1):
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                e[mask][v]=sum(e[rest][u] for u in range(n) if rows[u]>>v&1) if rest else 1
        h[mask]=sum(e[mask])
    return e,h

def profile(rows):
    n=len(rows); full=(1<<n)-1
    e,h=paths(rows)
    starts,_=paths([full^(1<<v)^rows[v] for v in range(n)])
    fs=[sum(e[m][r]*h[full^m] for m in range(1,full+1)) for r in range(n)]
    return h[full],starts[full],fs

def strong(rows):
    full=(1<<len(rows))-1
    for r in range(len(rows)):
        reached=1<<r
        while True:
            newer=reached
            for v,row in enumerate(rows):
                if reached>>v&1: newer|=row
            if newer==reached: break
            reached=newer
        if reached!=full: return False
    return True

def induced(rows,mask):
    keep=[v for v in range(len(rows)) if mask>>v&1]
    return [sum(1<<j for j,v in enumerate(keep) if rows[u]>>v&1) for u in keep]

def birth(rows,r):
    n=len(rows)
    return [row | ((1<<n) if v==r else 0) for v,row in enumerate(rows)]+[((1<<n)-1)^(1<<r)]

records=[]
for name,rows in [('R5',[6,12,24,17,3]),('C3_TT222',[14,12,56,48,35,3])]:
    h,s,f=profile(rows)
    need(strong(rows) and min(s)>=2 and min(x-h for x in f)>=2,'full forbidden criterion')
    cases=[]
    for mask in range(1,(1<<len(rows))-1):
        small=induced(rows,mask)
        hh,ss,ff=profile(small)
        if strong(small) and len(small)>1:
            need(1 in ss,'proper strong unique start')
            w=ss.index(1); target=birth(small,w)
            th,ts,tf=profile(target)
            need(tf[-1]-th==1,'actual strong ready witness')
            cases.append({'mask':mask,'rows':small,'H':hh,'unique_start':w,
                          'witness_rows':target,'witness_root':len(small),'witness_H':th})
        else:
            need(hh==1,'all other proper restrictions transitive')
            target=small if len(small)>1 else birth(small,0)
            th,ts,tf=profile(target)
            roots=[r for r,value in enumerate(tf) if value-th==1]
            need(roots,'actual transitive ready witness')
            cases.append({'mask':mask,'rows':small,'H':hh,'transitive':True,
                          'witness_rows':target,'witness_root':roots[0],'witness_H':th})
    if name=='C3_TT222':
        need(h==45 and s==[8,7,8,7,8,7] and [x-h for x in f]==[13,26,13,26,13,26], 'six core profile')
        five=[induced(rows,((1<<6)-1)^(1<<v)) for v in range(6)]
        need(all(any(row.bit_count()!=2 for row in t) for t in five),'regular five core absent')
        need(all(profile(t)[0]==11 for t in five),'five-deletion counts')
    records.append({'name':name,'rows':rows,'H':h,'starts':s,'readiness_defects':[x-h for x in f],
                    'proper_subset_witnesses':cases})
out={'status':'FINITE-EXACT full profiles and all92 nonempty proper induced subsets',
     'scope':'minimal forbidden cores use the proved all-order criterion; proper-subset embeddings are literal',
     'records':records,
     'proof_sha256':hashlib.sha256((HERE/'ROOTED_HEREDITY.md').read_bytes()).hexdigest()}
suffix='.optimized' if sys.flags.optimize else ''
(HERE/('forbidden_cores'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps([(r['name'],r['H'],len(r['proper_subset_witnesses'])) for r in records]))
