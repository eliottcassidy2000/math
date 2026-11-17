#!/usr/bin/env python3
"""Finite exact Chan–Pak birth decision by divisor-only strong search.
No network or third-party packages. Budget exhaustion is NEVER nonexistence.
"""
import argparse,csv,hashlib,json,math
from pathlib import Path

def busch_floor(n):
    if n<3:raise ValueError('nontrivial strong order must be >=3')
    q,r=divmod(n,3)
    return (3*5**(q-1),5**q,9*5**(q-1))[r]

def strong_cap(h):
    n=3
    while busch_floor(n)<=h:n+=1
    return n-1

def validate(rows):
    n=len(rows)
    if n<1 or any(type(r)!=int or r<0 or r>=1<<n or r>>i&1 for i,r in enumerate(rows)):
        raise ValueError('invalid tournament rows')
    if any(((rows[i]>>j)&1)+((rows[j]>>i)&1)!=1 for i in range(n) for j in range(i)):
        raise ValueError('rows do not orient every pair exactly once')

def paths(rows,cap=None):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                if rest:
                    x=sum(dp[rest][u] for u in range(n) if rest>>u&1 and rows[u]>>v&1)
                    dp[mask][v]=min(x,cap+1) if cap is not None else x
    x=sum(dp[-1]);return min(x,cap+1) if cap is not None else x

def strong(rows):
    n=len(rows)
    for start in (0,):
        for reverse in (False,True):
            seen=1<<start
            while True:
                new=seen
                for i in range(n):
                    if seen>>i&1:
                        new |= sum(((rows[j]>>i)&1)<<j for j in range(n)) if reverse else rows[i]
                if new==seen:break
                seen=new
            if seen!=(1<<n)-1:return False
    return True

def labelled(n,code):
    rows=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:rows[i]|=1<<j
            else:rows[j]|=1<<i
            k+=1
    return rows

def join(a,b):
    n=len(a);m=len(b)
    return [r|(((1<<m)-1)<<n) for r in a]+[r<<n for r in b]

def divisors(h):
    low=[];high=[]
    for d in range(1,math.isqrt(h)+1):
        if h%d==0:
            low.append(d)
            if d*d!=h:high.append(h//d)
    return low+high[::-1]

def load_atlas(path,allowed=None):
    if not path:return {},None
    p=Path(path);raw=p.read_bytes();out={}
    if p.suffix.lower()=='.json':
        data=json.loads(raw)
        entries=({'h':h,'outgoing_masks':v} for h,v in data.items()) if isinstance(data,dict) else data
    else:entries=csv.DictReader(raw.decode().splitlines(),delimiter='\t')
    for e in entries:
        h=int(e['h'])
        if allowed is not None and h not in allowed:continue
        rows=e['outgoing_masks']
        if isinstance(rows,str):rows=[int(x) for x in rows.split(',')]
        validate(rows)
        if paths(rows)!=h:raise ValueError(f'atlas witness has wrong H for {h}')
        if h not in out or len(rows)<len(out[h]):out[h]=rows
    return out,hashlib.sha256(raw).hexdigest()

def solve(h,max_labelled=1000000,atlas=None):
    if h<1:raise ValueError('h must be positive')
    base={'target':h,'method':'Busch strong-order caps and exact divisor-first SCC recursion',
          'budget':max_labelled,'tested_labelled':0,'certificate':[]}
    if h==1:return dict(base,status='minimum',order=1,outgoing_masks=[0],reason='unique transitive first birth')
    if h%2==0:return dict(base,status='unattainable',reason='every tournament has an odd Hamiltonian path count')
    ds=divisors(h)
    upper,sha=load_atlas(atlas,set(ds));base['atlas_sha256']=sha
    core_caps={1:0}
    for d in ds[1:]:
        core_caps[d]=max(strong_cap(f)+core_caps[d//f] for f in divisors(d)[1:])
    base['Busch_factor_core_cap']=core_caps[h]
    known={1:[]}
    for d in ds[1:]:
        best=upper.get(d);origin='verified atlas upper bound' if best is not None else None
        for f in divisors(d)[1:-1]:
            if known.get(f) is not None and known.get(d//f) is not None:
                candidate=join(known[f],known[d//f])
                if best is None or len(candidate)<len(best):best=candidate;origin=f'ordinal product {f}*{d//f}'
        cap=strong_cap(d);record={'divisor':d,'Busch_cap':cap,'orders':[]}
        for n in range(3,min(cap,len(best)-1 if best is not None else cap)+1):
            if math.factorial(n)<d:
                record['orders'].append({'order':n,'status':'excluded by H<=n!'})
                continue
            total=1<<(n*(n-1)//2);found=False
            for code in range(total):
                if max_labelled and base['tested_labelled']>=max_labelled:
                    record['orders'].append({'order':n,'status':'partial','tested':code,'total':total})
                    base['certificate'].append(record)
                    result=dict(base,status='incomplete',blocked_divisor=d,reason='labelled-search budget exhausted; no nonexistence claim')
                    candidate=upper.get(h)
                    if d==h and best is not None:candidate=best
                    if candidate is not None:result.update(witness_upper_order=len(candidate),outgoing_masks=candidate)
                    return result
                rows=labelled(n,code);base['tested_labelled']+=1
                if strong(rows) and paths(rows,d)==d:
                    best=rows;origin='literal strong witness';found=True
                    record['orders'].append({'order':n,'status':'witness','tested':code+1,'total':total})
                    break
            if found:break
            record['orders'].append({'order':n,'status':'exhausted','tested':total,'total':total})
        known[d]=best
        record.update(status='minimum' if best is not None else 'unattainable',minimum_order=len(best) if best is not None else None,source=origin)
        base['certificate'].append(record)
    witness=known[h]
    if witness is None:return dict(base,status='unattainable',reason='all strong divisors exhausted through exact Busch caps; no SCC factorization exists')
    validate(witness)
    if paths(witness)!=h:raise ArithmeticError('final witness failed independent recount')
    return dict(base,status='minimum',order=len(witness),outgoing_masks=witness)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('h',type=int)
    p.add_argument('--max-labelled',type=int,default=1000000,help='total labelled graphs tested; 0 means unlimited (may be impractical)')
    p.add_argument('--atlas',help='optional JSON h->rows or TSV h/order/outgoing_masks; witnesses verified, minima not trusted')
    p.add_argument('--output',help='save JSON result; otherwise stdout')
    a=p.parse_args()
    if a.max_labelled<0:p.error('--max-labelled must be nonnegative')
    try:r=solve(a.h,a.max_labelled,a.atlas)
    except (ValueError,KeyError,TypeError) as e:p.error(str(e))
    text=json.dumps(r,indent=2)+'\n'
    if a.output:Path(a.output).write_text(text)
    else:print(text,end='')
    return 2 if r['status']=='incomplete' else 0
if __name__=='__main__':raise SystemExit(main())
