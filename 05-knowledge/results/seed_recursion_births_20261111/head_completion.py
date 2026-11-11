"""Complete fixed-kernel decisions via the inherited strong H<=375 inventory.

The computation is new; the all-order completeness of that inventory is an
explicit inherited input. Exterior SCCs multiply every Hamilton coefficient.
"""
from pathlib import Path
from functools import lru_cache,reduce
from math import gcd
import hashlib,json
import head_search as hs

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SOURCE=HERE.parent/'session173_modular_clocks_20261103/counts527/inventory_375.json'


def join(left,right):
    n=len(left)
    return [row | (((1<<len(right))-1)<<n) for row in left]+[row<<n for row in right]


def value(coeff,a):
    return sum(c*hs.ib.onto(a,j) for j,c in enumerate(coeff,1))


def shift(coeff,s):
    result=[]
    for a in range(1,len(coeff)+1):
        weights=[hs.ib.onto(a,j) for j in range(1,a+1)]
        numerator=value(coeff,a+s)-sum(c*w for c,w in zip(result,weights))
        hs.ib.tc.require(numerator%weights[-1]==0,'shift coefficient integral')
        result.append(numerator//weights[-1])
    return hs.trim(result)


def main():
    raw=SOURCE.read_bytes().replace(b'\r\n',b'\n')
    source=json.loads(raw)
    records=[dict(n=int(n),code=int(row[0]),H=row[1]) for n,items in source['orders'].items() for row in items]
    records.sort(key=lambda r:(r['n'],r['code']))
    by_h={}
    for item in records:by_h.setdefault(item['H'],[]).append(item)
    @lru_cache(None)
    def witness(h):
        if h==1:return ()
        if h in by_h:
            item=by_h[h][0]
            return tuple(hs.ib.cc.decode(item['n'],item['code']))
        for d in sorted(by_h):
            if d>=h:break
            if h%d==0 and witness(h//d) is not None:
                return tuple(join(witness(d),witness(h//d)))
        return None
    queries={}
    needed=set()
    for name,coeff in hs.TARGETS.items():
        options=[]
        s=0
        while value(coeff,s+1)<=source['limit']:
            target=shift(coeff,s)
            content=reduce(gcd,target)
            factors=[g for g in range(1,content+1) if content%g==0 and witness(g) is not None]
            options.append(dict(shift=s,coefficients=target,factors=factors))
            for g in factors:needed.add(tuple(c//g for c in target))
            s+=1
        queries[name]=options
    cache,checked={},0
    for item in records:
        if item['H'] not in {target[0] for target in needed}:continue
        rows=hs.ib.cc.decode(item['n'],item['code'])
        hs.ib.tc.require(hs.ib.cc.hp(rows)==item['H'],'inventory code/count convention')
        for v in range(item['n']):
            coeff=tuple(hs.kernel(rows,v))
            checked+=1
            if coeff in needed and coeff not in cache:
                cache[coeff]=dict(**item,rows=rows,port=v,coefficients=list(coeff))
    decisions={}
    for name,options in queries.items():
        rows=[]
        for option in options:
            candidates=[]
            for g in option['factors']:
                target=tuple(c//g for c in option['coefficients'])
                if target not in cache:continue
                core=cache[target]
                exterior=list(witness(g))
                base=join(exterior,core['rows'])
                port=len(exterior)+core['port']
                hs.ib.tc.require(hs.kernel(base,port)==option['coefficients'],'joined witness kernel')
                candidates.append(dict(exterior_H=g,rows=base,port=port,order=len(base),core=core))
            chosen=min(candidates,key=lambda x:(x['order'],x['exterior_H'])) if candidates else None
            rows.append(dict(**option,witness=chosen))
        decisions[name]=rows
    heads={}
    for name,options in decisions.items():
        positive=next((option for option in options if option['witness'] is not None),None)
        if positive:
            prefix=[]
            for r in range(1,positive['shift']+1):
                h=value(hs.TARGETS[name],r)
                graph=witness(h)
                hs.ib.tc.require(graph is not None and hs.ib.cc.hp(graph)==h,'finite head witness')
                prefix.append(dict(parameter=r,H=h,rows=list(graph)))
            heads[name]=dict(start=positive['shift']+1,kernel=positive['coefficients'],
                             witness=positive['witness'],prefix=prefix)
    report=dict(inventory_source=str(SOURCE.relative_to(ROOT)).replace('\\','/'),
                inventory_sha256_lf=hashlib.sha256(raw).hexdigest(),inventory_limit=source['limit'],
                scope='global fixed-kernel decisions conditional on inherited complete strong inventory',
                strong_ports_checked=checked,decisions=decisions,completed_heads=heads,
                unresolved_heads=[name for name in hs.TARGETS if name not in heads])
    (HERE/'head_completion.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'strong_ports_checked':checked,'completed':{k:v['start'] for k,v in heads.items()},
                      'unresolved':report['unresolved_heads']},indent=2))


if __name__=='__main__':main()
