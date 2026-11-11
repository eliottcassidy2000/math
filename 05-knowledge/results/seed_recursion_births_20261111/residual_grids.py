"""All positive power-of-two constant binary slopes in the 95-port bank."""
from pathlib import Path
from math import gcd
import json,sys
import port_transport as pt
import head_completion as hc

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_coverage_compiler_20261111'))
import verify as independent
import construct as old
import realize35


def main():
    records=json.loads((HERE.parent/'endpoint_coverage_compiler_20261111/inflation_bank_7.json').read_text())['records']
    binary=[]
    for item in records:
        if len(pt.cp.components(item['rows']))!=1:continue
        coeff=list(item['slope_cover_coefficients'])
        while coeff and coeff[-1]==0:coeff.pop()
        if len(coeff)!=2:continue
        A=coeff[1];B=coeff[0]-2*A
        if B<=0 or B&(B-1):continue
        d=A;v=0
        while d%2==0:d//=2;v+=1
        binary.append(dict(n=item['n'],code=item['code'],rows=item['rows'],selected=item['selected'],
                           port=item['inflated_vertex'],A=A,B=B,d=d,valuation_B=B.bit_length()-1,
                           block_size_for_parameter_shift=B.bit_length()-1-v))
    item=next(x for x in binary if x['d']==35)
    values=[]
    for a in range(1,item['n']+1):
        labels,rows=pt.ib.inflate(item['rows'],item['port'],a)
        selected=[i for i,v in enumerate(labels) if item['selected']>>v&1]
        h,w,u=independent.endpoint_state(rows,selected)
        independent.require(h-w==1 and u==35*2**a+32,'d35 literal seed')
        values.append(u)
    independent.require(independent.invert(values)==[102,35]+[0]*5,'d35 full independent slope kernel')
    for r in [1,2]:
        for s in [6,7]:
            labels,rows=pt.ib.inflate(item['rows'],item['port'],r+5)
            selected=[i for i,v in enumerate(labels) if item['selected']>>v&1]
            out=old.tail(rows,selected,s-5)
            independent.require(independent.count(out)[0]==1+(35*2**r+1)*2**s,'literal d35 ray')
    # Every below-head shift with initial H<=375: two for s=1, one for s=2.
    queries=[]
    for s in range(1,6):
        q=0
        while 1+2**s+35*2**(1+q+s)<=375:
            A=35*2**(q+s);target=[1+2**s+2*A,A]
            content=gcd(*target)
            factors=[g for g in [1,3,5] if content%g==0]
            # These target contents are 1 or 5, so this is exhaustive.
            independent.require(content in [1,5],'all exterior divisors included')
            queries.append(dict(exponent=s,shift=q,kernel=target,factors=factors,witness=None))
            q+=1
    inventory=json.loads(hc.SOURCE.read_text());tested=0
    for n,entries in inventory['orders'].items():
        for record in entries:
            relevant=[(q,g) for q in queries for g in q['factors'] if record[1]*g==q['kernel'][0]]
            if not relevant:continue
            rows=pt.ib.cc.decode(int(n),int(record[0]))
            independent.require(pt.ib.cc.hp(rows)==record[1],'inventory decode')
            for v in range(int(n)):
                kernel=hc.hs.kernel(rows,v);tested+=1
                for q,g in relevant:
                    if q['witness'] is None and [g*x for x in kernel]==q['kernel']:
                        q['witness']=dict(core_rows=rows,core_port=v,exterior_H=g)
    verified_heads=0
    for q in queries:
        if q['witness'] is None:continue
        w=q['witness'];values=[]
        # Between repeated port occurrences are at least two old vertices:
        # v -> x -> v is impossible. Therefore j <= floor((n+1)/2).
        limit=(len(w['core_rows'])+1)//2
        for a in range(1,limit+1):
            _,rows=pt.ib.inflate(w['core_rows'],w['core_port'],a)
            values.append(independent.count(rows)[0])
        coeff=independent.invert(values)
        independent.require([w['exterior_H']*c for c in coeff]==q['kernel']+[0]*(limit-2),'independent d35 head kernel')
        verified_heads+=1
    report=dict(scope='all 95 strong stable ports, slopes A*2^a+2^j only',binary_records=binary,
                odd_residuals=sorted(set(x['d'] for x in binary)),d35_independent_seed_sizes=7,
                d35_literal_tail_controls=4,head_strong_ports_tested=tested,head_queries=queries,
                independently_verified_head_kernels=verified_heads)
    for r in [1,2]:
        for s in [1,6,7]:
            item=realize35.construct(r,s,report)
            independent.require(independent.count(item['rows'])[0]==item['expected_H'],'d35 emitted graph')
    target=HERE/'residual_grids.json'
    encoded=json.dumps(report,indent=2)+'\n'
    if target.exists() and 'independently_verified_head_kernels' in json.loads(target.read_text()):
        independent.require(json.loads(target.read_text())==report,'residual replay unchanged')
    target.write_text(encoded,encoding='utf8')
    print(json.dumps(dict(odd_residuals=report['odd_residuals'],head_ports=tested,
                         head_decisions=[dict(exponent=q['exponent'],shift=q['shift'],witness=q['witness']) for q in queries]),indent=2))


if __name__=='__main__':main()
