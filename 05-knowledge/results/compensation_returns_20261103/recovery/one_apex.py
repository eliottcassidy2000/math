"""Complete bounded core test proving arbitrary-order unit-deletion exclusion."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
RESULTS=PACKET.parent
COUNTER=PACKET/'audit/exact_counter.cpp'
OLD=RESULTS/'c4_port_returns_20261103/entry_selector/receipt.json'


def need(ok,why):
    if not ok:raise ValueError(why)


def rows_from_bits(bits):
    m=len(bits)
    return [sum(1<<(i+1) for i,b in enumerate(bits) if not b)]+[
        b+sum(1<<(j+1) for j in range(i+1,m)) for i,b in enumerate(bits)]


def split_formula(bits):
    if not bits:return 1
    m=len(bits)
    return (1-bits[0])+bits[-1]+sum(bits[i]*(1-bits[i+1]) for i in range(m-1))+sum(
        (1-bits[j])*bits[i]*(1<<(i-j-1)) for i in range(m) for j in range(i))


def partition_count(bits):
    m=len(bits);full=(1<<m)-1
    if not m:return 1
    result=0
    for before in range(1<<m):
        after=full^before
        if before and not bits[before.bit_length()-1]:continue
        if after and bits[(after&-after).bit_length()-1]:continue
        result+=1
    return result


def strong(rows):
    full=(1<<len(rows))-1
    for root in range(len(rows)):
        seen=1<<root;front=seen
        while front:
            bit=front&-front;front^=bit;fresh=rows[bit.bit_length()-1]&~seen;seen|=fresh;front|=fresh
        if seen!=full:return False
    return True


def deletion_word(rows,names,apex):
    keep=[i for i,a in enumerate(names) if a!=apex]
    order=sorted(keep,key=lambda i:-sum(rows[i]>>j&1 for j in keep))
    need(all(rows[a]>>b&1 for k,a in enumerate(order) for b in order[k+1:]),'actual deletion is transitive')
    w=names.index(apex)
    return [names[i] for i in order],[int(rows[i]>>w&1) for i in order]


def main():
    targets={1409,1427};words=[];census=[];typed=0
    for n in range(3,14):
        m=n-1;values=[]
        for middle in range(1<<(m-2)):
            bits=[0]+[(middle>>i)&1 for i in range(m-2)]+[1]
            rows=rows_from_bits(bits);h=split_formula(bits)
            need(strong(rows),'strong endpoint characterization')
            need(h==partition_count(bits),'independent before/after split count')
            need(h>=2+(1<<(n-3)),'elementary exponential core bound')
            words.append(dict(n=n,word=''.join(map(str,bits)),H=h,rows=rows));values.append(h)
            typed+=bits[1]==1
        census.append(dict(order=n,words=len(values),distinct_H=len(set(values)),minimum_H=min(values),maximum_H=max(values),target_hits=sorted(set(values)&targets)))
        if n>=5:
            need(min(values)==((1<<((n-3)//2))+1)*((1<<((n-2)//2))+1),'inherited sharp one-apex minimum')
    need(len(words)==2047 and typed==1024,'complete strong and typed word universes')
    need(not any(w['H'] in targets for w in words),'unit-deletion target exists')
    weak_controls=[]
    for m in range(1,6):
        for mask in range(1<<m):
            bits=[mask>>i&1 for i in range(m)]
            h=split_formula(bits);need(h==partition_count(bits),'weak boundary split formula')
            weak_controls.append(dict(n=m+1,rows=rows_from_bits(bits),H=h))
    source=json.loads(OLD.read_text());debts=[];controls=[]
    for record in source['records'][:2]:
        child=record['stages'][-1]['graph'];repair=record['neutral_triangle_tests'][0]
        need(repair['strong'],'declared strong debt source')
        g=repair['graph'];order,bits=deletion_word(g['rows'],g['names'],record['entry_deletion_mark'])
        need(split_formula(bits)==repair['H'],'debt source one-apex formula')
        prepared=record['stages'][1]['graph'];first_birth=record['stages'][2]['graph']
        _,prepared_bits=deletion_word(prepared['rows'],prepared['names'],record['entry_deletion_mark'])
        need(prepared_bits[-2:]==[1,1],'prepared terminal00 reader word')
        weight=sum((1-b)*(1<<(len(prepared_bits)-i-2)) for i,b in enumerate(prepared_bits[:-1]))
        prepared_h=split_formula(prepared_bits)
        need(bits==prepared_bits[:-1]+[0,0,1],'actual strong repair inserts two zero letters')
        need(repair['H']==prepared_h+3*weight+4,'exact strong-repair count law')
        need(record['C2']['target_H']-repair['H']==2*prepared_h-3*weight-2,'exact30/66 debt law')
        debts.append(dict(family=record['family'],target_H=record['C2']['target_H'],source_H=repair['H'],
                          required_compensation=record['C2']['target_H']-repair['H'],apex=record['entry_deletion_mark'],root=99,
                          prepared_H=prepared_h,prepared_word=''.join(map(str,prepared_bits)),weighted_zero_sum=weight,
                          first_birth_H=prepared_h+weight+2,
                          deletion_order=order,deletion_word=''.join(map(str,bits)),rows=g['rows'],names=g['names']))
        controls.extend([dict(kind='child',family=record['family'],graph=child,H=record['C2']['target_H'],q=3,apex=record['entry_deletion_mark']),
                         dict(kind='strong_q1_repair',family=record['family'],graph=g,H=repair['H'],q=1,apex=record['entry_deletion_mark']),
                         dict(kind='prepared',family=record['family'],graph=prepared,H=prepared_h,q=1,apex=record['entry_deletion_mark']),
                         dict(kind='first_C2_birth',family=record['family'],graph=first_birth,H=prepared_h+weight+2,q=1,apex=record['entry_deletion_mark'])])
    payload='\n'.join(' '.join(map(str,[i,r['n']]+r['rows'])) for i,r in enumerate(words+weak_controls))+'\n'
    with tempfile.TemporaryDirectory(prefix='one_apex_complete_') as tmp:
        exe=Path(tmp)/'counter.exe';subprocess.run(['g++','-O3','-std=c++20',str(COUNTER),'-o',str(exe)],check=True)
        values=[json.loads(line) for line in subprocess.check_output([str(exe),'h'],input=payload,text=True).splitlines()]
        need(len(values)==len(words)+len(weak_controls),'complete independent DP coverage')
        for i,(actual,expected) in enumerate(zip(values,words+weak_controls)):
            need(actual==dict(id=i,H=expected['H']),'literal DP differs from split formula')
        payload='\n'.join(' '.join(map(str,[i,len(r['graph']['rows'])]+r['graph']['rows'])) for i,r in enumerate(controls))+'\n'
        full=[json.loads(line) for line in subprocess.check_output([str(exe)],input=payload,text=True).splitlines()]
    for c,p in zip(controls,full):
        names=c['graph']['names'];need(p['H']==c['H'],'actual source/child count')
        need(p['E2'][names.index(99)]==p['H'],'actual zero-end root')
        need(p['deletion_H'][names.index(c['apex'])]==c['q'],'actual deletion resource')
        c['profile']=p
    need([d['required_compensation'] for d in debts]==[30,66],'exact requested debts')
    result=dict(status='PASS complete one-apex core exclusion',targets=sorted(targets),core_orders=list(range(3,14)),
                strong_words=2047,typed_root_words=typed,literal_weak_boundary_controls=len(weak_controls),
                all_word_counts_verified_by_literal_subset_DP=True,partition_count_controls=True,
                elementary_cutoff=dict(first_excluded_order=14,lower_bound_H=2050),census=census,
                target_hits=[],debt_sources=debts,actual_source_controls=controls,
                consequence='No tournament of any order with H1409 or H1427 has any vertex w with H(T-w)=1.',
                scope='Exact core reduction plus finite exhaustive certificate; targets remain attained with deletion count3.',
                source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD,COUNTER,Path(__file__)]})
    out=HERE/('one_apex_optimized.json' if not __debug__ else 'one_apex.json');out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('actual_source_controls','source_sha256','debt_sources')},indent=2))
    print(json.dumps({'debt_sources':[{k:v for k,v in d.items() if k not in ('rows','names')} for d in debts]},indent=2))


if __name__=='__main__':main()
