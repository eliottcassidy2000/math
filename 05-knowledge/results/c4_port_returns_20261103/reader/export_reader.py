"""Focused exact reader acceptance and canonical retained-name targets."""
from pathlib import Path
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
SELECTOR=RESULTS/'recursive_seed_switching_20261103'/'selector'
READER_PATH=RESULTS/'zero_ready_component_returns_20261103'/'inverse'/'reader.py'


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


iv=load('section_inverse',SELECTOR/'inverse.py')
reader=load('literal_retained_reader',READER_PATH)


def need(ok,message):
    if not ok:
        raise RuntimeError(message)


def source(family,t):
    B=[30,28,8,16,4] if family=='24' else [30,28,24,16,0]
    W={0,1,2} if family=='24' else {0,2,4}
    prefix=list(range(10,10+t))
    names=prefix+list(range(5))+[99]
    rows=[0]*len(names)
    for i in range(t):
        rows[i]=sum(1 << j for j in range(i+1,t+5))
        rows[-1]|=1 << i
    for i in range(5):
        rows[t+i]=B[i] << t
        if i in W:
            rows[t+i]|=1 << (len(names)-1)
        else:
            rows[-1]|=1 << (t+i)
    need(reader.strong(rows),'actual source strongness')
    H=(-13 if family=='24' else -11)+21*(1 << t)
    return dict(rows=rows,names=names,root=99,prefix=prefix,
                suffix=list(range(5)),ports=sorted(W),H=H)


def actual_target(seed,hit,old):
    component=reader.graph(tuple(seed['rows']),seed['root'],hit['word'])
    component_order=len(component)
    clock=hit['prefix']
    rows=list(component)
    for j in range(clock):
        rows.append((1 << component_order)-1 |
                    sum(1 << (component_order+k) for k in range(j+1,clock)))
    rows=list(reader.birth(tuple(rows),((1 << len(rows))-1) ^ (1 << seed['root'])))
    need(reader.strong(rows),'actual target strongness')
    need(len(rows)>=len(old['names']),'no old-name deletion capacity')
    need(component_order>=5,'five old suffix names fit in component')
    names=[None]*len(rows)
    names[-1]=99
    for i in range(5):
        names[i]=i
    kept_prefix=min(clock,len(old['prefix']))
    for i in range(kept_prefix):
        names[component_order+i]=old['prefix'][i]
    promoted=old['prefix'][kept_prefix:]
    free_component=[i for i in range(component_order) if names[i] is None]
    need(len(free_component)>=len(promoted),'promoted old prefix capacity')
    for index,name in zip(free_component,promoted):
        names[index]=name
    fresh=200
    for i in range(len(names)):
        if names[i] is None:
            while fresh in old['names'] or fresh in names:
                fresh+=1
            names[i]=fresh
            fresh+=1
    need(len(names)==len(set(names)),'injective target names')
    need(set(old['names'])<=set(names),'every old physical name retained')
    return dict(rows=rows,names=names,root=99,
                mark=names[seed['root']],component=names[:component_order],
                prefix=names[component_order:-1],promoted_prefix_names=promoted,
                fresh_names=[n for n in names if n not in old['names']],
                scope='Actual named target with old names embedded; no native conversion or common-frame claim.')


def all_hits(bank,target):
    hits=[]
    for index,seed in enumerate(bank):
        found,_=iv.inverse(seed['fields'],target)
        for hit in found:
            fields=reader.run(tuple(seed['fields']),hit['word'])
            need(fields[0]+(1 << hit['prefix'])*fields[3]==target,'independent reader replay')
            hits.append(dict(bank_index=index,order=len(seed['rows'])+len(hit['word'])+hit['prefix']+1,
                             fields=fields,**hit))
    hits.sort(key=lambda h:(h['order'],h['bank_index'],h['prefix'],h['word']))
    return hits


def main():
    bank=json.loads((SELECTOR/'bank.json').read_text())['seeds']
    need(len(bank)==219,'declared bank')
    expected={'24':[9,8,10,11,4,6,2,5,9,6], '44':[6,4,3,7,5,7,5,3,3,3]}
    records=[]
    for family,constant in (('24',35),('44',29)):
        counts=[]
        for t in range(1,11):
            target=63*(1 << t)-constant
            hits=all_hits(bank,target)
            need(hits,'requested finite target accepted')
            best=hits[0]
            seed=bank[best['bank_index']]
            old=source(family,t)
            need(3*old['H']+4==target,'actual source-to-requested-count relation')
            actual=actual_target(seed,best,old)
            need(len(actual['rows'])==best['order'],'minimum bank size')
            records.append(dict(family=family,t=t,target_H=target,hit_count=len(hits),
                minimum_order_in_fixed_bank=best['order'],all_hits=hits,
                selected=dict(**best,seed=seed),source=old,target=actual))
            counts.append(len(hits))
        need(counts==expected[family],'declared finite acceptance counts')
    common=[]
    for family in ('24','44'):
        family_records=[r for r in records if r['family']==family]
        for start in (1,2,3,4):
            shared=None
            for record in family_records[start-1:]:
                keys={(h['bank_index'],h['prefix']) for h in record['all_hits']}
                shared=keys if shared is None else shared & keys
            need(not shared,'recorded finite common-key observation')
            common.append(dict(family=family,start_t=start,end_t=10,shared_seed_clock_keys=sorted(shared)))
    pins={str(p.relative_to(RESULTS)).replace('\\','/'):
          hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
          for p in (SELECTOR/'bank.json',SELECTOR/'inverse.py',READER_PATH)}
    report=dict(status='FINITE_EXACT_TWENTY_AFFINE_READER_TARGETS',
        bank_scope='The inherited219 actual eight-field seed representatives; all word lengths and all nonempty reader prefix clocks.',
        naming='Old outerroot99, suffixnames0..4, and oldprefix10..t+9 survive injectively; old prefix names may move into the component.',
        minimum_scope='Minimum vertex count among all accepted words/clocks in this fixed219-key bank; not a global tournament or native-cost minimum.',
        closure_scope='Twenty successes do not prove all-t acceptance, essential infinite switching, native transport, or regenerated entry.',
        records=records,finite_common_key_observations=common,source_pins=pins,
        reproduction='python 05-knowledge/results/c4_port_returns_20261103/reader/export_reader.py')
    (HERE/'witnesses.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    for family in ('24','44'):
        rs=[r for r in records if r['family']==family]
        print(family,'hits',[r['hit_count'] for r in rs],'minimum orders',[r['minimum_order_in_fixed_bank'] for r in rs])
    print('PASS: all20 strong actual targets retain every source name, including outerroot99')


if __name__=='__main__':
    main()
