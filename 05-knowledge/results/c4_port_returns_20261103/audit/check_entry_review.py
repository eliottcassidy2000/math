"""Bounded independent audit of three named entry/C2 chains; no search."""
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from check_reader_review import need,graph,strong,birth,replay,expected_source,target_rows,phase

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
RESULTS=PACKET.parent
ENTRY=PACKET/'entry_selector'
SCOPE=PACKET/'native/reader_repairs'


def rows_of(g,names):
    return [sum(1<<j for j,b in enumerate(names) if b in g[a]) for a in names]


def induced(g,deleted):
    names=[a for a in g if a!=deleted]
    return dict(rows=rows_of(g,names),names=names)


def flip(g,a,b):
    g=deepcopy(g)
    if b in g[a]:g[a].remove(b);g[b].add(a)
    else:g[b].remove(a);g[a].add(b)
    return g


def c2_replay(record):
    stages=record['stages'];before=graph(stages[1]['graph']['rows'],stages[1]['graph']['names'])
    g=deepcopy(before);aux=record['auxiliary'];first,second=record['C2']['births']
    need(record['C2']['type']=='EXACT_DIRECT_PERMANENT_BIRTHS','typed direct births')
    need(first['name']==400 and set(first['incoming'])==set(g)-{aux},'C2 first sole-outport cut')
    birth(g,first);need(g[400]=={aux},'first outneighbor identity')
    need(strong(g),'first completed C2 graph strong')
    need(g==graph(stages[2]['graph']['rows'],stages[2]['graph']['names']),'first C2 stage equality')
    need(second['name']==401 and set(second['incoming'])==set(before),'C2 second cut')
    birth(g,second);need(g[401]=={400},'second outneighbor identity')
    need(strong(g),'C2 child strong')
    need(g==graph(stages[3]['graph']['rows'],stages[3]['graph']['names']),'child stage equality')
    need(all(g[a]&set(before)==before[a] for a in before),'C2 old arcs unchanged')
    need(set(record['source']['names'])<=set(g) and 99 in g,'child retains original names')
    return before,g


def main():
    receipt=json.loads((ENTRY/'receipt.json').read_text())
    observed_profiles=json.loads((ENTRY/'profiles.json').read_text())['records']
    source_first=json.loads((PACKET/'reader/witnesses.json').read_text())['records']
    source_tail=json.loads((PACKET/'reader/tail.json').read_text())['records']
    bank=json.loads((RESULTS/'recursive_seed_switching_20261103/selector/bank.json').read_text())['seeds']
    desired={};sigma_counts=Counter()
    for record in source_first+source_tail:
        for hit in record['all_hits']:
            word=hit['word'];seed=bank[hit['bank_index']]
            if not word:continue
            need(word.endswith('0'),'accepted nonempty strong-ending word')
            sigma=seed['fields'][1] if word=='0' else seed['fields'][5] if word[-2]=='0' else 0
            if sigma in (1,2):
                key=(record['family'],record['t'],hit['bank_index'],hit['prefix'],word)
                desired[key]=sigma;sigma_counts[sigma]+=1
    declared={}
    for group in receipt['qualifying_hits']:
        for hit in group['hits']:
            need(hit['seed']==bank[hit['bank_index']],'actual qualifying seed')
            declared[group['family'],group['t'],hit['bank_index'],hit['prefix'],hit['word']]=hit['terminal_start_defect']
    need(desired==declared and len(desired)==24 and sigma_counts=={1:23,2:1},'complete saved-hit entry filter')
    need(len({key[:2] for key in desired})==13,'thirteen qualifying target cases')
    expected_chains=[(155,469,773,1409),(157,475,771,1427),(323,973,1699,2923)]
    need(len(receipt['records'])==3,'exact three-chain universe')
    inputs=[];meta=[];totals=Counter();deletion_checks=[]
    def enqueue(g,tag):
        identifier=len(meta);inputs.append(' '.join(map(str,[identifier,len(g['rows'])]+g['rows'])))
        meta.append(tag);return identifier
    profile_lookup={(r['family'],r['t'],r['stage']):r for r in observed_profiles}
    scope_controls=json.loads((SCOPE/'entry_scope.json').read_text())['records']
    for i,r in enumerate(receipt['records']):
        key=r['family'],r['t'];native=r['native_repair'];target=r['prepared_target'];selected=r['selected']
        need(graph(r['source']['rows'],r['source']['names'])==expected_source(*key),'specified original source')
        need(target_rows(bank[selected['bank_index']],selected['word'],selected['prefix'])==target['rows'],'actual selected reader target')
        need((key+(selected['bank_index'],selected['prefix'],selected['word'])) in desired,'selected qualifying tuple')
        adapted=dict(source=r['source'],permanent_births=native['permanent_births'],programs=native['programs'],
                     target_rows=native['rows'],target_names=native['names'],declared_target=target,peak_order=native['peak_order'])
        replay(adapted);need(native['native_complements']==2*len(native['programs']),'native count metadata')
        need(native['peak_order']==12,'native peak twelve')
        totals['native_moves']+=native['native_complements'];totals['arc_gadgets']+=len(native['programs'])
        totals['temporary_helpers']+=sum(len(p['births']) for p in native['programs']);totals['permanent_preparation_births']+=len(native['permanent_births'])
        before,child=c2_replay(r);totals['direct_C2_births']+=2
        need(before==graph(native['rows'],native['names']),'prepared C2 source equals native endpoint')
        for stage in r['stages']:
            need(strong(graph(stage['graph']['rows'],stage['graph']['names'])),'strong declared stage')
            enqueue(stage['graph'],('stage',i,r,stage))
        mark=r['entry_deletion_mark'];q=bank[selected['bank_index']]['fields'][5]
        need(mark==target['mark'],'retained deletion mark')
        enqueue(induced(before,mark),('deletion',q,'prepared',i))
        enqueue(induced(child,mark),('deletion',3*q,'child',i))
        need(r['deletion_resource']=={'before':q,'after_C2':3*q},'q transport metadata')
        tests=r['neutral_triangle_tests'];need(len(tests)==(3 if i<2 else 0),'six-test universe')
        for test in tests:
            g=flip(child,*test['reversed_pair'])
            need(g==graph(test['graph']['rows'],test['graph']['names']),'literal triangle reversal')
            need(strong(g)==test['strong'],'hostile strongness')
            enqueue(test['graph'],('hostile',i,r,test))
            enqueue(induced(g,mark),('deletion',1,'hostile',i))
    for control in scope_controls:
        enqueue(dict(rows=control['rows']),('scope',control))
    hostile=deepcopy(receipt['records'][0]);hostile['C2']['births'][0]['incoming'].append(hostile['auxiliary'])
    try:c2_replay(hostile)
    except ValueError as exc:typed_hostile=str(exc)
    else:raise ValueError('Corrupted C2 cut was accepted')
    with tempfile.TemporaryDirectory(prefix='entry_review_') as tmp:
        exe=Path(tmp)/'profile.exe'
        subprocess.run(['g++','-std=c++20','-O3',str(HERE/'reader_profile.cpp'),'-o',str(exe)],check=True)
        output=subprocess.check_output([str(exe)],input='\n'.join(inputs)+'\n',text=True)
    profiles=[json.loads(line) for line in output.splitlines()];need(len(profiles)==len(meta),'all queued profiles')
    chains=[[] for _ in range(3)];hostile_deltas=[[] for _ in range(2)]
    for p,tag in zip(profiles,meta):
        kind=tag[0];h=p['H']
        if kind=='stage':
            _,i,r,stage=tag;g=stage['graph'];names=g['names'];name=stage['stage'];key=r['family'],r['t'],name
            old=profile_lookup[key];compare=dict(p,id=old['profile']['id']);need(compare==old['profile'],'independent full stage profile')
            need(names==old['names'],'full profile names');chains[i].append(h)
            start={str(a):p['S2'][j]-h for j,a in enumerate(names)};end={str(a):p['E2'][j]-h for j,a in enumerate(names)}
            need(start==old['start_defects'] and end==old['end_defects'],'all named stage defects')
            need({str(a):phase(h,p['E1'][j],p['S2'][j]) for j,a in enumerate(names)}==old['charts'],'all named phase charts')
            if name=='prepared_target':
                need(start[str(r['auxiliary'])]==r['required_start_defect'] and end['99']==0,'regenerated positive entry')
            if name=='C2_child':
                need(start['400']==0 and end['99']==0 and not any(v in (1,2) for v in start.values()),'child transported-zero boundary')
                need(h==3*r['C2']['source_H']+2*r['required_start_defect'],'typed C2 digit')
        elif kind=='deletion':
            need(h==tag[1],'literal marked-deletion count');deletion_checks.append(h)
        elif kind=='hostile':
            _,i,r,test=tag;root=test['graph']['names'].index(99);e=p['E2'][root];old=r['C2']['target_H']
            need((h,e,h-old,e-old)==(test['H'],test['root_E2'],test['H_delta'],test['root_E2_delta']),'hostile H/E debts')
            need(h!=old and e!=old and test['neutral'] is False,'nonneutral hostile repair')
            hostile_deltas[i].append(h-old)
        else:
            control=tag[1]
            for k in ('H','E2','S2'):need(p[k]==control[k],'independent terminal-law scope control')
            need([v-h for v in p['E2']]==control['end_defects'] and [v-h for v in p['S2']]==control['start_defects'],'scope polarity controls')
    need([tuple(v) for v in chains]==expected_chains,'exact three chains')
    need(hostile_deltas==[[-30,-940,-636],[-66,-952,-656]],'six specified debts')
    need(totals=={'native_moves':58,'arc_gadgets':29,'temporary_helpers':29,'permanent_preparation_births':5,'direct_C2_births':6},'exact operation budget')
    for path in (ENTRY/'manifest.json',):
        need(path.exists(),'producer stable manifest')
    pins=[ENTRY/'PROOF.md',ENTRY/'check_selector.py',ENTRY/'receipt.json',ENTRY/'profiles.json',SCOPE/'ENTRY_SCOPE.md',SCOPE/'entry_scope.json',
          HERE/'reader_profile.cpp',HERE/'check_reader_review.py',Path(__file__)]
    result=dict(status='PASS independent three-entry-chain review',saved_target_cases=32,qualifying_targets=13,qualifying_tuples=24,
                terminal_defect_counts=dict(sigma_counts),chains=chains,operation_totals=dict(totals),fresh_stage_full_profiles=12,
                exact_triangle_hostiles=6,literal_deletion_checks=len(deletion_checks),terminal_scope_controls=len(scope_controls),
                triangle_H_and_rootE2_deltas=hostile_deltas,corrupted_C2_cut_rejected=typed_hostile,
                scope='Three finite selected entries and children; no unbounded selector, neutral deletion reset, or reusable two-branch closure.',
                source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    out=HERE/('entry_review_optimized.json' if not __debug__ else 'entry_review.json');out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
