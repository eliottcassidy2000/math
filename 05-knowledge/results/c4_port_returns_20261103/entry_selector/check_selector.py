"""Three entry-aware reader selections, actual native repairs, and C2 children."""
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
RESULTS=PACKET.parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


export=load('entry_reader_export',PACKET/'reader'/'export_reader.py')
compiler=load('entry_native_compiler',PACKET/'native'/'reader_repairs'/'check_repairs.py')


def need(ok,message):
    if not ok:
        raise RuntimeError(message)


def phase(h,e,s2):
    need((s2-e)%2==0,'integral phase coordinates')
    x,y=(s2-e)//2,(s2+e)//2
    a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    return dict(H_X_Y=[h,x,y],F49_plus_charge=[a,b,c],
                F49_frobenius_plus_charge=[(a+4*b)%7,(-b)%7,c],
                F343_vector=[h%7,x%7,y%7],
                F343_frobenius=[(h+4*x+6*y)%7,(2*x)%7,(3*x+4*y)%7])


def compile_target(source,target,case):
    rows,names=list(source['rows']),list(source['names'])
    indices={name:i for i,name in enumerate(target['names'])}
    births=[]
    for new in target['fresh_names']:
        incoming=[v for v in names if target['rows'][indices[v]] >> indices[new] & 1]
        rows,names=compiler.add_named(rows,names,new,incoming)
        births.append(dict(name=new,incoming=incoming))
    expected=[sum(1 << j for j,b in enumerate(names)
                  if target['rows'][indices[a]] >> indices[b] & 1) for a in names]
    programs=[]
    for i,a in enumerate(names):
        for j in range(i+1,len(names)):
            if (rows[i]^expected[i]) >> j & 1:
                program=compiler.compile_arc(rows,names,a,names[j],50000+1000*case+2*len(programs))
                rows=compiler.replay(rows,names,program)
                programs.append(program)
    need(rows==expected,'exact native target adjacency')
    need(set(source['names'])<=set(names) and 99 in names,'all old names and root99 retained')
    need(export.reader.strong(rows),'completed actual target strong')
    return dict(permanent_births=births,programs=programs,rows=rows,names=names,
                native_complements=2*len(programs),
                temporary_births=sum(len(p['births']) for p in programs),
                peak_order=max([len(names)]+[len(names)+len(p['births']) for p in programs]))


def main():
    bank=json.loads((export.SELECTOR/'bank.json').read_text())['seeds']
    first=json.loads((PACKET/'reader'/'witnesses.json').read_text())['records']
    tail=json.loads((PACKET/'reader'/'tail.json').read_text())['records']
    qualifying=[]
    for record in first+tail:
        hits=[]
        for hit in record['all_hits']:
            seed=bank[hit['bank_index']]
            word=hit['word']
            sigma=(seed['fields'][1] if word=='0' else
                   seed['fields'][5] if word.endswith('00') else 0)
            if sigma in (1,2):
                value=export.reader.run(seed['fields'],hit['word'])
                need(value[0]+(1 << hit['prefix'])*value[3]==record['target_H'],'selected actual reader count')
                hits.append(dict(**hit,seed=seed,terminal_start_defect=sigma))
        if hits:
            qualifying.append(dict(family=record['family'],t=record['t'],target_H=record['target_H'],hits=hits))
    need(len(qualifying)==13 and sum(len(r['hits']) for r in qualifying)==24,'complete declared entry-aware hit subset')
    records=[]
    choices=(('24',3,0,'11100',2,469,1409,1),('44',3,0,'01100',2,475,1427,1),
             ('24',4,95,'0',3,973,2923,2))
    for case,(family,t,seed_index,word,clock,target_H,child_H,sigma) in enumerate(choices):
        record=next(r for r in first if r['family']==family and r['t']==t)
        hits=[h for h in record['all_hits'] if h['bank_index']==seed_index and h['prefix']==clock and h['word']==word]
        need(len(hits)==1,'chosen alternative is an accepted exact hit')
        hit=hits[0];seed=bank[seed_index]
        target=export.actual_target(seed,hit,record['source'])
        need(len(target['rows'])==11,'alternative target size')
        auxiliary=target['component'][-1]
        need(auxiliary==(201 if case<2 else 200),'stable physical last-reader name')
        native=compile_target(record['source'],target,case)
        rows,names=list(native['rows']),list(native['names'])
        stages=[dict(stage='source',graph=record['source'],expected_H=record['source']['H']),
                dict(stage='prepared_target',graph=dict(rows=rows,names=names,root=99),expected_H=target_H)]
        old_names=list(names)
        first_birth=dict(name=400,incoming=[v for v in names if v!=auxiliary],
                         meaning='new400 beats only the selected last reader')
        rows,names=compiler.add_named(rows,names,400,first_birth['incoming'])
        need(export.reader.strong(rows),'first C2 birth is strong')
        stages.append(dict(stage='C2_first_birth',graph=dict(rows=rows,names=names,root=99)))
        second_birth=dict(name=401,incoming=old_names,
                          meaning='new401 beats new400 and loses to every vertex of the prepared target')
        rows,names=compiler.add_named(rows,names,401,second_birth['incoming'])
        need(export.reader.strong(rows),'C2 child strong')
        need(set(record['source']['names'])<=set(names),'C2 retains all original source names')
        stages.append(dict(stage='C2_child',graph=dict(rows=rows,names=names,root=99),expected_H=child_H))
        def deletion_count(graph,mark):
            keep=[i for i,name in enumerate(graph['names']) if name!=mark]
            induced=[sum(1 << j for j,v in enumerate(keep) if graph['rows'][u] >> v & 1) for u in keep]
            return export.reader.tables(induced)[0][-1]
        entry_mark=target['mark']
        q_before=deletion_count(stages[1]['graph'],entry_mark)
        q_after=deletion_count(stages[-1]['graph'],entry_mark)
        need(q_before==seed['fields'][5] and q_after==3*q_before,'exact consumed deletion resource')
        neutral=[]
        if case<2:
            for a,b in ((auxiliary,401),(400,auxiliary),(401,400)):
                changed=compiler.flip(rows,names,a,b)
                graph=dict(rows=changed,names=list(names),root=99)
                q_repaired=deletion_count(graph,entry_mark)
                need(q_repaired==1,'triangle reversal restores unit deletion count')
                neutral.append(dict(reversed_pair=[a,b],graph=graph,
                    strong=export.reader.strong(changed),restored_deletion_count=q_repaired))
        records.append(dict(family=family,t=t,selected=dict(**hit,seed=seed),
            source=record['source'],prepared_target=target,auxiliary=auxiliary,
            required_start_defect=sigma,entry_deletion_mark=entry_mark,
            deletion_resource=dict(before=q_before,after_C2=q_after),neutral_triangle_tests=neutral,
            native_repair=native,C2=dict(type='EXACT_DIRECT_PERMANENT_BIRTHS',
                births=[first_birth,second_birth],old_arcs_unchanged=True,
                source_H=target_H,target_H=child_H,transported_zero_auxiliary=400),
            stages=stages))
    meta=[];lines=[]
    for record in records:
        for stage in record['stages']:
            graph=stage['graph'];identifier=len(meta)
            lines.append(' '.join(map(str,[identifier,len(graph['rows'])]+graph['rows'])))
            meta.append((record,stage))
    with tempfile.TemporaryDirectory(prefix='entry-selector-profile-') as directory:
        exe=Path(directory)/'profile.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(PACKET/'reader'/'endpoint_profile.cpp'),'-o',str(exe)],check=True)
        result=subprocess.run([str(exe)],input='\n'.join(lines)+'\n',capture_output=True,text=True,check=True)
    need(not result.stderr,'unexpected profile stderr')
    profiles=[json.loads(line) for line in result.stdout.splitlines()]
    need(len(profiles)==12,'twelve fresh full profiles')
    observations=[]
    for profile,(record,stage) in zip(profiles,meta):
        graph=stage['graph'];names=graph['names'];h=profile['H']
        if 'expected_H' in stage:
            need(h==stage['expected_H'],'actual stage Hamilton count')
        start={str(name):profile['S2'][i]-h for i,name in enumerate(names)}
        end={str(name):profile['E2'][i]-h for i,name in enumerate(names)}
        if stage['stage']=='prepared_target':
            need(start[str(record['auxiliary'])]==record['required_start_defect'],'fresh regenerated positive entry')
            need(end['99']==0,'prepared end-ready root99')
        elif stage['stage']=='C2_child':
            need(start['400']==0,'actual transported zero-start auxiliary')
            need(all(value not in (1,2) for value in start.values()),'C2 child does not silently regenerate positive entry')
            need(end['99']==0,'C2 child retains zero-end wrapper root99')
        observations.append(dict(family=record['family'],t=record['t'],stage=stage['stage'],names=names,
            profile=profile,start_defects=start,end_defects=end,
            charts={str(name):phase(h,profile['E1'][i],profile['S2'][i]) for i,name in enumerate(names)},
            scope='Fresh literal actual graph data, full ordered J matrix and passive named phase charts.'))
    neutral_input=[];neutral_meta=[]
    for record in records:
        for test in record['neutral_triangle_tests']:
            graph=test['graph'];identifier=len(neutral_meta)
            neutral_input.append(' '.join(map(str,[identifier,len(graph['rows'])]+graph['rows'])))
            neutral_meta.append((record,test))
    with tempfile.TemporaryDirectory(prefix='entry-neutral-profile-') as directory:
        exe=Path(directory)/'profile.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(PACKET/'reader'/'endpoint_profile.cpp'),'-o',str(exe)],check=True)
        result=subprocess.run([str(exe)],input='\n'.join(neutral_input)+'\n',capture_output=True,text=True,check=True)
    neutral_profiles=[json.loads(line) for line in result.stdout.splitlines()]
    need(len(neutral_profiles)==6,'exact six-test hostile universe')
    for profile,(record,test) in zip(neutral_profiles,neutral_meta):
        root=test['graph']['names'].index(99)
        old_h=record['C2']['target_H']
        test.update(H=profile['H'],root_E2=profile['E2'][root],
                    H_delta=profile['H']-old_h,root_E2_delta=profile['E2'][root]-old_h)
        test['neutral']=test['H_delta']==0 and test['root_E2_delta']==0
        need(not test['neutral'],'no unrecorded neutral one-arc repair')
    need(sum(r['native_repair']['native_complements'] for r in records)>0,'actual native paths recorded')
    paths=[export.SELECTOR/'bank.json',PACKET/'reader'/'witnesses.json',PACKET/'reader'/'tail.json',
           PACKET/'reader'/'export_reader.py',PACKET/'reader'/'endpoint_profile.cpp',
           PACKET/'native'/'reader_repairs'/'check_repairs.py',PACKET/'native'/'reader_repairs'/'ENTRY_SCOPE.md']
    pins={str(path.relative_to(RESULTS)).replace('\\','/'):
          hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for path in paths}
    report=dict(status='FINITE_EXACT_ENTRY_AWARE_SELECTION_AND_THREE_REGENERATED_C2_RETURNS',
        qualifying_target_count=len(qualifying),qualifying_hit_count=sum(len(r['hits']) for r in qualifying),
        qualifying_hits=qualifying,records=records,source_pins=pins,
        scope='All terminal-defect1/2 hits among the already computed32 target list; three actual old-name-preserving native repairs and direct C2 children. No all-t selector or repeated two-digit closure.')
    (HERE/'receipt.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    (HERE/'profiles.json').write_text(json.dumps(dict(status='FRESH_ALL_PORT_PROFILES',records=observations),indent=2)+'\n',encoding='utf-8',newline='\n')
    for record in records:
        native=record['native_repair']
        print(record['family'],'source',record['source']['H'],'prepared',record['C2']['source_H'],
              'child',record['C2']['target_H'],'auxiliary',record['auxiliary'],
              'native complements',native['native_complements'],'peak',native['peak_order'])
    print('PASS:24 entry-aware hits across13 of32 targets;12 fresh full profiles; three regenerated entries; six nonneutral triangle repairs')


if __name__=='__main__':
    main()
