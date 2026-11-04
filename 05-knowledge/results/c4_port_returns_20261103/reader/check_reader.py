"""Independent interval inverse and fresh actual-graph sidecar audit."""
import json
from pathlib import Path
import subprocess
import tempfile
from export_reader import HERE,SELECTOR,reader,need,all_hits


def phase(h,e,s2):
    need((s2-e)%2==0,'integral golden endpoint coordinates')
    x,y=(s2-e)//2,(s2+e)//2
    a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    return dict(H_X_Y=[h,x,y],F49_plus_charge=[a,b,c],
                F49_frobenius_plus_charge=[(a+4*b)%7,(-b)%7,c],
                F343_vector=[h%7,x%7,y%7],
                F343_frobenius=[(h+4*x+6*y)%7,(2*x)%7,(3*x+4*y)%7])


def main():
    bank=json.loads((SELECTOR/'bank.json').read_text())['seeds']
    data=json.loads((HERE/'witnesses.json').read_text())
    interval_integers=0
    inputs=[]
    meta=[]
    for record in data['records']:
        target=record['target_H']
        exact=all_hits(bank,target)
        declared={(h['bank_index'],h['prefix'],h['word']) for h in record['all_hits']}
        need({(h['bank_index'],h['prefix'],h['word']) for h in exact}==declared,'section inverse replay')
        independent=set()
        for index,seed in enumerate(bank):
            hits,tried=reader.inverse(seed['fields'],target)
            interval_integers+=tried
            independent.update((index,h['prefix'],h['word']) for h in hits)
        need(independent==declared,'independent complete interval inverse')
        for role in ('source','target'):
            graph=record[role]
            need(reader.strong(graph['rows']),'actual strong graph')
            identifier=len(meta)
            inputs.append(' '.join(map(str,[identifier,len(graph['rows'])]+graph['rows'])))
            meta.append((record,role,graph))
    with tempfile.TemporaryDirectory(prefix='reader-fresh-profiles-') as directory:
        exe=Path(directory)/'profile.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(HERE/'endpoint_profile.cpp'),'-o',str(exe)],check=True)
        result=subprocess.run([str(exe)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    need(not result.stderr,'unexpected endpoint profiler stderr')
    profiles=[json.loads(line) for line in result.stdout.splitlines()]
    need(len(profiles)==40,'two actual profiles per target')
    sidecars=[]
    entry_records=[]
    for value,(record,role,graph) in zip(profiles,meta):
        need(value['id']==len(sidecars),'profile order')
        expected=record['target_H'] if role=='target' else record['source']['H']
        need(value['H']==expected,'fresh actual Hamilton count')
        root=graph['names'].index(99)
        defect=value['E2'][root]-value['H']
        wanted=0 if role=='target' else (16 if record['family']=='24' else 12)
        need(defect==wanted,'fresh outer-root defect')
        charts={str(name):phase(value['H'],value['E1'][i],value['S2'][i])
                for i,name in enumerate(graph['names'])}
        sidecars.append(dict(family=record['family'],t=record['t'],role=role,
            names=graph['names'],old_names=record['source']['names'],outer_root=99,
            outer_root_defect=defect,profile=value,charts=charts,
            scope='Fresh observations on actual graphs; name roles are passive and no field-Frobenius graph operation is inferred.'))
        if role=='target':
            start_debts={str(name):value['S2'][i]-value['H'] for i,name in enumerate(graph['names'])}
            end_debts={str(name):value['E2'][i]-value['H'] for i,name in enumerate(graph['names'])}
            start_one_two=[int(name) for name,debt in start_debts.items() if debt in (1,2)]
            end_one_two=[int(name) for name,debt in end_debts.items() if debt in (1,2)]
            need(not start_one_two and not end_one_two,'all-port1/2 entry obstruction')
            need([int(name) for name,debt in end_debts.items() if debt==0]==[99],'unique zero-end gauge')
            entry_records.append(dict(family=record['family'],t=record['t'],H=value['H'],
                names=graph['names'],actual_prefix=graph['prefix'],component=graph['component'],
                old_names=record['source']['names'],old_prefix=record['source']['prefix'],
                promoted_prefix_names=graph['promoted_prefix_names'],fresh_names=graph['fresh_names'],
                start_defects=start_debts,end_defects=end_debts,
                start_zero=[int(name) for name,debt in start_debts.items() if debt==0],
                start_one_two=start_one_two,end_one_two=end_one_two,
                minimum_start_defect=min(start_debts.values())))
    (HERE/'sidecars.json').write_text(json.dumps(dict(status='FRESH_NAMED_GRAPH_SIDECARS',
        charts_definition='Inherited golden/F49-plus-charge/F343 decoder; every ordered J pair retained, including promoted old prefix names.',
        records=sidecars),indent=2)+'\n',encoding='utf-8',newline='\n')
    (HERE/'ENTRY_DEFECTS.json').write_text(json.dumps(dict(
        status='FINITE_EXACT_ALL_ACTUAL_TARGET_VERTICES',records=entry_records,
        consequence='No selected target or its converse has a vertex with start defect1 or2. Remapping vertices alone cannot enable an immediate whole-graph C2 return3H+2/4.',
        boundary='Outerroot99 is end-ready with defect0; this does not supply start defect1/2. Some other vertices supply the zero channel. No arbitrary-operation impossibility asserted.'),indent=2)+'\n',encoding='utf-8',newline='\n')
    receipt=dict(status='INDEPENDENT_FINITE_EXACT_READER_FAMILY_AUDIT',
        targets=20,bank_seeds=219,independent_interval_integers=interval_integers,
        complete_section_and_interval_hit_sets_agree=True,
        fresh_actual_profiles=40,source_and_target_orders_at_most=17,
        all_target_vertex_entry_audit=dict(targets=20,start_defect1_or2=0,end_defect1_or2=0,
                                           unique_end_zero_name=99),
        physical_name_embedding=True,native_transport_asserted=False,
        all_order_acceptance_asserted=False,
        scope='Exact finite acceptance, minimum bank sizes, full actual graph counts and named observers; no extension above the separately bounded t16 tail.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: all20 complete hit sets match independent interval inverse over',interval_integers,'candidate integers')
    print('PASS: 40 fresh actual-graph profiles, full named ordered J matrices, and regenerated endpoint/phase charts')


if __name__=='__main__':
    main()
