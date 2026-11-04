"""One bounded q5 carrier: TT8 -> strong C4, with all actual names retained."""
from pathlib import Path
from itertools import permutations
import importlib.util
import json
import subprocess
import tempfile
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('q3_formula',HERE.parent/'q3'/'check_q3.py')
q3=importlib.util.module_from_spec(spec);spec.loader.exec_module(q3)
base=q3.base;need=base.need
def dump(name,value,compact=False):
    (HERE/name).write_text(json.dumps(value,indent=None if compact else 2,
        separators=(',',':') if compact else None)+'\n',encoding='utf-8',newline='\n')
def main():
    source_path=base.OLD/'entry_selector'/'receipt.json'
    originals=json.loads(source_path.read_text())['records'][:2]
    children=[next(s['graph'] for s in r['stages'] if s['stage']=='C2_child') for r in originals]
    names=children[0]['names'];a=names.index(0);root=names.index(99);keep=[i for i in range(len(names)) if i!=a]
    cap_names=[200,201,400,401];prefix_names=[99,10,11,2,1,3,4,12]
    cap=[names.index(v) for v in cap_names];prefix=[names.index(v) for v in prefix_names];cap_rows=[8,1,3,6]
    carrier=list(children[0]['rows'])
    for u in keep:carrier[u]&=1<<a
    for i,u in enumerate(prefix):
        for v in prefix[i+1:]+cap:carrier[u]|=1<<v
    for i,u in enumerate(cap):
        for j,v in enumerate(cap):carrier[u]|=(cap_rows[i]>>j&1)<<v
    nonstar=[[names[u],names[v]] for u in keep for v in keep if u<v and (carrier[u]^children[0]['rows'][u])>>v&1]
    literal_cap_paths=[list(p) for p in permutations(range(4)) if all(cap_rows[p[i]]>>p[i+1]&1 for i in range(3))]
    need(len(literal_cap_paths)==5 and base.strong(cap_rows),'independent cap permutation paths and strongness')
    profiler=base.OLD/'reader'/'endpoint_profile.cpp'
    with tempfile.TemporaryDirectory(prefix='q5-terminal-cap-') as directory:
        exe=Path(directory)/'profile.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(profiler),'-o',str(exe)],check=True)
        def profiles(rows_list):
            inp='\n'.join(' '.join(map(str,[i,len(rows)]+rows)) for i,rows in enumerate(rows_list))+'\n'
            raw=subprocess.run([str(exe)],input=inp,text=True,capture_output=True,check=True)
            need(not raw.stderr,'unexpected q5 profile stderr')
            answer=[json.loads(line) for line in raw.stdout.splitlines()]
            need(len(answer)==len(rows_list),'complete q5 literal profiles');return answer
        dp,cp=profiles([q3.induce(carrier,keep),cap_rows])
        need(dp['H']==5 and cp['H']==5 and cp['P2']==32,'actual q5 deleted carrier and cap')
        census=[];hits=[]
        for mask in range(1<<len(keep)):
            rows=base.star(carrier,a,keep,mask);word=[mask>>j&1 for j in range(len(keep))]
            h=q3.count(dp,word);strong=base.strong(rows)
            expected_strong=bool(rows[a]>>root&1 and any(rows[v]>>a&1 for v in cap))
            need(strong==expected_strong,'exact strongness guard versus reachability')
            rootdef=5*(rows[a]>>names.index(10)&1) if strong else None
            item=dict(mask=mask,H=h,strong=strong,root_end_defect=rootdef);census.append(item)
            if strong and h in (1409,1427):hits.append(dict(record=item,rows=rows))
        target_profiles=profiles([v['rows'] for v in hits])
        # Actual-source positive controls and deterministic weak/strong/guard
        # controls use the independent literal profiler, not just the formula.
        control_masks=sorted(set([0,4095,1<<keep.index(root)]+[1<<i for i in range(12)]+[4095^(1<<i) for i in range(12)]))
        control_rows=[base.star(carrier,a,keep,mask) for mask in control_masks]
        controls=profiles(control_rows+[c['rows'] for c in children])
    witnesses=[]
    for hit,p in zip(hits,target_profiles):
        need(p['H']==hit['record']['H'],'q5 quadratic formula versus actual target H')
        sigma=[v-p['H'] for v in p['S2']];end=[v-p['H'] for v in p['E2']]
        need(end[root]==hit['record']['root_end_defect'],'literal target root defect')
        entries=[dict(name=names[i],sigma=v) for i,v in enumerate(sigma) if v in (1,2)]
        compatible=bool(entries and end[root]==0)
        hit['record'].update(entries=entries,compatible=compatible)
        old=children[0 if p['H']==1409 else 1]['rows']
        flips=[[names[u],names[v]] for u in range(len(names)) for v in range(u+1,len(names)) if (old[u]^hit['rows'][u])>>v&1]
        witnesses.append(dict(H=p['H'],mask=hit['record']['mask'],graph=dict(rows=hit['rows'],names=names,root=99,inner_mark=0),
            source=children[0 if p['H']==1409 else 1],entries=entries,root_end_defect=end[root],compatible=compatible,q=5,
            end_zero_names=[names[i] for i,v in enumerate(end) if v==0],
            all_flips=flips,nonstar_flips=nonstar,star_flips=[p for p in flips if 0 in p],
            incoming=[names[v] for v in keep if hit['rows'][v]>>a&1],outgoing=[names[v] for v in keep if hit['rows'][a]>>v&1],
            profile=p,start_defects=sigma,end_defects=end,
            charts={str(v):base.phase(p['H'],p['E1'][i],p['S2'][i]) for i,v in enumerate(names)}))
    for mask,p in zip(control_masks,controls):
        need(p['H']==census[mask]['H'],'literal hostile formula control')
        if census[mask]['strong']:need(p['E2'][root]-p['H']==census[mask]['root_end_defect'],'literal hostile root guard')
    need([p['H'] for p in controls[-2:]]==[1409,1427],'unchanged actual child positive controls')
    summary=[]
    for h in (1409,1427):
        group=[w for w in witnesses if w['H']==h]
        summary.append(dict(H=h,strong_hits=len(group),with_root_end_zero=sum(w['root_end_defect']==0 for w in group),
            with_start_entry=sum(bool(w['entries']) for w in group),compatible=sum(w['compatible'] for w in group)))
    need(len(witnesses)==1 and witnesses[0]['mask']==3021 and not witnesses[0]['end_zero_names'],'unique strong1409 hostile has no end-zero gauge')
    repaired_mask=3021^(1 << keep.index(names.index(10)))
    need(census[repaired_mask]['H']==975 and census[repaired_mask]['root_end_defect']==0,'single edge root-guard repair loses count')
    dump('census.json',dict(status='FINITE_EXACT_TT8_C4_COMPLETE_STAR',names=names,star_bit_order=[names[v] for v in keep],
        deleted_profile=dp,carrier_rows=carrier,cap_names=cap_names,cap_rows=cap_rows,prefix_names=prefix_names,
        nonstar_flips=nonstar,records=census),True)
    dump('witnesses.json',dict(status='FRESH_FULL_Q5_TARGET_PROFILES',records=witnesses))
    dump('controls.json',dict(status='INDEPENDENT_LITERAL_Q5_CONTROLS',cap_profile=cp,cap_literal_paths=literal_cap_paths,
        control_masks=control_masks,profiles=controls))
    dump('receipt.json',dict(status='FINITE_EXACT_SINGLE_Q5_CARRIER',all_masks=len(census),strong=sum(v['strong'] for v in census),
        root_end_zero=sum(v['strong'] and v['root_end_defect']==0 for v in census),summary=summary,
        unique_count_hostile=dict(mask=3021,H=1409,root_end_defect=5,end_zero_names=[],start_one_two=[],
            flip=[0,10],after_flip_mask=repaired_mask,after_flip_H=975,after_flip_root_end_defect=0),
        full_target_profiles=len(target_profiles),literal_controls=len(controls)+2,
        dependencies={str(source_path.relative_to(base.RESULTS)):base.digest(source_path),
            str(profiler.relative_to(base.RESULTS)):base.digest(profiler),'../q3/check_q3.py':base.digest(HERE.parent/'q3'/'check_q3.py'),
            '../check_star.py':base.digest(HERE.parent/'check_star.py')},
        scope='Exactly TT8->the specified strong C4, all4096 stars. No arbitrary q5 classification, neutral native conversion, or reusable return closure.'))
    print('PASS: TT8->C4 all4096 stars;',sum(v['strong'] for v in census),'strong;',json.dumps(summary,separators=(',',':')))
    print('PASS:',len(target_profiles),'full target profiles;',len(controls)+2,'literal controls including both retained actual children')
if __name__=='__main__':main()
