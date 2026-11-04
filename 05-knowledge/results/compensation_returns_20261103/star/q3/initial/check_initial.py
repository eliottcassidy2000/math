"""Last fixed q3 carrier: C3 -> TT9; classify target-count port compatibility."""
from pathlib import Path
import importlib.util
import json
import subprocess
import tempfile
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('q3_parent',HERE.parent/'check_q3.py')
q3=importlib.util.module_from_spec(spec);spec.loader.exec_module(q3)
base=q3.base;need=base.need
def dump(name,value,compact=False):
    (HERE/name).write_text(json.dumps(value,indent=None if compact else 2,
        separators=(',',':') if compact else None)+'\n',encoding='utf-8',newline='\n')
def main():
    source_path=base.OLD/'entry_selector'/'receipt.json'
    originals=json.loads(source_path.read_text())['records'][:2]
    children=[next(s['graph'] for s in r['stages'] if s['stage']=='C2_child') for r in originals]
    names=children[0]['names'];a=names.index(0);keep=[i for i in range(len(names)) if i!=a]
    cap=[names.index(v) for v in (201,400,401)];prefix=[v for v in keep if v not in cap]
    initial=list(children[0]['rows']);nonstar=[]
    for u in prefix:
        for v in cap:
            need(initial[u] >> v & 1,'actual TT9 precedes cap')
            initial[u]^=1<<v;initial[v]^=1<<u;nonstar.append([names[u],names[v]])
    profiler=base.OLD/'reader'/'endpoint_profile.cpp'
    with tempfile.TemporaryDirectory(prefix='initial-q3-cap-') as directory:
        exe=Path(directory)/'profile.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(profiler),'-o',str(exe)],check=True)
        def profiles(rows_list):
            inp='\n'.join(' '.join(map(str,[i,len(rows)]+rows)) for i,rows in enumerate(rows_list))+'\n'
            raw=subprocess.run([str(exe)],input=inp,text=True,capture_output=True,check=True)
            need(not raw.stderr,'unexpected initial profile stderr')
            answer=[json.loads(line) for line in raw.stdout.splitlines()]
            need(len(answer)==len(rows_list),'complete initial profiles');return answer
        dp=profiles([q3.induce(initial,keep)])[0];need(dp['H']==3,'initial cap q3')
        census=[];hits=[]
        for mask in range(1<<len(keep)):
            rows=base.star(initial,a,keep,mask);word=[mask>>j&1 for j in range(len(keep))]
            h=q3.count(dp,word);strong=base.strong(rows)
            item=dict(mask=mask,H=h,strong=strong);census.append(item)
            if strong and h in (1409,1427):hits.append(dict(record=item,rows=rows))
        fresh=profiles([v['rows'] for v in hits])
        # Literal positive controls are the unchanged actual source children.
        controls=profiles([c['rows'] for c in children])
    witnesses=[]
    for hit,p in zip(hits,fresh):
        need(p['H']==hit['record']['H'],'initial quadratic formula versus literal subset DP')
        sigma=[v-p['H'] for v in p['S2']];end=[v-p['H'] for v in p['E2']]
        entries=[dict(name=names[i],sigma=v) for i,v in enumerate(sigma) if v in (1,2)]
        endzero=[names[i] for i,v in enumerate(end) if v==0]
        hit['record'].update(entries=entries,end_zero_names=endzero,compatible=bool(entries and endzero))
        witnesses.append(dict(H=p['H'],mask=hit['record']['mask'],graph=dict(rows=hit['rows'],names=names,inner_mark=0),
            entries=entries,end_zero_names=endzero,compatible=bool(entries and endzero),q=3,
            incoming=[names[v] for v in keep if hit['rows'][v]>>a&1],outgoing=[names[v] for v in keep if hit['rows'][a]>>v&1],
            profile=p,start_defects=sigma,end_defects=end,
            charts={str(v):base.phase(p['H'],p['E1'][i],p['S2'][i]) for i,v in enumerate(names)}))
    summary=[]
    for h in (1409,1427):
        group=[w for w in witnesses if w['H']==h]
        summary.append(dict(H=h,strong_hits=len(group),with_start_entry=sum(bool(w['entries']) for w in group),
            with_end_zero=sum(bool(w['end_zero_names']) for w in group),compatible=sum(w['compatible'] for w in group)))
    dump('census.json',dict(status='FINITE_EXACT_INITIAL_CAP_COMPLETE_STAR',names=names,star_bit_order=[names[v] for v in keep],
        deleted_profile=dp,initial_rows=initial,nonstar_relocation_pairs=nonstar,records=census),True)
    dump('witnesses.json',dict(status='FRESH_FULL_TARGET_PROFILES',records=witnesses))
    dump('receipt.json',dict(status='FINITE_EXACT_INITIAL_CAP_COMPATIBILITY',all_masks=len(census),strong=sum(v['strong'] for v in census),
        summary=summary,full_target_profiles=len(fresh),original_child_controls=controls,
        dependencies={str(source_path.relative_to(base.RESULTS)):base.digest(source_path),
            str(profiler.relative_to(base.RESULTS)):base.digest(profiler),'../check_q3.py':base.digest(HERE.parent/'check_q3.py'),
            '../../check_star.py':base.digest(HERE.parent.parent/'check_star.py')},
        scope='Initial C3->TT9 carrier, any prospective root end-zero vertex and any positive start1/2 entry. No label remapping or native repair is asserted.'))
    print('PASS: initial C3->TT9 all4096 stars;',sum(v['strong'] for v in census),'strong;',json.dumps(summary,separators=(',',':')))
    print('PASS:',len(fresh),'full target profiles and two retained source controls')
if __name__=='__main__':main()
