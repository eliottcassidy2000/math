"""Two exact q3 deleted carriers: original cap and one terminal relocation."""
from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('star_parent',HERE.parent/'check_star.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
need=base.need

def dump(name,value,compact=False):
    (HERE/name).write_text(json.dumps(value,indent=None if compact else 2,
        separators=(',',':') if compact else None)+'\n',encoding='utf-8',newline='\n')

def induce(rows,keep):
    return [sum((rows[u] >> v & 1) << j for j,v in enumerate(keep)) for u in keep]

def count(profile,word):
    return sum(word[j]*profile['S1'][j]+(1-word[j])*profile['E1'][j] for j in range(len(word)))+sum(
        (1-word[i])*word[j]*profile['J'][i][j] for i in range(len(word)) for j in range(len(word)))

def main():
    receipt=base.OLD/'entry_selector'/'receipt.json'
    records=json.loads(receipt.read_text())['records'][:2]
    children=[next(s['graph'] for s in r['stages'] if s['stage']=='C2_child') for r in records]
    names=children[0]['names'];a=names.index(0);rr=names.index(99);z=names.index(200)
    keep=[i for i in range(len(names)) if i!=a]
    cap=[names.index(v) for v in (201,400,401)]
    source=children[0]['rows']
    need(names==children[1]['names'] and all((source[u]^children[1]['rows'][u])&~(1<<a)==0 for u in keep),'common actual q3 deletion carrier')
    moved=list(source)
    for v in cap:
        need(moved[z] >> v & 1,'actual200 precedes cap')
        moved[z]^=1 << v;moved[v]^=1 << z
    carriers=[dict(name='original',rows=source,nonstar_flips=[]),
        dict(name='relocated200',rows=moved,nonstar_flips=[[200,names[v]] for v in cap])]
    profiler=base.OLD/'reader'/'endpoint_profile.cpp'
    with tempfile.TemporaryDirectory(prefix='q3-actual-star-') as directory:
        exe=Path(directory)/'profile.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(profiler),'-o',str(exe)],check=True)
        def profiles(rows_list):
            inp='\n'.join(' '.join(map(str,[i,len(rows)]+rows)) for i,rows in enumerate(rows_list))+'\n'
            run=subprocess.run([str(exe)],input=inp,text=True,capture_output=True,check=True)
            need(not run.stderr,'unexpected q3 profiler stderr')
            ans=[json.loads(line) for line in run.stdout.splitlines()]
            need(len(ans)==len(rows_list),'complete q3 profiles');return ans
        deleted_profiles=profiles([induce(c['rows'],keep) for c in carriers])
        floor=deleted_profiles[1]['J'][keep.index(z)][keep.index(rr)]
        need(floor==(1 << 7)*12==1536,'positive relocated endpoint-joint floor')
        all_records=[];typed=[]
        for ci,(carrier,dp) in enumerate(zip(carriers,deleted_profiles)):
            need(dp['H']==3,'exact retained q3')
            for mask in range(1 << len(keep)):
                word=[mask >> j & 1 for j in range(len(keep))]
                rows=base.star(carrier['rows'],a,keep,mask)
                h=count(dp,word);strong=base.strong(rows)
                if ci==1 and strong:need(h>=floor,'every strong relocated star exceeds the target counts')
                # The second prefix name10 is still the unique possible first
                # vertex in D-root, which has three Hamilton paths.
                rootdef=3*int(rows[a] >> names.index(10)&1) if strong else None
                record=dict(carrier=carrier['name'],mask=mask,H=h,strong=strong,root_end_defect=rootdef)
                all_records.append(record)
                if strong and rootdef==0:
                    typed.append(dict(carrier_index=ci,mask=mask,rows=rows,record=record))
        actual=profiles([v['rows'] for v in typed])
        controls=profiles([c['rows'] for c in children])
    witnesses=[];summaries=[]
    for item,p in zip(typed,actual):
        record=item['record'];need(p['H']==record['H'],'quadratic deleted-joint formula versus full actual DP')
        need(p['E2'][rr]==p['H'],'literal retained-root end guard')
        sigma=[s-p['H'] for s in p['S2']]
        entries=[dict(name=names[i],sigma=s) for i,s in enumerate(sigma) if s in (1,2)]
        record['entries']=entries
        record['sigma200']=sigma[z]
        ci=item['carrier_index']
        if ci==0:need(not entries,'original q3 cap carrier has no positive1/2 entry')
        else:need(sigma[z]==sum(item['rows'][v] >> a & 1 for v in cap),'relocated terminal exact entry law')
        if p['H'] in (1409,1427):
            target=0 if p['H']==1409 else 1
            old=children[target]['rows']
            flips=[[names[u],names[v]] for u in range(len(names)) for v in range(u+1,len(names)) if (old[u]^item['rows'][u]) >> v & 1]
            starflips=[pair for pair in flips if 0 in pair]
            need([pair for pair in flips if 0 not in pair]==sorted([sorted(pair,key=names.index) for pair in carriers[ci]['nonstar_flips']],key=lambda pair:(names.index(pair[0]),names.index(pair[1]))),'exact old nonstar support')
            witnesses.append(dict(carrier=record['carrier'],target_H=p['H'],mask=item['mask'],source=children[target],
                graph=dict(rows=item['rows'],names=names,root=99,inner_mark=0),
                all_flips=flips,star_flips=starflips,star_flip_count=len(starflips),nonstar_flips=carriers[ci]['nonstar_flips'],
                incoming=[names[v] for v in keep if item['rows'][v] >> a & 1],
                outgoing=[names[v] for v in keep if item['rows'][a] >> v & 1],
                q=3,entries=entries,start_defects=sigma,end_defects=[e-p['H'] for e in p['E2']],profile=p,
                charts={str(v):base.phase(p['H'],p['E1'][i],p['S2'][i]) for i,v in enumerate(names)}))
    for ci,c in enumerate(carriers):
        group=[v for v in all_records if v['carrier']==c['name']]
        admitted=[v for v in group if v['strong'] and v['root_end_defect']==0]
        summaries.append(dict(carrier=c['name'],all_masks=len(group),strong=sum(v['strong'] for v in group),typed=len(admitted),
            entry_states=sum(bool(v['entries']) for v in admitted),
            minimum_typed_H=min(v['H'] for v in admitted),
            minimum_positive_entry_H=min((v['H'] for v in admitted if v['entries']),default=None),
            targets=[dict(H=h,total=sum(v['H']==h for v in admitted),
                sigma1=sum(v['H']==h and any(e['sigma']==1 for e in v['entries']) for v in admitted),
                sigma2=sum(v['H']==h and any(e['sigma']==2 for e in v['entries']) for v in admitted)) for h in (1409,1427)]))
    witnesses.sort(key=lambda v:(v['carrier'],v['target_H'],v['star_flip_count'],v['mask']))
    dump('census.json',dict(status='FINITE_EXACT_TWO_DELETED_CARRIERS',names=names,star_bit_order=[names[v] for v in keep],
        deleted_profiles=deleted_profiles,carriers=carriers,records=all_records),True)
    dump('witnesses.json',dict(status='FINITE_EXACT_ACTUAL_NAMED_TARGET_STARS',records=witnesses))
    dump('receipt.json',dict(status='FINITE_EXACT_Q3_RELAXATION',summaries=summaries,full_typed_profiles=len(actual),
        relocated_strong_count_floor=floor,
        retained_controls=controls,dependencies={str(receipt.relative_to(base.RESULTS)):base.digest(receipt),
        str(profiler.relative_to(base.RESULTS)):base.digest(profiler),'../check_star.py':base.digest(HERE.parent/'check_star.py')},
        scope='Both complete fixed-deletion star languages. The relocated carrier changes exactly three nonstar arcs and preservesq3; no q1 restoration or native program claimed.'))
    print('PASS:',json.dumps(summaries,separators=(',',':')))
    print('PASS:',len(actual),'fresh full typed profiles;',len(witnesses),'exact target-count witnesses with all names,joints,phases')

if __name__=='__main__':main()
