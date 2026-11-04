"""Independent literal replay of the saved initial-C3 q3 language."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from check_q3_review import need, phase, star, strong

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
RESULTS=PACKET.parent
SOURCE=PACKET/'star/q3/initial'


def main():
    census=json.loads((SOURCE/'census.json').read_text())
    receipt=json.loads((SOURCE/'receipt.json').read_text())
    witnesses=json.loads((SOURCE/'witnesses.json').read_text())['records']
    originals=json.loads((RESULTS/'c4_port_returns_20261103/entry_selector/receipt.json').read_text())['records'][:2]
    children=[r['stages'][-1]['graph'] for r in originals]
    names=census['names'];apex=names.index(0);keep=[i for i in range(len(names)) if i!=apex]
    cap=[names.index(v) for v in (201,400,401)];prefix=[i for i in keep if i not in cap]
    need(len(names)==13 and all(g['names']==names for g in children),'old physical-name universe')
    rows=list(children[0]['rows']);flips=[]
    for u in prefix:
        for v in cap:
            need(rows[u]>>v&1,'old prefix precedes cap')
            rows[u]^=1<<v;rows[v]^=1<<u;flips.append([names[u],names[v]])
    need(rows==census['initial_rows'] and flips==census['nonstar_relocation_pairs'],'all27 actual nonstar reversals')
    need(census['star_bit_order']==[names[i] for i in keep],'explicit physical mask convention')
    records=census['records'];by_mask={r['mask']:r for r in records}
    need(len(records)==len(by_mask)==4096 and set(by_mask)==set(range(4096)),'complete unique mask universe')
    graphs=[];hits=[];strong_count=0
    for mask in range(4096):
        graph=star(rows,apex,keep,mask);r=by_mask[mask];isstrong=strong(graph)
        need(isstrong==r['strong'],'independent directed strongness');strong_count+=isstrong
        graphs.append(graph)
        if isstrong and r['H'] in (1409,1427):hits.append((mask,graph))
    need(strong_count==1792 and len(hits)==9,'strong and target counts')
    deleted=[sum((rows[u]>>v&1)<<j for j,v in enumerate(keep)) for u in keep]
    with tempfile.TemporaryDirectory(prefix='initial_q3_independent_') as tmp:
        exe=Path(tmp)/'counter.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(HERE/'exact_counter.cpp'),'-o',str(exe)],check=True)
        payload='\n'.join(' '.join(map(str,[i,len(g)]+g)) for i,g in enumerate(graphs))+'\n'
        counts=[json.loads(line) for line in subprocess.check_output([str(exe),'h'],input=payload,text=True).splitlines()]
        need(len(counts)==4096,'all4096 literal counts returned')
        for mask,p in enumerate(counts):need(p==dict(id=mask,H=by_mask[mask]['H']),'literal count versus saved quadratic')
        full_graphs=[g for _,g in hits]+[deleted]+[g['rows'] for g in children]
        payload='\n'.join(' '.join(map(str,[i,len(g)]+g)) for i,g in enumerate(full_graphs))+'\n'
        profiles=[json.loads(line) for line in subprocess.check_output([str(exe)],input=payload,text=True).splitlines()]
        bad=list(graphs[0]);bad[0]|=1
        hostile=subprocess.run([str(exe),'h'],input=' '.join(map(str,[0,len(bad)]+bad))+'\n',text=True,capture_output=True)
        # Uncaught C++ exceptions can have no captured diagnostic on Windows;
        # rejection is the nonzero status with no accepted JSON profile.
        need(hostile.returncode!=0 and not hostile.stdout.strip(),'self-loop hostile rejected')
    need(len(profiles)==12,'nine target, deleted and two retained source profiles')
    for expected,p in zip([census['deleted_profile']]+receipt['original_child_controls'],profiles[-3:]):
        need({k:(expected['id'] if k=='id' else p[k]) for k in expected}==expected,'deleted/source full profile')
    witness_map={w['mask']:w for w in witnesses};need(set(witness_map)=={m for m,_ in hits},'exact saved witness set')
    counts=Counter()
    for (mask,g),p in zip(hits,profiles):
        w=witness_map[mask];r=by_mask[mask]
        need(w['graph']['names']==names and w['graph']['rows']==g,'actual names and rows')
        need({k:(w['profile']['id'] if k=='id' else p[k]) for k in w['profile']}==w['profile'],'all literal witness observers')
        start=[v-p['H'] for v in p['S2']];end=[v-p['H'] for v in p['E2']]
        entries=[dict(name=names[i],sigma=s) for i,s in enumerate(start) if s in (1,2)]
        endzero=[names[i] for i,s in enumerate(end) if s==0]
        need(start==w['start_defects'] and end==w['end_defects'],'all named defects')
        need(entries==w['entries']==r['entries']==[] and endzero==w['end_zero_names']==r['end_zero_names'] and endzero,'exact noentry/endzero boundary')
        need(not w['compatible'] and not r['compatible'],'no target has both guards')
        need(p['deletion_H'][apex]==w['q']==3,'actual q3 deletion resource')
        need({str(v):phase(p['H'],p['E1'][i],p['S2'][i]) for i,v in enumerate(names)}==w['charts'],'fresh all-name phase charts')
        need([names[v] for v in keep if g[v]>>apex&1]==w['incoming'],'actual incoming cut')
        need([names[v] for v in keep if g[apex]>>v&1]==w['outgoing'],'actual outgoing cut')
        counts[p['H']]+=1
    need(counts=={1409:6,1427:3},'complete exact-count fibres')
    pins=[SOURCE/f for f in ('PROOF.md','check_initial.py','census.json','receipt.json','witnesses.json')]+[HERE/'exact_counter.cpp',HERE/'check_q3_review.py',Path(__file__)]
    result=dict(status='PASS independent initial-cap q3 review',all_stars_literal_H_checked=4096,strong_stars=strong_count,
        named_target_full_profiles=9,deleted_and_control_profiles=3,target_counts=dict(counts),target_start_entry_count=0,
        target_with_end_zero_count=9,hostile_self_loop_rejected=True,
        scope='Fixed initial C3->TT9 deleted carrier; no native program or arbitrary-carrier exclusion inferred.',
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    out=HERE/('initial_review_optimized.json' if not __debug__ else 'initial_review.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
