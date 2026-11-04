"""Review the saved fixed-q5 formula and its sole strong target witness."""
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
SOURCE=PACKET/'star/q5'


def main():
    c=json.loads((SOURCE/'census.json').read_text())
    ws=json.loads((SOURCE/'witnesses.json').read_text())['records']
    need(len(ws)==1,'one saved strong target');w=ws[0]
    names=c['names'];a=names.index(0);r=names.index(99);keep=[i for i in range(len(names)) if i!=a]
    need(names==w['graph']['names']==w['source']['names'] and len(names)==13,'all original names retained')
    need(c['star_bit_order']==[names[i] for i in keep],'explicit mask convention')
    need(c['prefix_names']==[99,10,11,2,1,3,4,12] and c['cap_names']==[200,201,400,401] and c['cap_rows']==[8,1,3,6],'declared ordered carrier')
    cap=[names.index(v) for v in c['cap_names']];prefix=[names.index(v) for v in c['prefix_names']]
    rows=c['carrier_rows']
    for i,u in enumerate(prefix):
        need(all(rows[u]>>v&1 for v in prefix[i+1:]+cap),'literal transitive prefix before cap')
    for i,u in enumerate(cap):
        need(sum((rows[u]>>v&1)<<j for j,v in enumerate(cap))==c['cap_rows'][i],'literal cap rows')
    nonstar=[[names[u],names[v]] for u in keep for v in keep if u<v and (rows[u]^w['source']['rows'][u])>>v&1]
    need(nonstar==c['nonstar_flips']==w['nonstar_flips'],'exact nonstar changes')
    target=star(rows,a,keep,w['mask']);need(target==w['graph']['rows'],'actual target matrix from saved cut')
    repaired=list(target);ten=names.index(10);repaired[a]^=1<<ten;repaired[ten]^=1<<a
    deleted=[sum((rows[u]>>v&1)<<j for j,v in enumerate(keep)) for u in keep]
    graphs=[deleted,c['cap_rows'],target,repaired,w['source']['rows']]
    with tempfile.TemporaryDirectory(prefix='q5_witness_review_') as tmp:
        exe=Path(tmp)/'counter.exe';subprocess.run(['g++','-O3','-std=c++20',str(HERE/'exact_counter.cpp'),'-o',str(exe)],check=True)
        payload='\n'.join(' '.join(map(str,[i,len(g)]+g)) for i,g in enumerate(graphs))+'\n'
        profiles=[json.loads(s) for s in subprocess.check_output([str(exe)],input=payload,text=True).splitlines()]
    need(len(profiles)==5,'five literal independent profiles')
    d,cp,p,repair,old=profiles
    need({k:(c['deleted_profile']['id'] if k=='id' else d[k]) for k in c['deleted_profile']}==c['deleted_profile'],'full independent deleted coefficients')
    need(d['H']==cp['H']==5 and cp['P2']==32 and strong(c['cap_rows']),'actual q5 cap')
    records={v['mask']:v for v in c['records']};need(len(records)==4096 and set(records)==set(range(4096)),'complete saved universe')
    hits=[];ns=nr=0
    for mask,record in records.items():
        bits=[mask>>i&1 for i in range(12)]
        value=sum((1-x)*e+x*s for x,e,s in zip(bits,d['E1'],d['S1']))
        value+=sum((1-bits[i])*bits[j]*d['J'][i][j] for i in range(12) for j in range(12))
        need(value==record['H'],'saved count from freshly audited deleted coefficients')
        graph=star(rows,a,keep,mask);isstrong=strong(graph)
        need(isstrong==record['strong'],'independent saved strongness')
        if isstrong:
            ns+=1;gap=5*int(graph[a]>>ten&1);need(gap==record['root_end_defect'],'root guard from fixed source SCC')
            nr+=gap==0
            if value in (1409,1427):hits.append((mask,value))
    need(ns==1920 and nr==960 and hits==[(3021,1409)],'complete fixed-carrier target fibre')
    need({k:(w['profile']['id'] if k=='id' else p[k]) for k in w['profile']}==w['profile'],'literal full target profile')
    sigma=[v-p['H'] for v in p['S2']];end=[v-p['H'] for v in p['E2']]
    need(p['H']==1409 and p['deletion_H'][a]==5 and end[r]==5,'actual target H/q/root gap')
    need(sigma==w['start_defects'] and end==w['end_defects'] and all(v not in (1,2) for v in sigma) and 0 not in end,'every named port failure')
    need(w['entries']==w['end_zero_names']==[] and not w['compatible'],'saved rejection fields')
    need({str(v):phase(p['H'],p['E1'][i],p['S2'][i]) for i,v in enumerate(names)}==w['charts'],'fresh phase observers')
    flips=[[names[i],names[j]] for i in range(13) for j in range(i+1,13) if (w['source']['rows'][i]^target[i])>>j&1]
    need(flips==w['all_flips'] and [v for v in flips if 0 in v]==w['star_flips'],'all actual source-relative changes')
    need([names[v] for v in keep if target[v]>>a&1]==w['incoming'] and [names[v] for v in keep if target[a]>>v&1]==w['outgoing'],'physical star cuts')
    need(strong(repaired) and repair['H']==975 and repair['E2'][r]==975 and repair['deletion_H'][a]==5,'literal single-edge root repair price')
    need(old['H']==1409 and old['deletion_H'][a]==3 and old['E2'][r]==1409,'unchanged actual old child control')
    pins=[SOURCE/f for f in ('PROOF.md','check_q5.py','census.json','witnesses.json')]+[HERE/'exact_counter.cpp',HERE/'check_q3_review.py',Path(__file__)]
    result=dict(status='PASS independent fixed-q5 witness review',literal_full_profiles=5,saved_formula_cuts_replayed=4096,
        strong_states=ns,root_end_zero_states=nr,strong_targets=[dict(mask=m,H=h) for m,h in hits],target_q=5,
        target_root_defect=5,target_end_zero_ports=0,target_start_one_two_ports=0,root_repair_H=975,
        scope='One fixed TT8->C4 carrier; all saved cuts replayed by proved formula with independently counted coefficients, not4096 literal full profiles.',
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    (HERE/('q5_review_optimized.json' if not __debug__ else 'q5_review.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
