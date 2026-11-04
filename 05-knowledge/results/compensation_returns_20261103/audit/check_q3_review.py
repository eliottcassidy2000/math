"""Independent literal replay of the two saved q3 star languages."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
RESULTS=PACKET.parent
SOURCE=PACKET/'star/q3'


def need(ok,why):
    if not ok:raise ValueError(why)


def strong(rows):
    full=(1<<len(rows))-1
    for root in range(len(rows)):
        seen=1<<root;stack=[root]
        while stack:
            fresh=rows[stack.pop()]&~seen;seen|=fresh
            while fresh:
                bit=fresh&-fresh;fresh-=bit;stack.append(bit.bit_length()-1)
        if seen!=full:return False
    return True


def star(rows,apex,keep,mask):
    adjacency=[set(j for j in range(len(rows)) if row>>j&1) for row in rows]
    for bit,v in enumerate(keep):
        adjacency[apex].discard(v);adjacency[v].discard(apex)
        if mask>>bit&1:adjacency[apex].add(v)
        else:adjacency[v].add(apex)
    return [sum(1<<j for j in neighbors) for neighbors in adjacency]


def phase(h,e,s):
    need((s-e)%2==0,'integral observer');x,y=(s-e)//2,(s+e)//2
    a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    return dict(H_X_Y=[h,x,y],F49_plus_charge=[a,b,c],F49_frobenius_plus_charge=[(a+4*b)%7,-b%7,c],
                F343_vector=[h%7,x%7,y%7],F343_frobenius=[(h+4*x+6*y)%7,2*x%7,(3*x+4*y)%7])


def main():
    census=json.loads((SOURCE/'census.json').read_text());receipt=json.loads((SOURCE/'receipt.json').read_text())
    witnesses=json.loads((SOURCE/'witnesses.json').read_text())['records']
    old=json.loads((RESULTS/'c4_port_returns_20261103/entry_selector/receipt.json').read_text())['records'][:2]
    children=[r['stages'][-1]['graph'] for r in old]
    names=census['names'];a=names.index(0);root=names.index(99);keep=[i for i in range(len(names)) if i!=a]
    need(census['star_bit_order']==[names[i] for i in keep],'explicit mask physical order')
    need(len(names)==13 and all(g['names']==names for g in children),'same old-name universe')
    carriers=census['carriers'];need([c['name'] for c in carriers]==['original','relocated200'],'exact two carriers')
    need(carriers[0]['rows']==children[0]['rows'],'original actual child carrier')
    moved=list(children[0]['rows']);z=names.index(200);cap=[names.index(v) for v in (201,400,401)]
    for v in cap:moved[z]^=1<<v;moved[v]^=1<<z
    need(moved==carriers[1]['rows'],'exact three nonstar relocation flips')
    records=census['records'];need(len(records)==8192,'all two times4096 masks')
    bykey={(r['carrier'],r['mask']):r for r in records}
    need(len(bykey)==8192 and set(bykey)=={(c['name'],m) for c in carriers for m in range(4096)},'complete unique saved mask universe')
    typed=[];graphs=[];strong_counts=Counter()
    for c in carriers:
        for mask in range(4096):
            rows=star(c['rows'],a,keep,mask);record=bykey[c['name'],mask]
            isstrong=strong(rows);need(isstrong==record['strong'],'independent directed strongness')
            strong_counts[c['name']]+=isstrong
            if isstrong:
                # root99 has unique possible inport0. The three paths of
                # D-root all start at physical10, hence this exact gap.
                need(rows[a]>>root&1,'strong root must have inport0')
                gap=3*int(rows[a]>>names.index(10)&1)
                need(gap==record['root_end_defect'],'root-guard formula')
                if gap==0:typed.append((c,mask,rows,record))
            graphs.append((rows,record))
    need(len(typed)==1408,'complete typed universe')
    with tempfile.TemporaryDirectory(prefix='q3_independent_') as tmp:
        exe=Path(tmp)/'counter.exe';subprocess.run(['g++','-O3','-std=c++20',str(HERE/'exact_counter.cpp'),'-o',str(exe)],check=True)
        payload='\n'.join(' '.join(map(str,[i,len(rows)]+rows)) for i,(rows,_) in enumerate(graphs))+'\n'
        h_values=[json.loads(s) for s in subprocess.check_output([str(exe),'h'],input=payload,text=True).splitlines()]
        need(len(h_values)==8192,'all literal counts')
        for i,(p,(_,r)) in enumerate(zip(h_values,graphs)):need(p==dict(id=i,H=r['H']),'literal count differs from saved quadratic')
        special=[]
        for c in carriers:
            special.append([sum((c['rows'][u]>>v&1)<<j for j,v in enumerate(keep)) for u in keep])
        special.extend(g['rows'] for g in children)
        payload='\n'.join(' '.join(map(str,[i,len(rows)]+rows)) for i,rows in enumerate([v[2] for v in typed]+special))+'\n'
        profiles=[json.loads(s) for s in subprocess.check_output([str(exe)],input=payload,text=True).splitlines()]
    need(len(profiles)==1412,'full typed and source profiles')
    for index,p in enumerate(profiles[-4:]):
        expected=(census['deleted_profiles']+receipt['retained_controls'])[index]
        need({k:(expected['id'] if k=='id' else p[k]) for k in expected}==expected,'independent deleted/control profile')
    witness_map={(w['carrier'],w['mask']):w for w in witnesses};observed_witnesses=set();entry_counts=Counter();target_counts=Counter()
    for (c,mask,rows,r),p in zip(typed,profiles):
        need(p['deletion_H'][a]==3 and p['E2'][root]==p['H']==r['H'],'actual typed H/q/root response')
        sigmas=[v-p['H'] for v in p['S2']]
        entries=[dict(name=names[i],sigma=s) for i,s in enumerate(sigmas) if s in (1,2)]
        need(entries==r['entries'] and sigmas[z]==r['sigma200'],'all-name next-entry classification')
        entry_counts[c['name']]+=bool(entries)
        if c['name']=='original':need(not entries,'original carrier entry obstruction')
        else:need(sigmas[z]==sum(rows[v]>>a&1 for v in cap),'relocated terminal entry law')
        if p['H'] in (1409,1427):
            key=c['name'],mask;need(key in witness_map,'missing saved target witness');observed_witnesses.add(key);w=witness_map[key]
            need(w['graph']['rows']==rows and w['graph']['names']==names,'actual witness names/rows')
            need({k:(w['profile']['id'] if k=='id' else p[k]) for k in w['profile']}==w['profile'],'full witness profile')
            need(sigmas==w['start_defects'] and [v-p['H'] for v in p['E2']]==w['end_defects'],'witness defects')
            need({str(v):phase(p['H'],p['E1'][i],p['S2'][i]) for i,v in enumerate(names)}==w['charts'],'fresh witness phase charts')
            source=children[0 if p['H']==1409 else 1]
            flips=[[names[i],names[j]] for i in range(len(names)) for j in range(i+1,len(names)) if (source['rows'][i]^rows[i])>>j&1]
            need(flips==w['all_flips'] and [pair for pair in flips if 0 in pair]==w['star_flips'],'exact changed named pairs')
            target_counts[c['name'],p['H']]+=1
    need(observed_witnesses==set(witness_map) and len(observed_witnesses)==9,'exact complete target witness set')
    need(entry_counts=={'original':0,'relocated200':384},'finite entry counts')
    need(target_counts=={('original',1409):6,('original',1427):3},'finite target counts')
    pins=[SOURCE/'PROOF.md',SOURCE/'check_q3.py',SOURCE/'census.json',SOURCE/'receipt.json',SOURCE/'witnesses.json',HERE/'exact_counter.cpp',Path(__file__)]
    result=dict(status='PASS independent literal q3 review',carriers=2,all_stars_literal_H_checked=8192,
                typed_full_profiles=1408,deleted_and_control_profiles=4,all_named_target_witnesses=9,
                strong_counts=dict(strong_counts),entry_counts=dict(entry_counts),target_counts=[dict(carrier=c,H=h,count=n) for (c,h),n in sorted(target_counts.items())],
                scope='Two fixed deleted carriers only; no new native programs, arbitrary-order exclusion, or returning-entry closure.',
                source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    out=HERE/('q3_review_optimized.json' if not __debug__ else 'q3_review.json');out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
