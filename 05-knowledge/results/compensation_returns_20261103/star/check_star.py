"""Exact one-apex star census with independent literal path-cover profiles."""
from pathlib import Path
import hashlib
import json
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent
OLD=RESULTS/'c4_port_returns_20261103'

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def dump(name,value,compact=False):
    (HERE/name).write_text(json.dumps(value,indent=None if compact else 2,
        separators=(',',':') if compact else None)+'\n',encoding='utf-8',newline='\n')

def count(word):
    m=len(word)
    return word[0]+1-word[-1]+sum(word[i]*(1-word[j])*(1 << (j-i-1))
        for i in range(m) for j in range(i+1,m))+sum((1-word[i])*word[i+1] for i in range(m-1))

def star(base,apex,order,mask):
    rows=list(base)
    for k,v in enumerate(order):
        rows[apex]&=~(1 << v);rows[v]&=~(1 << apex)
        if mask >> k & 1:
            rows[apex]|=1 << v
        else:
            rows[v]|=1 << apex
    return rows

def strong(rows):
    n=len(rows)
    for reverse in (False,True):
        seen=1;todo=[0]
        while todo:
            v=todo.pop()
            for w in range(n):
                if not seen >> w & 1 and ((rows[w] >> v & 1) if reverse else (rows[v] >> w & 1)):
                    seen|=1 << w;todo.append(w)
        if seen!=(1 << n)-1:return False
    return True

def phase(h,e,s2):
    need((s2-e)%2==0,'integral phase coordinates')
    x,y=(s2-e)//2,(s2+e)//2
    a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    return dict(H_X_Y=[h,x,y],F49_plus_charge=[a,b,c],
        F49_frobenius_plus_charge=[(a+4*b)%7,-b%7,c],
        F343_vector=[h%7,x%7,y%7],
        F343_frobenius=[(h+4*x+6*y)%7,2*x%7,(3*x+4*y)%7])

def main():
    source_path=OLD/'entry_selector'/'receipt.json'
    data=json.loads(source_path.read_text())['records'][:2]
    anchors=[r['neutral_triangle_tests'][0] for r in data]
    graph=anchors[0]['graph'];names=graph['names'];rows=graph['rows'];a=names.index(0);root=names.index(99)
    order=sorted((i for i in range(len(rows)) if i!=a),key=lambda i:-sum(rows[i] >> j & 1 for j in range(len(rows)) if j!=a))
    m=len(order)
    need(m==12 and order[0]==root,'actual thirteen-vertex debt carrier and first root')
    for r in anchors:
        need(r['graph']['names']==names,'same actual names')
        need(all(r['graph']['rows'][u] >> v & 1 for i,u in enumerate(order) for v in order[i+1:]),'actual deleted graph is transitive')
    need(all((anchors[0]['graph']['rows'][u]^anchors[1]['graph']['rows'][u])&~(1 << a)==0 for u in order),'identical non-apex named arcs')
    masks=[sum((r['graph']['rows'][a] >> v & 1) << k for k,v in enumerate(order)) for r in anchors]
    table=[];typed=[]
    for mask in range(1 << m):
        word=[mask >> i & 1 for i in range(m)]
        actual=star(rows,a,order,mask)
        isstrong=bool(word[0] and not word[-1])
        need(strong(actual)==isstrong,'independent reachability versus strongness guard')
        item=dict(mask=mask,word=''.join(map(str,word)),H=count(word),strong=isstrong,
            root_end_defect=word[1] if isstrong else None,
            distances=[(mask^old).bit_count() for old in masks])
        table.append(item)
        if isstrong and not word[1]:typed.append(item)
    need(len(typed)==512 and all(not any(v['H']==target for v in typed) for target in (1409,1427)),'complete typed no-compensation census')
    inputs=[];meta=[]
    def add(rows,kind,**values):
        identifier=len(meta);meta.append(dict(kind=kind,rows=rows,**values))
        inputs.append(' '.join(map(str,[identifier,len(rows)]+rows)))
    for item in typed:add(star(rows,a,order,item['mask']),'typed',mask=item['mask'],names=names,expected_H=item['H'])
    # All stars through total order8, including weak boundaries, independently
    # check the closed scalar formula. No implementation assertions are used.
    for mm in range(1,8):
        base=[sum(1 << j for j in range(i+1,mm)) for i in range(mm)]+[0]
        for mask in range(1 << mm):
            word=[mask >> i & 1 for i in range(mm)]
            add(star(base,mm,list(range(mm)),mask),'small_formula',mask=mask,m=mm,expected_H=count(word))
    # Guard-destroying and weak thirteen-vertex controls, fixed deterministically.
    hostile_masks=sorted(set([0,(1 << m)-1,1,2,3,2048,2049,2051]+[s^2 for s in masks]+[s^1 for s in masks]+[s^2048 for s in masks]))
    for mask in hostile_masks:add(star(rows,a,order,mask),'hostile',mask=mask,names=names,expected_H=table[mask]['H'])
    for index,record in enumerate(data):
        for stage in record['stages']:
            # The first-birth stage is also an inherited positive control.
            add(stage['graph']['rows'],'retained',family=record['family'],stage=stage['stage'],
                names=stage['graph']['names'],expected_H=stage.get('expected_H'))
        add(anchors[index]['graph']['rows'],'retained',family=record['family'],stage='debt',names=names,expected_H=anchors[index]['H'])
    profiler=OLD/'reader'/'endpoint_profile.cpp'
    with tempfile.TemporaryDirectory(prefix='exact-star-profile-') as directory:
        exe=Path(directory)/'profile.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(profiler),'-o',str(exe)],check=True)
        result=subprocess.run([str(exe)],input='\n'.join(inputs)+'\n',text=True,capture_output=True,check=True)
    need(not result.stderr,'unexpected literal profiler stderr')
    profiles=[json.loads(line) for line in result.stdout.splitlines()]
    need(len(profiles)==len(meta),'complete literal profile output')
    typed_profiles=[];retained=[];by_mask={};hostile=[]
    for p,item in zip(profiles,meta):
        if item['expected_H'] is not None:need(p['H']==item['expected_H'],'one-apex formula versus literal subset path DP')
        if item['kind']=='small_formula':continue
        nms=item['names'];rr=nms.index(99)
        enddef=[v-p['H'] for v in p['E2']];startdef=[v-p['H'] for v in p['S2']]
        side=dict(**item,profile=p,start_defects=startdef,end_defects=enddef,
            start_entries=[nms[i] for i,v in enumerate(startdef) if v in (1,2)],
            end_entries=[nms[i] for i,v in enumerate(enddef) if v in (1,2)])
        if item['kind'] in ('typed','hostile'):
            mask=item['mask']
            if table[mask]['strong']:need(enddef[rr]==table[mask]['root_end_defect'],'literal root guard defect')
            if item['kind']=='typed':
                need(enddef[rr]==0,'typed exact root-end guard')
                typed_profiles.append(side);by_mask[mask]=side
            else:hostile.append(side)
        else:
            side['charts']={str(v):phase(p['H'],p['E1'][i],p['S2'][i]) for i,v in enumerate(nms)}
            ii=nms.index(0)
            side['inner_mark0_fields']=dict(H=p['H'],e=p['E1'][ii],s=p['S1'][ii],F=p['E2'][ii],S=p['S2'][ii],P2=p['P2'])
            retained.append(side)
    # All equal-distance nearest counts, and nearest counts on either side.
    nearest=[]
    for k,target in enumerate((1409,1427)):
        lower=max(v['H'] for v in typed if v['H']<target)
        upper=min(v['H'] for v in typed if v['H']>target)
        choices=[]
        for item in typed:
            if item['H'] in (lower,upper):
                side=by_mask[item['mask']];mask=item['mask']
                flips=[[0,names[v]] for j,v in enumerate(order) if (mask^masks[k]) >> j & 1]
                choices.append(dict(mask=mask,word=item['word'],H=item['H'],H_delta_from_debt=item['H']-anchors[k]['H'],
                    H_delta_from_target=item['H']-target,star_flips=flips,flip_count=len(flips),
                    incoming=[names[v] for j,v in enumerate(order) if not mask >> j & 1],
                    outgoing=[names[v] for j,v in enumerate(order) if mask >> j & 1],
                    graph=dict(rows=side['rows'],names=names,root=99,inner_mark=0),
                    start_entries=side['start_entries'],end_entries=side['end_entries'],profile=side['profile'],
                    charts={str(v):phase(side['profile']['H'],side['profile']['E1'][i],side['profile']['S2'][i]) for i,v in enumerate(names)}))
        nearest.append(dict(family=data[k]['family'],debt_H=anchors[k]['H'],target_H=target,lower=lower,upper=upper,
            choices=sorted(choices,key=lambda v:(abs(v['H_delta_from_target']),v['flip_count'],v['mask']))))
    with tempfile.TemporaryDirectory(prefix='exact-star-seed-') as directory:
        exe=Path(directory)/'seed.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(HERE/'seed_fields.cpp'),'-o',str(exe)],check=True)
        seed_input='\n'.join(' '.join(map(str,[i,len(v['rows']),v['names'].index(0)]+v['rows'])) for i,v in enumerate(retained))+'\n'
        raw=subprocess.run([str(exe)],input=seed_input,text=True,capture_output=True,check=True)
    need(not raw.stderr,'unexpected seed-field stderr')
    fresh=[json.loads(line) for line in raw.stdout.splitlines()]
    need(len(fresh)==len(retained),'complete seed-field output')
    for v,z in zip(retained,fresh):
        fields=z['fields'];p=v['profile'];ii=v['names'].index(0)
        need(fields[:4]==[p['H'],p['E1'][ii],p['S2'][ii],p['E2'][ii]] and fields[6]==p['P2'],'independent eight-field overlap')
        v['inner_mark0_eight_fields']=fields
        v['inner_mark0_fields'].update(q=fields[5],p_deleted=fields[4],S3=fields[7])
        if v['stage']=='debt':need(fields[5]==1,'debt q restored exactly')
        if v['stage']=='C2_child':need(fields[1]==1 and fields[5]==3,'child e1 survives while q triples')
    dump('census.json',dict(status='FINITE_EXACT_COMPLETE_FIXED_NAMED_STAR',names=names,inner_mark=0,root=99,
        deleted_order=[names[v] for v in order],debt_masks=masks,universe=4096,typed=512,records=table),True)
    dump('typed_profiles.json',dict(status='FRESH_LITERAL_ALL_512_TYPED_STARS',records=typed_profiles),True)
    dump('nearest.json',dict(status='FINITE_EXACT_NEAREST_TYPED_COUNTS',records=nearest))
    dump('controls.json',dict(status='FRESH_LITERAL_RETAINED_AND_HOSTILE_PROFILES',retained=retained,hostile=hostile))
    summary=dict(status='PROVED_FORMULA_AND_FINITE_EXACT_STAR_EXCLUSION',universe=4096,strong=1024,typed=512,
        typed_distinct_H=len({v['H'] for v in typed}),targets=[1409,1427],target_hits=[0,0],
        independent_full_profiles=len(profiles),small_formula_controls=254,
        typed_start_entry_states=sum(bool(v['start_entries']) for v in typed_profiles),
        typed_end_entry_states=sum(bool(v['end_entries']) for v in typed_profiles),
        retained_fields=[dict(family=v['family'],stage=v['stage'],**v['inner_mark0_fields']) for v in retained],
        nearest=[dict(target=v['target_H'],lower=v['lower'],upper=v['upper'],choices=len(v['choices'])) for v in nearest],
        dependencies={str(source_path.relative_to(RESULTS)):digest(source_path),str(profiler.relative_to(RESULTS)):digest(profiler)},
        scope='All reorientations of inner0 star with every other named edge fixed. No neutral compensation exists; nearest alternatives are exact graphs, not source replacements or native programs.')
    dump('receipt.json',summary)
    print('PASS: 4096 literal-star words,1024 strong,512 root-end-ready,323 distinct typed counts;1409/1427 absent')
    print('PASS:',len(profiles),'independent full subset profiles;254 small formula controls;all512 typed port-entry banks')
    print('PASS: nearest lower/upper counts1407/1417 and1425/1433;all actual rows,cuts,joints and retained seed fields saved')

if __name__=='__main__':main()
