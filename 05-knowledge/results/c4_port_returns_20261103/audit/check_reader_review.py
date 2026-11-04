"""Independent named-set native replay, finite reader inverse, fresh profiles."""
from collections import Counter
from copy import deepcopy
import hashlib
from itertools import permutations
import json
from pathlib import Path
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
PACKET=HERE.parent
RESULTS=PACKET.parent
NATIVE=PACKET/'native/reader_repairs'
READER=PACKET/'reader'
BANK=RESULTS/'recursive_seed_switching_20261103/selector/bank.json'


def need(ok,why):
    if not ok: raise ValueError(why)


def graph(rows,names):
    need(len(rows)==len(names)==len(set(names)), 'unique named vertices')
    n=len(names)
    for i,row in enumerate(rows):
        need(row>=0 and row < 1<<n and not row>>i&1,'valid row')
        for j in range(i+1,n): need((row>>j&1)+(rows[j]>>i&1)==1,'tournament pair')
    return {a:{b for j,b in enumerate(names) if rows[i]>>j&1} for i,a in enumerate(names)}


def strong(g):
    for a in g:
        seen={a};todo=[a]
        while todo:
            for b in g[todo.pop()]-seen: seen.add(b);todo.append(b)
        if len(seen)!=len(g): return False
    return True


def birth(g,item):
    x=item['name'];incoming=item['incoming']
    need(x not in g,'fresh birth name')
    need(len(incoming)==len(set(incoming)) and set(incoming)<=set(g),'valid birth cut')
    g[x]=set(g)-set(incoming)-{x}
    for a in incoming:g[a].add(x)


def native(g,word):
    need(len(word)==len(g) and set(word)==set(g),'spanning guard names')
    need(all(b in g[a] for a,b in zip(word,word[1:])),'directed guard')
    # Reverse every arc, then restore the directed path arcs.
    answer={a:{b for b in g if a in g[b]} for a in g}
    for a,b in zip(word,word[1:]): answer[b].remove(a);answer[a].add(b)
    return answer


def replay(record):
    g=graph(record['source']['rows'],record['source']['names'])
    old=set(g);need(99 in old and strong(g),'strong original rooted source')
    peak=len(g)
    for item in record['permanent_births']:birth(g,item);peak=max(peak,len(g))
    permanent=set(g)
    for program in record['programs']:
        before=deepcopy(g);a,b=program['old_arc'];need(b in g[a],'actual oriented old arc')
        born=[]
        for item in program['births']:birth(g,item);born.append(item['name']);peak=max(peak,len(g))
        need(len(program['guards'])==2,'two native moves per gadget')
        for word in program['guards']:g=native(g,word)
        need(program['deletions']==born[::-1],'only new helper deletions')
        for x in program['deletions']:
            need(x not in permanent,'no old or permanent deletion')
            del g[x]
            for adjacent in g.values():adjacent.discard(x)
        before[a].remove(b);before[b].add(a)
        need(g==before,'single arc effect after helper deletion')
    need(old<=set(g) and set(g)==permanent,'all original names retained')
    need(g==graph(record['target_rows'],record['target_names']),'saved endpoint equality')
    need(g==graph(record['declared_target']['rows'],record['declared_target']['names']),'declared named target equality')
    need(strong(g),'strong final target')
    need(peak==record['peak_order'],'peak order tally')
    return peak


def literal_seed_fields(seed):
    rows=seed['rows'];root=seed['root'];n=len(rows);N=1<<n
    h=[0]*N;e=[[0]*n for _ in range(N)];s=[[0]*n for _ in range(N)];h[0]=1
    for mask in range(1,N):
        for word in permutations([i for i in range(n) if mask>>i&1]):
            if all(rows[a]>>b&1 for a,b in zip(word,word[1:])):
                h[mask]+=1;e[mask][word[-1]]+=1;s[mask][word[0]]+=1
    p=[0]*N
    for mask in range(N):
        part=mask
        while True:
            p[mask]+=h[part]*h[mask^part]
            if not part:break
            part=(part-1)&mask
    full=N-1
    return [h[full],e[full][root],sum(s[m][root]*h[full^m] for m in range(N)),
            sum(e[m][root]*h[full^m] for m in range(N)),p[full^(1<<root)],
            h[full^(1<<root)],p[full],sum(s[m][root]*p[full^m] for m in range(N))]


def run(fields,word):
    h,e,s,f,p,q,P,S=fields
    for bit in word:
        if bit=='1':h,e,s,f,p,q,P,S=h,0,2*s,f,2*p,q,2*P,3*S
        else:h,e,s,f,p,q,P,S=h+s-e,q,2*s-q,f+p,2*p,q,2*(P-f+S),3*S-p
    return [h,e,s,f,p,q,P,S]


def inverse(seed,target):
    h,e,s,f,p,q,*_=seed;clock=1;hits=set();tried=0
    while h+(1<<clock)*f<=target:
        residual=target-h-(1<<clock)*f;K=s-q+(1<<clock)*p;d=q-e
        if residual==0:hits.add((clock,''))
        lower=max(1,(residual-d+K+q-1)//(K+q));upper=residual//K
        for r in range(lower,upper+1):
            word=''.join('0' if r>>j&1 else '1' for j in range(r.bit_length()))
            # Evaluate the one-apex word statistic directly in chronological bits.
            base=1;binary=last=0
            for bit in map(int,word):base+=(1-bit)*(binary+last);binary,last=2*binary+bit,bit
            tried+=1
            if K*r+q*(base-1)+d*(r&1)==residual:hits.add((clock,word))
        clock+=1
    return hits,tried


def target_rows(seed,word,clock):
    rows=list(seed['rows']);mark=seed['root']
    def append(out):
        n=len(rows)
        for i in range(n):
            if not out>>i&1:rows[i]|=1<<n
        rows.append(out)
    for bit in word:append(0 if bit=='1' else 1<<mark)
    component=len(rows)
    for _ in range(clock):append((1<<component)-1)
    append(((1<<len(rows))-1)^(1<<mark))
    return rows


def expected_source(family,t):
    prefix=list(range(10,10+t));names=prefix+list(range(5))+[99]
    g={name:set() for name in names};suffix=[30,28,8,16,4] if family=='24' else [30,28,24,16,0]
    ports={0,1,2} if family=='24' else {0,2,4}
    for i,a in enumerate(prefix):g[a].update(prefix[i+1:]+list(range(5)))
    g[99].update(prefix)
    for a in range(5):
        g[a].update(b for b in range(5) if suffix[a]>>b&1)
        if a in ports:g[a].add(99)
        else:g[99].add(a)
    return g


def phase(h,e,s):
    need((s-e)%2==0,'integral decoded chart')
    x,y=(s-e)//2,(s+e)//2;a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
    return dict(H_X_Y=[h,x,y],F49_plus_charge=[a,b,c],F49_frobenius_plus_charge=[(a+4*b)%7,(-b)%7,c],
                F343_vector=[h%7,x%7,y%7],F343_frobenius=[(h+4*x+6*y)%7,2*x%7,(3*x+4*y)%7])


def main():
    repairs=json.loads((NATIVE/'receipt.json').read_text());witnesses=json.loads((READER/'witnesses.json').read_text())
    tail=json.loads((READER/'tail.json').read_text());sidecars=json.loads((READER/'sidecars.json').read_text())
    entries=json.loads((READER/'ENTRY_DEFECTS.json').read_text())
    entry_index={(r['family'],r['t']):r for r in entries['records']}
    bank=json.loads(BANK.read_text())['seeds'];need(len(bank)==219,'fixed bank size')
    for seed in bank:
        need(strong(graph(seed['rows'],list(range(len(seed['rows']))))),'strong actual seed')
        need(literal_seed_fields(seed)==seed['fields'],'literal eight-field seed audit')
    witness_map={(r['family'],r['t']):r for r in witnesses['records']}
    need(set(witness_map)=={(f,t) for f in ('24','44') for t in range(1,11)},'exact twenty-case scope')
    tried=0;hits_checked=0
    for record in witnesses['records']+tail['records']:
        declared={(h['bank_index'],h['prefix'],h['word']) for h in record['all_hits']};found=set()
        for index,seed in enumerate(bank):
            hits,work=inverse(seed['fields'],record['target_H']);tried+=work
            found.update((index,clock,word) for clock,word in hits)
        need(found==declared and len(found)==record['hit_count'],'complete finite inverse hit set')
        orders=[]
        for hit in record['all_hits']:
            seed=bank[hit['bank_index']];fields=run(seed['fields'],hit['word'])
            need(fields==hit['fields'] and fields[0]+(1<<hit['prefix'])*fields[3]==record['target_H'],'literal affine state replay')
            order=len(seed['rows'])+len(hit['word'])+hit['prefix']+1
            need(order==hit['order'],'word order');orders.append(order);hits_checked+=1
        need(min(orders)==record['minimum_order_in_fixed_bank'],'bank-scoped minimum order')
    need({(r['family'],r['t']) for r in tail['records']}=={(f,t) for f in ('24','44') for t in range(11,17)},'bounded tail scope')
    for family in ('24','44'):
        for start in range(1,5):
            shared=None
            for t in range(start,11):
                keys={(h['bank_index'],h['prefix']) for h in witness_map[family,t]['all_hits']}
                shared=keys if shared is None else shared&keys
            need(not shared,'finite common seed-clock observation')
    need(len(repairs['records'])==20,'twenty native programs')
    costs=Counter();profile_inputs=[];profile_meta=[]
    for r in repairs['records']:
        w=witness_map[r['family'],r['t']];need(r['source']==w['source'] and r['declared_target']==w['target'],'actual source and target provenance')
        need(graph(r['source']['rows'],r['source']['names'])==expected_source(r['family'],r['t']),'specified original sibling source')
        selected=w['selected'];seed=bank[selected['bank_index']]
        need(target_rows(seed,selected['word'],selected['prefix'])==w['target']['rows'],'literal reader target realization')
        need(len(w['target']['rows'])==w['minimum_order_in_fixed_bank'],'selected minimum target')
        peak=replay(r);costs['peak_order']=max(costs['peak_order'],peak)
        costs['permanent_births']+=len(r['permanent_births']);costs['arc_gadgets']+=len(r['programs'])
        for p in r['programs']:
            costs['temporary_births']+=len(p['births']);costs['native_moves']+=len(p['guards']);costs['temporary_deletions']+=len(p['deletions'])
            costs['bridge_gadgets']+=p['case']=='bridge_then_append'
        for role in ('source','target'):
            g=w[role];identifier=len(profile_meta)
            profile_inputs.append(' '.join(map(str,[identifier,len(g['rows'])]+g['rows'])))
            profile_meta.append((r,w,role,g))
    hostile=[]
    for kind in ('guard','birth_name','old_deletion','target_arc'):
        r=deepcopy(repairs['records'][0])
        if kind=='guard':r['programs'][0]['guards'][0].reverse()
        elif kind=='birth_name':r['programs'][0]['births'][0]['name']=r['source']['names'][0]
        elif kind=='old_deletion':r['programs'][0]['deletions'][0]=r['source']['names'][0]
        else:r['declared_target']['rows'][0]^=2;r['declared_target']['rows'][1]^=1
        try:replay(r)
        except ValueError as exc:hostile.append({'corruption':kind,'rejection':str(exc)})
        else:raise ValueError('Corrupted certificate accepted: '+kind)
    with tempfile.TemporaryDirectory(prefix='reader_review_') as tmp:
        exe=Path(tmp)/'profile.exe'
        subprocess.run(['g++','-std=c++20','-O3',str(HERE/'reader_profile.cpp'),'-o',str(exe)],check=True)
        output=subprocess.check_output([str(exe)],input='\n'.join(profile_inputs)+'\n',text=True)
    profiles=[json.loads(line) for line in output.splitlines()];need(len(profiles)==40,'all fresh profiles')
    indexed={(r['family'],r['t'],r['role']):r for r in sidecars['records']}
    for observed,(r,w,role,g) in zip(profiles,profile_meta):
        old=indexed[r['family'],r['t'],role]
        need(observed==old['profile'],'independent full endpoint and joint profile')
        need(old['names']==g['names'],'profile name order')
        h=observed['H'];expected=w['target_H'] if role=='target' else w['source']['H']
        need(h==expected,'actual Hamilton count')
        need(w['target_H']==3*w['source']['H']+4,'actual required digit')
        native_side=r[role+'_sidecar'];need(native_side['H']==h,'native source/target count')
        for i,name in enumerate(g['names']):
            chart=phase(h,observed['E1'][i],observed['S2'][i]);need(chart==old['charts'][str(name)],'all named phase charts')
            if str(name) in native_side['ports']:
                expected_chart=dict(E2=observed['E2'][i],**chart)
                need(expected_chart==native_side['ports'][str(name)],'native selected port chart')
        for key,value in native_side['ordered_J'].items():
            a,b=map(int,key.split(','));need(value==observed['J'][g['names'].index(a)][g['names'].index(b)],'native selected ordered joint')
        root=g['names'].index(99);defect=observed['E2'][root]-h
        need(defect==(0 if role=='target' else 16 if r['family']=='24' else 12),'outer-root defect scope')
        if role=='target':
            entry=entry_index[r['family'],r['t']]
            start={str(name):observed['S2'][i]-h for i,name in enumerate(g['names'])}
            end={str(name):observed['E2'][i]-h for i,name in enumerate(g['names'])}
            need(start==entry['start_defects'] and end==entry['end_defects'],'all-vertex entry defects')
            need(not any(v in (1,2) for v in list(start.values())+list(end.values())),'no immediate positive C2 entry, including converse')
            need([int(k) for k,v in end.items() if v==0]==[99],'unique end-zero root')
    need((NATIVE/'receipt.json').read_bytes()==(NATIVE/'receipt_optimized.json').read_bytes(),'producer normal/-O equality')
    pins=[BANK,READER/'PROOF.md',READER/'witnesses.json',READER/'sidecars.json',READER/'tail.json',READER/'ENTRY_DEFECTS.json',NATIVE/'PROOF.md',NATIVE/'check_repairs.py',NATIVE/'receipt.json',HERE/'reader_profile.cpp',Path(__file__)]
    receipt=dict(status='PASS independent reader and native review',programs=20,finite_targets=32,bank_seeds=219,
        complete_inverse_candidate_integers=tried,accepted_words_replayed=hits_checked,actual_profiles_recomputed=40,
        costs=dict(costs),hostile_rejections=hostile,scope='All-t acceptance, branch regeneration and physical field actions remain unproved.',
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    destination=HERE/('reader_review_optimized.json' if not __debug__ else 'reader_review.json')
    destination.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
