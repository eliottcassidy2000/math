"""Three finite rounds with literal names, inverse seed records and native lifts."""
from pathlib import Path
import sys,json,hashlib,itertools,time
HERE=Path(__file__).resolve().parent
PACKET=HERE.parent; RESULTS=PACKET.parent
sys.path.insert(0,str(PACKET/'inverse'))
import reader
from reader import need

def seed_bank():
    small,labelled=reader.small_bank();bank={tuple(z['fields']):z for z in small};six=0
    for line in (RESULTS/'metagraph_primegame_20261102/low_sources.txt').read_text().splitlines():
        z=list(map(int,line.split()));n=z[0]
        if n!=6:continue
        rows=tuple(z[2:2+n]);need(reader.strong(rows),'Inherited six-vertex strong seed');six+=1
        for w,f in enumerate(reader.all_fields(rows)):
            bank.setdefault(tuple(f),dict(rows=rows,root=w,fields=f))
    need(len(bank)==219 and six==35,'Exact inherited seed-bank size')
    return list(bank.values()),dict(eight_field_keys=len(bank),labelled_strong_through5=labelled,inherited_six_vertex_representatives=six)

def transitive(rows,vs):
    chain=sorted(vs,key=lambda v:-sum(rows[v]>>u&1 for u in vs))
    return all(rows[u]>>v&1 for i,u in enumerate(chain) for v in chain[i+1:])

def structural_feedback(seed,root,spare):
    n=len(seed);required={root,spare};others=[i for i in range(n) if i not in required]
    for k in range(len(others)+1):
        for extra in itertools.combinations(others,k):
            choice=required|set(extra)
            if transitive(seed,[i for i in range(n) if i not in choice]):return sorted(choice)
    raise RuntimeError('No feedback set')

def build_named(parent,hit,H,offset,depth,branch):
    seed=hit['seed'];word=hit['word'];prefix_size=hit['prefix'];seed_rows=tuple(seed['rows']);w=seed['root'];m=len(seed_rows)
    output_order=m+len(word)+prefix_size+1
    names=list(parent['names']);need(output_order>=len(names),'All physical names fit')
    new=[];fresh=max(names)+1
    while len(names)<output_order:
        need(fresh<1000,'Permanent names stay outside temporary namespace')
        new.append(fresh);names.append(fresh);fresh+=1
    available=[v for v in names if v not in (39,31,6,50,51)]
    prefix=[50,51]+available[:prefix_size-2];available=available[prefix_size-2:]
    seed_names=[None]*m;seed_names[w]=31
    spare=next(i for i in range(m) if i!=w);seed_names[spare]=6
    for i in range(m):
        if seed_names[i] is None:seed_names[i]=available.pop(0)
    tail_names=available;need(len(tail_names)==len(word),'Every old/new name assigned once')
    component_names=seed_names+tail_names
    component=reader.graph(seed_rows,w,word)
    edges=set()
    for i,a in enumerate(component_names):
        for j,b in enumerate(component_names):
            if component[i]>>j&1:edges.add((a,b))
    for i,a in enumerate(prefix):
        for b in prefix[i+1:]+component_names:edges.add((a,b))
    for v in names:
        if v!=39:edges.add((31,39) if v==31 else (39,v))
    rows=tuple(sum(1<<j for j,b in enumerate(names) if (a,b) in edges) for a in names)
    need(reader.strong(component) and reader.strong(rows),'Strong component and full wrapper')
    loc={v:i for i,v in enumerate(names)};ploc={v:i for i,v in enumerate(parent['names'])}
    changes=[[a,b] for i,a in enumerate(parent['names']) for b in parent['names'][i+1:]
             if bool(parent['rows'][ploc[a]]>>ploc[b]&1)!=bool(rows[loc[a]]>>loc[b]&1)]
    fb=structural_feedback(seed_rows,w,spare)
    record=dict(H=H,parent_H=parent['H'],offset=offset,depth=depth,branch=branch,
        names=names,rows=rows,root=39,mark=31,spare=6,prefix=prefix,
        component_names=component_names,component_rows=component,
        seed_reader=dict(rows=seed_rows,root=w,names=seed_names,word=word,fields=seed['fields']),
        component_eight_fields=hit['fields'],new_names=new,
        old_name_identity=[[v,v] for v in parent['names']],old_arc_differences=changes,
        structural_feedback_names=[39]+[seed_names[i] for i in fb])
    need(hit['fields'][0]+(1<<prefix_size)*hit['fields'][3]==H,'Exact inverse replay')
    return record

def main():
    bank,description=seed_bank();fixed=json.loads((PACKET/'inverse/named_targets.json').read_text())
    source=dict(H=170807,names=fixed['source_names'],rows=fixed['source_rows'],branch='')
    records=[];sources={170807:source}
    for x in fixed['targets']:
        x=dict(x);x.update(parent_H=170807,depth=1,branch=str(x['offset']),old_name_identity=[[v,v] for v in source['names']])
        records.append(x);sources[x['H']]=x
    scouts=[]
    for depth in (2,3):
        parents=[x for x in records if x['depth']==depth-1]
        for parent in parents:
            for offset in (2,4):
                target=3*parent['H']+offset;hits=[];tried=0
                for seed in bank:
                    found,count=reader.inverse(seed['fields'],target);tried+=count
                    hits += [dict(seed=seed,**z) for z in found]
                admissible=[z for z in hits if z['prefix']>=2 and len(z['seed']['rows'])+len(z['word'])+z['prefix']+1>=len(parent['names'])]
                need(admissible,'At least one name-retaining prefix>=2 witness for '+str(target))
                def key(z):
                    seed=z['seed'];sp=next(i for i in range(len(seed['rows'])) if i!=seed['root'])
                    fb=structural_feedback(seed['rows'],seed['root'],sp)
                    return (len(seed['rows'])+len(z['word'])+z['prefix']+1,len(fb),len(z['word']),z['word'],z['prefix'],seed['rows'],seed['root'])
                best=min(admissible,key=key)
                rec=build_named(parent,best,target,offset,depth,parent['branch']+str(offset))
                records.append(rec);sources[target]=rec
                scouts.append(dict(H=target,total_hits=len(hits),admissible_hits=len(admissible),candidate_integers=tried,selected_order=len(rec['names'])))
                print('Named',target,'depth',depth,'order',len(rec['names']),'changes',len(rec['old_arc_differences']),flush=True)
    need(len(records)==14 and all(len([x for x in records if x['parent_H']==h])==2 for h in [170807]+[x['H'] for x in records if x['depth']<3]),'Complete depth-three binary tree')
    pins=[PACKET/'inverse/reader.py',PACKET/'inverse/named_targets.json',RESULTS/'metagraph_primegame_20261102/low_sources.txt',Path(__file__)]
    report=dict(status='FINITE_EXACT_NAMED_THREE_ROUNDS',scope='Complete binary3H+2/+4 tree through depth3 from actualH170807; arbitrary-depth closure OPEN.',seed_bank=description,source=source,targets=records,inverse_search=scouts,source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    (HERE/'named_tree.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
if __name__=='__main__':main()
