"""Two guarded named endpoints after the repair; no additional target layer."""
from pathlib import Path
import json,hashlib,itertools
from check import previous,bank_accept,reader,need
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]

def feedback(rows,names,seed_names,required):
    loc={v:i for i,v in enumerate(names)};others=[v for v in seed_names if v not in required]
    for k in range(len(others)+1):
        for extra in itertools.combinations(others,k):
            chosen=set(required)|set(extra);vs=[i for i,v in enumerate(names) if v not in chosen]
            chain=sorted(vs,key=lambda i:-sum(rows[i]>>j&1 for j in vs))
            if all(rows[u]>>v&1 for i,u in enumerate(chain) for v in chain[i+1:]):return [v for v in names if v in chosen]
    raise RuntimeError('Feedback set absent')

def main():
    source=json.loads((HERE/'named_repair.json').read_text())['target'];bank,description=previous.seed_bank();targets=[]
    for offset in (2,4):
        H=3*source['H']+offset;hits,stats=bank_accept(bank,H)
        guarded=[z for z in hits if z['prefix']>=2 and len(z['seed']['rows'])+len(z['word'])+z['prefix']+1>=len(source['names'])]
        need(guarded,'Actual old names fit a clock>=2 witness')
        exact_clock=[z for z in guarded if z['prefix']==len(source['prefix'])]
        z=min(exact_clock or guarded,key=lambda z:(len(z['seed']['rows'])+len(z['word'])+z['prefix']+1,len(z['word']),z['word'],z['prefix']))
        seed=z['seed'];seed_rows=tuple(seed['rows']);w=seed['root'];word=z['word'];k=z['prefix'];n=len(seed_rows)+len(word)+k+1
        names=list(source['names']);new=list(range(max(names)+1,max(names)+1+n-len(names)));names+=new
        available=[v for v in names if v not in (39,4,6,50,51)]
        prefix=[50,51]+available[:k-2];available=available[k-2:]
        if k==4:need(prefix==source['prefix'],'Same literal four-prefix')
        seed_names=[None]*len(seed_rows);seed_names[w]=4;spare=next(i for i in range(len(seed_rows)) if i!=w);seed_names[spare]=6
        for i in range(len(seed_rows)):
            if seed_names[i] is None:seed_names[i]=available.pop(0)
        need(len(available)==len(word),'Exact name capacity');component_names=seed_names+available
        component=reader.graph(seed_rows,w,word);edges=set()
        for i,a in enumerate(component_names):
            for j,b in enumerate(component_names):
                if component[i]>>j&1:edges.add((a,b))
        for i,a in enumerate(prefix):
            for b in prefix[i+1:]+component_names:edges.add((a,b))
        for v in names:
            if v!=39:edges.add((4,39) if v==4 else (39,v))
        rows=tuple(sum(1<<j for j,b in enumerate(names) if (a,b) in edges) for a in names)
        need(reader.strong(component) and reader.strong(rows),'Literal strong endpoint')
        loc={v:i for i,v in enumerate(names)}
        differences=[[a,b] for i,a in enumerate(source['names']) for b in source['names'][i+1:] if bool(source['rows'][i]>>loc[b]&1)!=bool(rows[i]>>loc[b]&1)]
        fb=feedback(rows,names,seed_names,{39,4});jointfb=feedback(rows,names,seed_names,{39,4,31,6})
        targets.append(dict(H=H,parent_H=source['H'],offset=offset,names=names,rows=rows,root=39,mark=4,previous_mark=4,retained_old_mark=31,spare=6,prefix=prefix,
            seed_reader=dict(rows=seed_rows,root=w,names=seed_names,fields=seed['fields'],word=word),component_names=component_names,component_rows=component,component_eight_fields=z['fields'],
            old_name_identity=[[v,v] for v in source['names']],new_names=new,deleted_old_names=[],old_arc_differences=differences,structural_feedback_names=fb,all_retained_port_feedback_names=jointfb,
            bank_hits=len(hits),name_capacity_clock_ge2_hits=len(guarded),same_clock_hits=len(exact_clock),bank_index=z['bank_index']))
        print(H,'order',n,'clock',k,'new',new,'changes',len(differences),'sameclock_hits',len(exact_clock),flush=True)
    pins=[HERE/'named_repair.json',HERE/'inverse.py',HERE/'check.py',HERE/'followup_named.py']
    report=dict(status='FINITE_EXACT_NAMED_READER_TEMPLATES',scope='Both selected next targets have actual named graph templates retaining every old name. Native programs and independent full graph counts are separate obligations; no further target layer.',source=source,targets=targets,source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    (HERE/'followup_named.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
