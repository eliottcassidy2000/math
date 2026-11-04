"""Exact independent interval controls, local menus, and selected recursive policy trace."""
from pathlib import Path
import sys,json,hashlib,importlib.util
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
OLD=RESULTS/'zero_ready_component_returns_20261103'
sys.path.insert(0,str(OLD/'inverse'));import reader
from inverse import need,excess,evaluate,sections,inverse
from growth import menu_report
spec=importlib.util.spec_from_file_location('previous_builder',OLD/'continuation/build.py');previous=importlib.util.module_from_spec(spec);spec.loader.exec_module(previous)

def canonical_graph(seed,word,clock):
    rows=reader.graph(tuple(seed['rows']),seed['root'],word);n=len(rows);out=list(rows)
    for j in range(clock):out.append((1<<n)-1|sum(1<<(n+k) for k in range(j+1,clock)))
    out=reader.birth(tuple(out),((1<<len(out))-1)^(1<<seed['root']))
    need(reader.strong(out),'Actual canonical strong witness')
    return dict(rows=out,names=list(range(len(out))),root=len(out)-1,mark=seed['root'],prefix=list(range(n,n+clock)),scope='Literal canonical names; not a parent-to-child physical alignment.')

def bank_accept(bank,T):
    hits=[];stats={k:0 for k in ('states','bisections','evaluations','gcd_prunes','negative_prunes')};depth=0;clocks=0
    for i,seed in enumerate(bank):
        found,cs=inverse(seed['fields'],T)
        for z in found:
            output=reader.run(seed['fields'],z['word']);need(output[0]+(1<<z['prefix'])*output[3]==T,'Independent state replay')
            hits.append(dict(bank_index=i,seed=seed,fields=output,**z))
        clocks+=len(cs)
        for c in cs:
            for key in stats:stats[key]+=c[key]
            depth=max(depth,c['max_depth'])
    hits.sort(key=lambda z:(len(z['seed']['rows'])+len(z['word'])+z['prefix']+1,z['prefix']<2,len(z['word']),z['word'],z['prefix']))
    return hits,dict(clocks=clocks,max_section_depth=depth,**stats)

def main():
    need(all(excess(r)==reader.word_h(''.join('0' if r>>j&1 else '1' for j in range(r.bit_length())))-1 for r in range(1,2049)),'Low-bit excess versus inherited streaming evaluator')
    tests=0
    for K in range(1,9):
        for q in range(1,6):
            for d in sorted({0,q//2,q}):
                for D in range(129):
                    got,_=sections(K,q,d,D);want=[r for r in range(D//K+1) if evaluate(K,q,d,r)==D]
                    need(got==want,'All section solutions versus literal integer enumeration');tests+=1
    need(evaluate(1,1,0,6)==10 and evaluate(1,1,0,7)==7,'Global monotonicity hostile')
    bank,description=previous.seed_bank();(HERE/'bank.json').write_text(json.dumps(dict(description=description,seeds=bank),separators=(',',':'))+'\n')
    source_tree=json.loads((OLD/'continuation/named_tree.json').read_text());old_targets=[z['H'] for z in source_tree['targets']]
    inherited=[]
    for T in old_targets:
        hits,stats=bank_accept(bank,T);need(hits,'Inherited three-round acceptance');inherited.append(dict(H=T,hits=len(hits)))
    chosen=[13835447,13835449,13835525,13835527,41506343,41506345];records=[];cache={}
    for T in chosen:
        hits,stats=bank_accept(bank,T);cache[T]=hits
        admissible=[z for z in hits if z['prefix']>=2 and len(z['seed']['rows'])+len(z['word'])+z['prefix']+1>=23]
        saved=[]
        for z in (admissible or hits)[:1]:
            saved.append(dict(**z,actual_graph=canonical_graph(z['seed'],z['word'],z['prefix'])))
        records.append(dict(H=T,hits=len(hits),same_name_capacity_clock_ge2_hits=len(admissible),statistics=stats,selected=saved))
    need(not cache[13835447] and len(cache[13835449])==5,'First fixed-bank failure and sibling control')
    need(len(cache[41506343])==len(cache[41506345])==9,'Selected next pair returns to original bank')
    old_hits=[];tried=0
    for seed in bank:
        got,n=reader.inverse(seed['fields'],13835447);old_hits+=got;tried+=n
    need(not old_hits and tried==1473275,'Independent complete interval inversion confirms fixed-bank miss')
    growth=menu_report();(HERE/'growth.json').write_text(json.dumps(growth,indent=2)+'\n')
    new_seed=dict(rows=(38,124,24,241,97,4,165,55),root=1,fields=(269,19,679,338,940,59,4520,5246))
    need(reader.all_fields(new_seed['rows'])[1]==new_seed['fields'],'Literal new8seed field count')
    repair,repair_stats=inverse(new_seed['fields'],13835447)
    need(any(z['word']=='1111000100' and z['prefix']==4 for z in repair),'New seed accepts exact repair')
    alone=[]
    for T in (41506343,41506345):
        h,c=inverse(new_seed['fields'],T);need(not h,'New8seed alone does not repeat');alone.append(dict(H=T,hits=h))
    boundary=[]
    for z in growth['first_menu']['witnesses']:
        seed=dict(rows=z['seed_rows'],root=z['root'],fields=z['seed_fields'])
        boundary.append(dict(**z,actual_graph=canonical_graph(seed,z['word'],z['prefix'])))
    pins=[HERE/'inverse.py',HERE/'growth.py',HERE/'enlarge.py',HERE/'check.py',HERE/'named.py',HERE/'named_repair.json',OLD/'inverse/reader.py',OLD/'continuation/build.py',OLD/'continuation/named_tree.json',RESULTS/'metagraph_primegame_20261102/low_sources.txt']
    report=dict(status='PASS',scope='Fixed219-bank counterexample at fourth all+2 child; one local seed enlargement repairs it, and the next selected pair returns to the old bank. No whole higher layer or arbitrary-depth closure claim.',section_controls=tests,streaming_excess_controls=2048,global_monotonicity_hostile=dict(K=1,q=1,d=0,values=[[6,10],[7,7]]),bank_description=description,inherited_three_round_targets=inherited,selected_targets=records,
        independent_miss=dict(H=13835447,hits=0,candidate_integers=tried,method='Inherited complete interval inversion'),new_seed=new_seed,new_seed_repair_hits=repair,new_seed_alone_next_pair=alone,clock1_boundary_graphs=boundary,
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
