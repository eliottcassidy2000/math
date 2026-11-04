"""Three small local seed menus at one failed recursive target."""
from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
sys.path.insert(0,str(RESULTS/'zero_ready_component_returns_20261103/inverse'))
import reader
from inverse import inverse,need
from enlarge import PARENT_ROWS,PARENT_ROOT,TARGET,options

def witness(operation,rows,root,fields,hit):
    output=reader.run(fields,hit['word']);component=reader.graph(rows,root,hit['word'])
    need(output[0]+(1<<hit['prefix'])*output[3]==TARGET,'Accepted local reader replay')
    return dict(operation=operation,root=root,seed_rows=rows,seed_fields=fields,seed_strong=reader.strong(rows),component_strong=reader.strong(component),fields=output,output_order=len(component)+hit['prefix']+1,**hit)

def inflate(rows,v):
    n=len(rows);slots=[(i,0) for i in range(n)]+[(v,1),(v,2)];out=[0]*len(slots)
    for i,(a,aa) in enumerate(slots):
        for j,(b,bb) in enumerate(slots):
            if (a!=b and rows[a]>>b&1) or (a==b and (bb-aa)%3==1):out[i]|=1<<j
    return tuple(out)

def menu_report():
    first=[];first_count=0
    for operation in options():
        rows=operation['rows']
        for w,f in enumerate(reader.all_fields(rows)):
            hits,_=inverse(f,TARGET);first_count+=1
            first += [witness(operation,rows,w,f,z) for z in hits]
    inflations=[];inflation_count=0
    for vertex in range(6):
        rows=inflate(PARENT_ROWS,vertex);need(reader.strong(rows),'Strong cyclic inflation')
        for w,f in enumerate(reader.all_fields(rows)):
            hits,_=inverse(f,TARGET,minimum_clock=2);inflation_count+=1
            inflations += [witness(dict(kind='C3_inflation',vertex=vertex),rows,w,f,z) for z in hits]
    second=[];second_count=0;restricted_count=0
    for out1 in (11,37):
        one=reader.birth(PARENT_ROWS,out1)
        for out2 in range(128):
            rows=reader.birth(one,out2)
            if not reader.strong(rows):continue
            for w,f in enumerate(reader.all_fields(rows)):
                hits,_=inverse(f,TARGET,minimum_clock=2);second_count+=1
                if w in (5,6,7):restricted_count+=1
                second += [witness(dict(kind='two_births',out_masks=[out1,out2]),rows,w,f,z) for z in hits]
    need(first_count==538 and len(first)==2,'Exact first local menu')
    need(inflation_count==48 and not inflations,'C3 inflation clock>=2 hostile')
    need(second_count==2016 and restricted_count==756 and len(second)==3,'Exact two-birth menu')
    need(not any(z['root'] in (5,6,7) for z in second),'Old mark and newborn-only return guard fails')
    return dict(target=TARGET,parent_seed_rows=PARENT_ROWS,parent_root=PARENT_ROOT,
      first_menu=dict(scope='15 single arc reversals or64 one-birth cuts; every resulting mark; every clock>=1 and every ending0 word; empty word included.',rooted_candidates=first_count,witnesses=first),
      inflation_menu=dict(scope='C3 inflation of one of6 seed vertices; all8 marks; clocks>=2 and every ending0 word; empty included.',rooted_candidates=inflation_count,witnesses=inflations),
      second_menu=dict(scope='Two births: first out-mask11 or37, second any of128 cuts; retain strong8-seeds, all8 marks; every clock>=2 and every ending0 word; empty included.',rooted_candidates=second_count,oldmark_or_newborn_rooted_candidates=restricted_count,witnesses=second))

if __name__=='__main__':
    z=menu_report();(HERE/'growth.json').write_text(json.dumps(z,indent=2)+'\n')
    print('Local menus:',[(k,z[k]['rooted_candidates'],len(z[k]['witnesses'])) for k in ('first_menu','inflation_menu','second_menu')])
