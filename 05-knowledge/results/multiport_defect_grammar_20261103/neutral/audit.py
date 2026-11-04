"""Independent backward counts and Warshall SCC filtering of the repair menu."""
from itertools import combinations
from pathlib import Path
import json
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
OLD=HERE.parent.parent/'zero_ready_strong_switch_20261103'/'literal'

def need(ok,why):
    if not ok:raise RuntimeError(why)

def reach(rows,deleted=None):
    n=len(rows); allowed=((1<<n)-1) ^ (0 if deleted is None else 1<<deleted)
    r=[(row|(1<<v))&allowed if v!=deleted else 0 for v,row in enumerate(rows)]
    for k in range(n):
        for v in range(n):
            if r[v]>>k&1:r[v]|=r[k]
    return r,allowed

def zero_marks(rows):
    reach_all,full=reach(rows)
    if any(x!=full for x in reach_all):return None
    result=[]
    for v,row in enumerate(rows):
        incoming=full^row^(1<<v)
        if incoming.bit_count()==1:
            w=incoming.bit_length()-1
            closure,allowed=reach(rows,v)
            if closure[w]!=allowed:result.append(v)
    return result

def main():
    seed=json.loads((OLD/'realized_generation2_seed4_depth3_all_prefixes.json').read_text())[1]
    names=seed['suffix_names'];base=seed['suffix_rows'];edges=list(combinations(range(14),2))
    cases=[];stats=[]
    for budget in range(4):
        checked=strong=repairable=0
        for edit in combinations(range(len(edges)),budget):
            checked+=1;rows=base[:]
            for e in edit:
                a,b=edges[e];rows[a]^=1<<b;rows[b]^=1<<a
            roots=zero_marks(rows)
            if roots is None:continue
            strong+=1
            if not roots:continue
            repairable+=1;cases.append((edit,rows,roots))
        stats.append({'exact_edges':budget,'examined':checked,'strong':strong,'repairable':repairable})
    parent=json.loads((OLD/'realized_depth3.json').read_text())[4]
    queries=[base,parent['suffix_rows']]+[row for _,row,_ in cases]
    with tempfile.TemporaryDirectory(prefix='neutral-backward-') as temp:
        exe=Path(temp)/'count.exe'
        subprocess.run(['g++','-O3','-std=c++20',str(OLD/'count.cpp'),'-o',str(exe)],check=True)
        raw=subprocess.check_output([str(exe)],input=''.join(' '.join(map(str,[len(row),*row]))+'\n' for row in queries),text=True)
    values=[list(map(int,line.split())) for line in raw.splitlines()]
    need(len(values)==len(queries),'all count queries answered')
    failed,positive=values[:2];r32=names.index(32)
    need(failed[0]==80319 and failed[1+r32]==83279,'actual defect2960')
    need(min(f-failed[0] for f in failed[1:])==2960,'no zero-ready child root')
    need(positive[0]==positive[1+r32]==46073,'actual positive inner collar')
    target=seed['target_H'];minimum=len(seed['prefix_names']);nearest=10**30;hits=[]
    for (edit,rows,roots),counts in zip(cases,values[2:]):
        h,*fs=counts
        need([r for r,f in enumerate(fs) if f==h]==roots,'independent DP checks component roots')
        for w,f in enumerate(fs):
            s=minimum
            while True:
                value=h+(1<<s)*f;nearest=min(nearest,abs(value-target))
                if value==target:hits.append({'edges':edit,'mark':w,'prefix':s})
                if value>target:break
                s+=1
    discovery=json.loads((HERE/'depth3.json').read_text())
    need(sum(s['examined'] for s in stats)==discovery['examined']==125672,'complete edit universe')
    need(sum(s['strong'] for s in stats)==discovery['strong']==119417,'independent strong count')
    need(len(cases)==discovery['with_zero_ready_mark']==4324,'independent inner repair count')
    need(hits==discovery['hits']==[] and nearest==discovery['nearest_distance']==80,'independent exclusion')
    result={'status':'FINITE-EXACT_INDEPENDENT', 'exact_budget_statistics':stats,
            'counted_repairable_graphs':len(cases),'target_H':target,'minimum_prefix':minimum,
            'hits':hits,'nearest_distance':nearest,'positive_control_H':46073,'hostile_minimum_defect':2960,
            'methods':['independent Warshall closure','backward endpoint deletion recurrence',
                       'exact monotone prefix equations on every actual mark']}
    (HERE/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
