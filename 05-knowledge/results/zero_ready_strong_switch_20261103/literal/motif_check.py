"""Exact recurring edit support and the first interior-collar failure."""
from pathlib import Path
import json,subprocess,tempfile
from verify import need,strong,induced
HERE=Path(__file__).resolve().parent

def collar(record):
    names=record['suffix_names'];t=record['suffix_rows']
    r,w,p=(names.index(v) for v in (32,33,40))
    c=[v for v in range(len(names)) if v not in (r,p)]
    return ([v for v in range(len(t)) if t[v]>>r&1]==[w]
            and all(t[p]>>v&1 for v in c) and strong(induced(t,c)))

def main():
    selected=[];queries=[]
    for seed,indices in ((1,(3,5)),(4,(0,1))):
        rows=json.loads((HERE/f'realized_generation2_seed{seed}_depth3_all_prefixes.json').read_text())
        for index in indices:
            record=rows[index]
            edges=[edge for edge in record['edits'] if 39 not in edge]
            degree={v:sum(v in edge for edge in edges) for edge in edges for v in edge}
            need(len(edges)==3 and sorted(degree.values())==[1,1,1,1,2],'P3 disjoint K2 support')
            center=next(v for v,d in degree.items() if d==2)
            fan=[edge for edge in edges if center in edge]
            leaves={v for edge in fan for v in edge if v!=center}
            isolated=next(edge for edge in edges if center not in edge)
            need(not leaves&set(isolated),'fan disjoint from isolated edge')
            survives=collar(record)
            need(survives==((seed,index)!=(4,1)),'exact inner collar failure boundary')
            r=record['suffix_names'].index(32)
            inports=[record['suffix_names'][v] for v,row in enumerate(record['suffix_rows']) if row>>r&1]
            selected.append(dict(seed=seed,index=index,H=record['target_H'],fan_center=center,
                fan_leaves=sorted(leaves),isolated_edge=isolated,inner_collar=survives,
                inner_root_32_inports=inports))
            queries.append(record)
    for i in (1,4):
        need(collar(json.loads((HERE/'realized_depth3.json').read_text())[i]),'selected first generation collar')
    with tempfile.TemporaryDirectory(prefix='motif-count-') as temp:
        exe=Path(temp)/'count.exe'
        subprocess.run(['g++','-O2','-std=c++20',str(HERE/'count.cpp'),'-o',str(exe)],check=True)
        raw=subprocess.check_output([str(exe)],input=''.join(' '.join(map(str,[14,*r['suffix_rows']]))+'\n' for r in queries),text=True)
    for item,record,line in zip(selected,queries,raw.splitlines()):
        h,*endpoint=map(int,line.split())
        roots=[name for name,f in zip(record['suffix_names'],endpoint) if f==h]
        need(roots==([32] if item['inner_collar'] else []),'all possible inner zero-ready roots')
        item.update(B_H=h,B_endpoint_2=endpoint,B_zero_ready_roots=roots)
    output=dict(status='FINITE_EXACT',selected=selected,
        limitation='Edit support forgets directions, the inport and the surrounding relation; it is not a count recurrence.')
    (HERE/'motifs.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,sort_keys=True))

if __name__=='__main__':main()
