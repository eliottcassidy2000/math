"""Standalone exact replay, independent of all imported historical checkers."""
from pathlib import Path
import copy,json,subprocess,tempfile,sys
HERE=Path(__file__).resolve().parent

def need(condition,label):
    if not condition:raise RuntimeError(label)

def valid(t):
    n=len(t);full=(1<<n)-1
    need(all(0<=row<=full and not(row>>i&1) for i,row in enumerate(t)),'row type')
    need(all((t[i]>>j&1)+(t[j]>>i&1)==1 for i in range(n) for j in range(i)),'pair orientation')

def strong(t):
    full=(1<<len(t))-1
    for reverse in (False,True):
        reached=1
        while True:
            new=reached
            for v in range(len(t)):
                if reached>>v&1:new|=(full^t[v]^(1<<v)) if reverse else t[v]
            if new==reached:break
            reached=new
        if reached!=full:return False
    return True

def induced(t,vertices):
    return tuple(sum((t[u]>>v&1)<<i for i,v in enumerate(vertices)) for u in vertices)

def birth(t,cut):
    n=len(t);full=(1<<n)-1
    need(0<=cut<=full,'birth cut range')
    return tuple(t[v]|(((cut>>v)&1)<<n) for v in range(n))+(full^cut,)

def native(t,path):
    need(sorted(path)==list(range(len(t))),'Hamiltonian path covers actual names')
    need(all(t[u]>>v&1 for u,v in zip(path,path[1:])),'Hamiltonian path directed in actual intermediate graph')
    full=(1<<len(t))-1
    result=[full^row^(1<<v) for v,row in enumerate(t)]
    for u,v in zip(path,path[1:]):result[u]^=1<<v;result[v]^=1<<u
    valid(result);return tuple(result)

def frame(path,port):
    j=path.index(port)
    return (path[j-1] if j else None,path[j+1] if j+1<len(path) else None)

def replay(record):
    names=record['source_names'][:];t=tuple(record['source_rows']);valid(t)
    need(strong(t),'strong source')
    source=set(names);deleted=[]
    for item in record['trace']:
        n=len(t);helper=item['helper'];need(helper not in names,'fresh temporary helper')
        local=names+[helper];before=t;t=birth(t,item['cut'])
        paths=[[local.index(v) for v in p] for p in item['paths']]
        need(all(frame(paths[0],names.index(v))==frame(paths[1],names.index(v)) for v in record['common_frame_ports']),'exact common frames')
        for path in paths:t=native(t,path)
        t=induced(t,list(range(n)));deleted.append(helper)
        expected=list(before);a,b=(names.index(v) for v in item['arc'])
        need(before[a]>>b&1,'named directed repair arc')
        expected[a]^=1<<b;expected[b]^=1<<a
        need(t==tuple(expected)==tuple(item['rows']),'exact one-edge effect')
        need(strong(t),'completed repair strong')
    for item in record['prefix_births']:
        need(item['name'] not in names,'fresh permanent prefix name')
        t=birth(t,item['cut']);names.append(item['name'])
        need(t==tuple(item['rows']) and strong(t),'strong exact prefix birth')
    need(t==tuple(record['target_rows']) and names==record['target_names'],'literal target')
    need(source<=set(names) and not source&set(deleted),'all old vertices survive')
    need(record['deleted_old_names']==[],'deletion receipt')
    r=names.index(record['target_root']);w=names.index(record['target_inport'])
    need([v for v in range(len(t)) if t[v]>>r&1]==[w],'root has exact sole inport')
    aa=[names.index(v) for v in record['prefix_names']];bb=[names.index(v) for v in record['suffix_names']]
    need(set(aa+bb+[r])==set(range(len(t))) and len(aa)+len(bb)+1==len(t),'grammar partition')
    need(all(t[a]>>b&1 for a in aa for b in bb),'prefix dominates suffix')
    need(all(t[a]>>b&1 for i,a in enumerate(aa) for b in aa[i+1:]),'transitive prefix order')
    suffix=induced(t,bb);need(suffix==tuple(record['suffix_rows']) and strong(suffix),'actual marked strong component')
    return t,suffix

def main():
    path=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'realized_suffix_depth3_moving1_all_prefixes.json'
    records=json.loads(path.read_text());queries=[];states=[]
    for record in records:
        t,b=replay(record);states.append((t,b))
        queries.extend((record['source_rows'],t,b))
    control_start=len(queries)
    queries.extend(((6,4,0),(2,4,1)))
    hostile=[]
    for mutation in ('duplicate_path_vertex','wrong_inport','altered_final_row'):
        broken=copy.deepcopy(records[0])
        if mutation=='duplicate_path_vertex':broken['trace'][0]['paths'][0][0]=broken['trace'][0]['paths'][0][1]
        elif mutation=='wrong_inport':broken['target_inport']=broken['target_root']
        else:broken['target_rows'][0]^=1
        try:replay(broken)
        except RuntimeError:hostile.append(mutation)
        else:raise RuntimeError('hostile input accepted: '+mutation)
    with tempfile.TemporaryDirectory(prefix='literal-check-') as temp:
        exe=Path(temp)/'count.exe'
        subprocess.run(['g++','-std=c++20','-O2',str(HERE/'count.cpp'),'-o',str(exe)],check=True)
        raw=subprocess.check_output([str(exe)],input=''.join(' '.join(map(str,[len(t),*t]))+'\n' for t in queries),text=True)
    counts=[list(map(int,line.split())) for line in raw.splitlines()]
    need(len(counts)==len(queries),'count query completeness')
    need(counts[control_start:]==[[1,1,2,4],[3,3,3,3]],'transitive and cyclic positive count controls')
    receipt=[]
    for i,(record,(t,b)) in enumerate(zip(records,states)):
        source,target,suffix=counts[3*i:3*i+3]
        h=suffix[0];f=suffix[1+record['suffix_names'].index(record['target_inport'])]
        r=record['target_names'].index(record['target_root'])
        need(source[0]==record['source_H'],'independent source H')
        need(target[0]==record['target_H']==3*source[0]+record['offset'],'independent target H')
        need(h==record['suffix_h'] and f==record['suffix_F'],'independent suffix h,F')
        need(h+(1<<len(record['prefix_names']))*f==target[0],'zero-ready formula')
        need(target[1+r]==target[0],'zero-ready root E2=H')
        receipt.append(dict(index=i,H=target[0],order=len(t),h=h,F=f,root=record['target_root'],inport=record['target_inport'],prefix=len(record['prefix_names']),native_complements=2*len(record['trace']),old_names_preserved=len(record['source_names'])))
    print(json.dumps({'status':'PASS','records':receipt,'hostile_rejections':hostile},sort_keys=True))

if __name__=='__main__':main()
