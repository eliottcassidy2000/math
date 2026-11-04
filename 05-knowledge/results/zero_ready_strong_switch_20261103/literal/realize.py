"""Realize a suffix witness by isolated edge repairs and prefix births."""
from pathlib import Path
from itertools import permutations
import importlib.util,json,sys

HERE=Path(__file__).resolve().parent
OLD=HERE.parents[1]/'triangle_return_170807_20261103'
spec=importlib.util.spec_from_file_location('old_realize',OLD/'realize.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)

def need(ok,why):
    if not ok: raise RuntimeError(why)

def induced(T,vs):
    return tuple(sum(((T[a]>>b)&1)<<j for j,b in enumerate(vs)) for a in vs)

def build(hit,seed=None):
    inherited=json.loads((OLD/'interface.json').read_text())
    zero=json.loads((OLD/'zero_ready.json').read_text())[1]
    if seed is None:
        source=tuple(inherited['target']);names=inherited['target_names'][:]
        bnames=zero['suffix_names'];rname=39;oldw=31;prefix=[50,51];source_H=170807
    else:
        source=tuple(seed['target_rows']);names=seed['target_names'][:]
        bnames=seed['suffix_names'];rname=seed['target_root'];oldw=seed['target_inport']
        prefix=seed['prefix_names'][:];source_H=seed['target_H']
    source_names=names[:];source_prefix=prefix[:];neww=bnames[hit['mark']]
    edits=[tuple(bnames[v] for v in edge) for edge in hit['edges']]
    if neww!=oldw:edits.extend(((rname,oldw),(rname,neww)))
    common=[50,51]
    for order in permutations(edits):
        T=source;trace=[]
        for i,edge in enumerate(order):
            a,b=(names.index(x) for x in edge)
            if not (T[a]>>b)&1:a,b=b,a
            ports=tuple(names.index(x) for x in common)
            path=old.guard(T,a,b,ports)
            if path is None:break
            helper=1000+i;n=len(T);full=(1<<n)-1;cut=full^(1<<a)
            X=old.birth(T,cut);k=path.index(a)
            p=(*path,n);q=(*path[:k+1],n,*path[k+1:])
            Y=old.native(X,p);Z=old.native(Y,q);U=induced(Z,range(n))
            expected=list(T);expected[a]^=1<<b;expected[b]^=1<<a
            need(U==tuple(expected),'exact isolated repair')
            if not old.strong(U):break
            def frame(p,v):
                j=p.index(v);return (p[j-1] if j else None,p[j+1] if j+1<len(p) else None)
            need(all(frame(p,v)==frame(q,v) for v in ports),'common unchanged prefix frames')
            local=names+[helper]
            trace.append(dict(arc=[names[a],names[b]],helper=helper,cut=cut,
                              paths=[[local[v] for v in p],[local[v] for v in q]],rows=U))
            T=U
        if len(trace)==len(edits):break
    need(len(trace)==len(edits),'strong native order exists')
    suffix=induced(T,[names.index(v) for v in bnames])
    need(suffix==tuple(hit['rows']),'full suffix witness is exact')
    births=[]
    need(hit['prefix_size']>=len(prefix),'no old prefix deletion')
    for j in range(hit['prefix_size']-len(prefix)):
        newname=max(names)+1
        cut=sum(1<<names.index(v) for v in [rname,*prefix])
        T=old.birth(T,cut);names.append(newname);prefix.append(newname)
        need(old.strong(T),'prefix birth remains strong')
        births.append(dict(name=newname,cut=cut,rows=T))
    need([names[i] for i in range(len(T)) if T[i]>>names.index(rname)&1]==[neww],'unique actual inport')
    need(set(source_names)<=set(names),'all old names survive')
    retained=list(dict.fromkeys([39,6,50,51,31,oldw,neww]))
    changed=sorted(set(retained)&set(v for edge in edits for v in edge))
    return dict(offset=hit['offset'],source_H=source_H,target_H=3*source_H+hit['offset'],
                source_rows=source,source_names=source_names,source_prefix=source_prefix,
                source_root=rname,source_inport=oldw,target_rows=T,target_names=names,
                target_root=39,target_inport=neww,prefix_names=prefix,
                suffix_names=bnames,suffix_rows=suffix,suffix_h=hit['h'],suffix_F=hit['F'],
                edits=[list(e) for e in edits],trace=trace,prefix_births=births,
                common_frame_ports=common,changed_retained_ports=changed,
                retained_ports=retained,deleted_old_names=[],
                native_complements=2*len(trace),
                scope='Finite named branch into zero-ready grammar; no uniform repeated return asserted.')

def main():
    file=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'suffix_depth2_moving1.json'
    data=json.loads(file.read_text())
    seed=None
    if data.get('input_source','').startswith('realized_depth3.json index'):
        seed_index=int(data['input_source'].split('index')[1])
        seed=json.loads((HERE/'realized_depth3.json').read_text())[seed_index]
    out=[build(hit,seed) for hit in data['hits']]
    if seed is not None or data.get('all_prefixes'):
        path=HERE/('realized_'+file.name)
    else:path=HERE/('realized_depth'+str(data['maximum_flips'])+'.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps([{'target':v['target_H'],'native_complements':v['native_complements'],
                      'births':len(v['prefix_births']),'order':len(v['target_rows'])} for v in out]))

if __name__=='__main__':main()
