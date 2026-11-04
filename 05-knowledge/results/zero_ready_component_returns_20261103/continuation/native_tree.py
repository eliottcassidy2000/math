"""Full compact native programs for each actual edge of the finite named tree."""
from pathlib import Path
import sys,json,hashlib,gc,time
from functools import lru_cache
HERE=Path(__file__).resolve().parent;PACKET=HERE.parent
sys.path.insert(0,str(PACKET/'inverse'))
import native
from reader import need,birth

def digest(rows):return hashlib.sha256(json.dumps(list(rows),separators=(',',':')).encode()).hexdigest()

class SearchBudget(Exception):pass

def bounded_path(rows,u,v,budget=20000):
    full=(1<<len(rows))-1;visited=0
    @lru_cache(None)
    def walk(mask,last):
        nonlocal visited
        visited+=1
        if visited>budget:raise SearchBudget
        if mask==full:return ()
        available=((1<<v)&~mask) if last==u else ((full if last<0 else rows[last])&~mask&~(1<<v))
        while available:
            bit=available&-available;available-=bit;w=bit.bit_length()-1
            if w==u and mask>>v&1:continue
            tail=walk(mask|bit,w)
            if tail is not None:return (w,)+tail
        return None
    try:return walk(0,-1)
    except SearchBudget:return None

def compile_both(rows,a,b):
    u,v=(a,b) if rows[a]>>b&1 else (b,a)
    for inverse in (False,True):
        T=native.old.flip(rows,u,v) if inverse else rows
        x,y=(v,u) if inverse else (u,v)
        path=bounded_path(T,x,y)
        if path is not None:
            U,V,W,p,q,case=native.old.traced_gadget(T,x,y,path)
            return (W,V,U,q,p,case,True) if inverse else (U,V,W,p,q,case,False)
    return native.compile_arc(rows,a,b)

def replay(source,program,target):
    rows=tuple(source['rows']);names=list(source['names'])
    for z in program['permanent_births']:
        need(z['name'] not in names,'fresh permanent physical name')
        rows=birth(rows,z['out_mask']);names.append(z['name'])
    for z in program['changes']:
        need(z['helper'] not in names,'fresh helper')
        loc=names+[z['helper']];U=birth(rows,z['helper_out_mask'])
        p,q=(tuple(loc.index(v) for v in path) for path in z['paths'])
        V=native.old.native(U,p);W=native.old.native(V,q)
        after=tuple(x&((1<<len(rows))-1) for x in W[:-1])
        a,b=(names.index(v) for v in z['changed_pair'])
        need(after==native.old.flip(rows,a,b),'one old arc and no other endpoint change')
        need(digest(after)==z['after_delete_sha256'],'exact intermediate endpoint hash')
        rows=after
    need(names==target['names'] and rows==tuple(target['rows']),'exact named target endpoint')
    return dict(target_H=target['H'],native_moves=2*len(program['changes']),permanent_births=len(program['permanent_births']),temporary_births=len(program['changes']),deleted_old_names=[],peak_order=len(names)+bool(program['changes']))

def main():
    tree=json.loads((HERE/'named_tree.json').read_text());sources={tree['source']['H']:tree['source']};sources.update({x['H']:x for x in tree['targets']})
    inherited=json.loads((PACKET/'inverse/native.json').read_text());old_programs={p['target_H']:p for p in inherited['programs']}
    programs=[];receipts=[]
    saved=json.loads((HERE/'native_tree.json').read_text()) if (HERE/'native_tree.json').exists() else {'programs':[]}
    completed={p['target_H']:p for p in saved['programs']}
    for index,target in enumerate(tree['targets']):
        source=sources[target['parent_H']];names=list(source['names']);rows=tuple(source['rows']);new=[];moves=[]
        if target['H'] in completed:
            done=completed[target['H']];new=done['permanent_births'];moves=done['changes']
        elif target['depth']==1:
            old=old_programs[target['H']]
            new=[dict(name=z['name'],out_mask=z['out_mask']) for z in old['permanent_births']]
            moves=[dict(changed_pair=z['changed_pair'],helper=z['helper'],helper_out_mask=z['helper_out_mask'],paths=z['paths'],case=z['case'],inverse=z['inverse'],after_delete_sha256=digest(z['after_delete'])) for z in old['changes']]
        else:
            loc={v:i for i,v in enumerate(target['names'])};final=target['rows']
            for v in target['new_names']:
                out=sum(1<<i for i,u in enumerate(names) if final[loc[v]]>>loc[u]&1)
                rows=birth(rows,out);names.append(v);new.append(dict(name=v,out_mask=out))
            need(names==target['names'],'retained named order')
            for number,(a,b) in enumerate(target['old_arc_differences']):
                ai,bi=names.index(a),names.index(b)
                U,V,W,p,q,case,inverse=compile_both(rows,ai,bi)
                helper=1000+1000*index+number;local=names+[helper]
                need(native.old.native(U,p)==V and native.old.native(V,q)==W,'both native path guards replay')
                after=tuple(x&((1<<len(rows))-1) for x in W[:-1])
                need(after==native.old.flip(rows,ai,bi),'single specified arc differs')
                moves.append(dict(changed_pair=[a,b],helper=helper,helper_out_mask=U[-1],paths=[[local[v] for v in path] for path in (p,q)],case=case,inverse=inverse,after_delete_sha256=digest(after)))
                rows=after
                gc.collect()
            need(rows==tuple(final),'compiled exact named endpoint')
        program=dict(parent_H=source['H'],target_H=target['H'],branch=target['branch'],permanent_births=new,changes=moves)
        receipt=replay(source,program,target);programs.append(program);receipts.append(receipt)
        report=dict(status='FINITE_EXACT_NATIVE_PROGRAMS' if len(programs)==14 else 'PARTIAL_NATIVE_PROGRAMS',scope='Every physical old name survives every branch. Endpoints strong; native intermediates need not preserve H, strongness or port stars. Complete through depth3 only; arbitrary-depth closure OPEN.',programs=programs,receipts=receipts,
          source_sha256={str(p.relative_to(PACKET)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'named_tree.json',PACKET/'inverse/native.py',PACKET/'inverse/native.json',Path(__file__))})
        (HERE/'native_tree.json').write_text(json.dumps(report,separators=(',',':'))+'\n')
        print('Compiled/replayed',target['H'],'native moves',receipt['native_moves'],'order',len(target['names']),flush=True)
    need(len(programs)==14,'Every finite-tree edge compiled')
if __name__=='__main__':main()
