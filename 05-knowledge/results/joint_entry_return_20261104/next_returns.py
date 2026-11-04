"""Literal two-child continuation from the inherited named H41506343 source.

Uses the inherited direct native compiler and exact gap counter. The compiled
counter keeps its n<=30 guard; longer three-command probes use Python integers.
"""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import importlib.util,json,sys,hashlib,os,shutil,subprocess,tempfile
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parent
OLD=RESULTS/'recursive_seed_switching_20261103'
sys.path.insert(0,str(OLD/'transport'))
import fast_gap
if os.name=='nt':
    def windows_executable():
        binary=Path(tempfile.gettempdir())/('math-gap-'+hashlib.sha256(fast_gap.SOURCE.read_bytes()).hexdigest()[:20]+'.exe')
        if not binary.exists():
            compiler=shutil.which('g++') or str(Path.home()/'scoop/apps/gcc/current/bin/g++.exe')
            subprocess.run([compiler,'-O3','-std=c++17',str(fast_gap.SOURCE),'-o',str(binary)],check=True)
        return binary
    fast_gap.executable=windows_executable
spec=importlib.util.spec_from_file_location('inherited_followup',OLD/'transport/followup.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
import gap,phase
need=base.need
python_count=gap.count
@lru_cache(None)
def count(rows,feedback,bins=1,starts=(),ends=()):
    return (fast_gap.count if len(rows)<=30 else python_count)(rows,feedback,bins,starts,ends)
gap.count=count;phase.count=count

def feedback(rows,names,seed_names,required):
    loc={v:i for i,v in enumerate(names)}
    others=[v for v in seed_names if v not in required]
    for k in range(len(others)+1):
        for extra in combinations(others,k):
            chosen=set(required)|set(extra);vs=[loc[v] for v in names if v not in chosen]
            chain=sorted(vs,key=lambda i:-sum(rows[i]>>j&1 for j in vs))
            if all(rows[u]>>v&1 for i,u in enumerate(chain) for v in chain[i+1:]):
                return [v for v in names if v in chosen]
    raise RuntimeError('No supplied feedback set')

def main():
    inp=OLD/'selector/followup_named.json'
    source=next(z for z in json.loads(inp.read_text())['targets'] if z['H']==41506343)
    probes=json.loads((HERE/'probe.json').read_text())['records'];records=[]
    for digit in (2,4):
        probe=next(z for z in probes if z['parent']==source['H'] and z['digit']==digit)
        # Prefer small cores, while preserving the whole actual parent's names.
        chosen=min(probe['hits'],key=lambda z:(len(z['seed']['rows']),z['word']))
        seed=chosen['seed'];k=chosen['prefix'];w=seed['root'];word=chosen['word']
        n=len(seed['rows'])+len(word)+k+1
        names=list(source['names']);new=list(range(max(names)+1,max(names)+1+n-len(names)));names+=new
        need(n>=len(source['names']),'No old-name deletion')
        reserved={39,4,6,31,50};prefix=[50]
        if k>=2:prefix.append(51);reserved.add(51)
        available=[v for v in names if v not in reserved]
        prefix+=available[:k-len(prefix)];available=available[k-1-int(k>=2):]
        seed_names=[None]*len(seed['rows']);seed_names[w]=4
        for v in (6,31):seed_names[seed_names.index(None)]=v
        for i in range(len(seed_names)):
            if seed_names[i] is None:seed_names[i]=available.pop(0)
        need(len(available)==len(word),'Named reader capacity')
        bn=seed_names+available;B=base.reader.graph(tuple(seed['rows']),w,word)
        canonical=bn+prefix+[39]
        rows=base.induced(base.wrapper(B,w,k),[canonical.index(v) for v in names])
        target=dict(H=probe['H'],rows=list(rows),names=names,root=39,mark=4,retained_old_mark=31,spare=6,
                    prefix=prefix,component_names=bn,component_rows=list(B),new_names=new,
                    seed_reader=dict(rows=seed['rows'],root=w,names=seed_names,word=word,fields=seed['fields']),
                    component_eight_fields=chosen['fields'])
        need(base.strong(rows) and base.strong(B),'Strong literal endpoint')
        current=tuple(source['rows']);current_names=list(source['names']);loc={v:i for i,v in enumerate(names)}
        births=[];changes=[]
        for v in new:
            out=sum(1<<i for i,u in enumerate(current_names) if rows[loc[v]]>>loc[u]&1)
            current=base.birth(current,out);current_names.append(v);births.append(dict(name=v,out_mask=out))
        for a,b in combinations(source['names'],2):
            ai,bi=loc[a],loc[b]
            if bool(current[ai]>>bi&1)==bool(rows[ai]>>bi&1):continue
            U,V,W,p,q,case,inverse,guard=base.direct_compile(current,ai,bi)
            need(base.native.old.native(U,p)==V and base.native.old.native(V,q)==W,'Both actual path guards')
            after=tuple(x&((1<<n)-1) for x in W[:-1])
            helper=51000+1000*digit+len(changes);labels=names+[helper]
            changes.append(dict(changed_pair=[a,b],helper=helper,helper_out_mask=U[-1],
                                paths=[[labels[v] for v in P] for P in (p,q)],case=case,inverse=inverse,
                                pair_guard=[names[v] for v in guard],after_delete_sha256=base.digest(after)))
            current=after
        need(current==rows,'Exact target after the native program')
        program=dict(parent_H=source['H'],target_H=probe['H'],permanent_births=births,changes=changes)
        receipt=base.replay(source,program,target)
        print('Native target',probe['H'],'order',n,'clock',k,'moves',receipt['native_moves'],flush=True)
        fb_names=feedback(rows,names,seed_names,{39,4});fb=tuple(loc[v] for v in fb_names)
        bf=tuple(bn.index(v) for v in fb_names if v in bn)
        seed_fields=base.reader.all_fields(tuple(seed['rows']))[w]
        need(list(seed_fields)==seed['fields'],'Independent seed DP')
        nine=gap.nine_fields(B,bf,w);h,s,q,F,e,P2,S2,E3,Q2=nine
        eight=[h,e,S2,F,Q2,q,P2,count(B,bf,3,starts=((0,w),))]
        need(eight==chosen['fields']==list(base.reader.run(seed_fields,word)),'Independent full component versus reader')
        lifted=base.lift(nine,k);need(lifted['H']==target['H']==3*source['H']+digit,'Requested ternary scalar')
        need(python_count(rows,fb)==target['H'],'Independent Python full Hamilton count')
        joint=[];charts={}
        for a,b in ((39,6),(31,6),(4,6),(50,51)):
            values=gap.pair_fields(rows,fb,loc[a],loc[b]);obs=phase.observer17(rows,fb,(loc[a],loc[b]))
            need(values[0]==obs[0]==target['H'] and values[1]==lifted['P2'],'Actual count and joint observer agree')
            for v,p in ((a,values[2:7]),(b,values[7:12])):
                end,_,_,_,start2=p;need((start2-end)%2==0,'Integral phase coordinate')
                golden=[values[0],(start2-end)//2,(start2+end)//2]
                charts[str(v)]=dict(golden=golden,charts=base.phase(golden))
            joint.append(dict(ports=[a,b],fields=values,observer17=obs))
            print('Joint target',probe['H'],'ports',a,b,flush=True)
        need(charts['39']['golden']==lifted['golden'],'Full root phase agrees with lift')
        records.append(dict(target=target,program=program,receipt=receipt,feedback_names=fb_names,
                            component_nine=nine,lifted=lifted,joint_pairs=joint,named_port_charts=charts))
        status='FINITE_EXACT_TWO_CHILD_NATIVE_RETURNS' if len(records)==2 else 'PARTIAL'
        pins=[inp,HERE/'probe.json',Path(__file__),OLD/'transport/followup.py',OLD/'transport/fast_gap.py',OLD/'transport/fast_gap.cpp']
        out=dict(status=status,scope='Only the two declared children of this named H41506343 source. All old names retained; the +2 prefix clock drops to one and physical51 joins the component. No arbitrary-depth return law.',source=source,records=records,
                 source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
        (HERE/'next_returns.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
