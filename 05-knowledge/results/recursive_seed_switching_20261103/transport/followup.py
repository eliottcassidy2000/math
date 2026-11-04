"""Compile and independently count the two actual next graphs after13835447.

Saved programs are replayed on every run. Partial checkpoints are explicitly
labelled; only the final receipt includes both complete count/port records.
"""
from pathlib import Path
import gc,hashlib,json,sys,time
from itertools import combinations
HERE=Path(__file__).resolve().parent;PACKET=HERE.parent;RESULTS=PACKET.parent
OLD=RESULTS/'zero_ready_component_returns_20261103'
sys.path.insert(0,str(OLD/'continuation'))
from native_tree import digest,replay
import native
from reader import need,birth
sys.path.insert(0,str(OLD/'interface'))
from gap import count,nine_fields,pair_fields,induced
from check import phase,lift,wrapper
from phase import observer17,QUERY_WORDS,strong
import reader

def pair_path(rows,a,b):
    """A word with a,b adjacent, legal in rows or its single-pair reversal.

    Build a Hamilton word Q on the other vertices by first-outneighbor
    insertion. Before the first Q vertex beaten by either a or b, all Q
    vertices dominate both. Insert the pair there with a beating endpoint
    last. If no such Q vertex exists, append the pair. Only the pair's own
    orientation may disagree with rows. No search or physical relabelling
    is needed: select that orientation and apply the inherited gadget.
    """
    Q=[]
    for v in range(len(rows)):
        if v in (a,b):continue
        j=next((j for j,u in enumerate(Q) if rows[v]>>u&1),len(Q))
        Q.insert(j,v)
    for j,v in enumerate(Q):
        if rows[a]>>v&1 or rows[b]>>v&1:
            block=[a,b] if rows[b]>>v&1 else [b,a]
            return tuple(Q[:j]+block+Q[j:])
    return tuple(Q+[a,b])

def direct_compile(rows,a,b):
    order=pair_path(rows,a,b)
    x,y=(a,b) if order.index(a)<order.index(b) else (b,a)
    inverse=not bool(rows[x]>>y&1)
    oriented=native.old.flip(rows,a,b) if inverse else rows
    need(native.old.legal(oriented,order),'constructed one-orientation guard')
    U,V,W,p,q,case=native.old.traced_gadget(oriented,x,y,order)
    return (W,V,U,q,p,case,True,order) if inverse else (U,V,W,p,q,case,False,order)

def small_guard_controls():
    cases=[]
    for n in (2,3,4):
        pairs=list(combinations(range(n),2))
        for code in range(1<<len(pairs)):
            rows=[0]*n
            for j,(a,b) in enumerate(pairs):
                u,v=(a,b) if code>>j&1 else (b,a);rows[u]|=1<<v
            rows=tuple(rows)
            for a,b in pairs:
                U,V,W,p,q,case,inverse,order=direct_compile(rows,a,b)
                need(native.old.native(U,p)==V and native.old.native(V,q)==W,'small native guards')
                need(tuple(row&((1<<n)-1) for row in W[:n])==native.old.flip(rows,a,b),'small exact arc endpoint')
                cases.append([n,code,a,b,U[-1],list(p),list(q),list(order)])
    need(len(cases)==410,'all labelled tournaments/pairs through order4')
    return cases

def main():
    inp=PACKET/'selector/followup_named.json';data=json.loads(inp.read_text());source=data['source'];targets=data['targets'];saved=HERE/'followup.json'
    need(source==json.loads((PACKET/'selector/named_repair.json').read_text())['target'],'exact actual repaired parent reuse')
    need(len(targets)==2 and {z['offset'] for z in targets}=={2,4},'only two declared next children')
    input_hash=hashlib.sha256(inp.read_bytes()).hexdigest()
    cache=json.loads(saved.read_text()) if saved.exists() else {}
    completed={p['target_H']:p for p in cache.get('programs',[]) if p.get('guard_constructor')=='direct_pair_block_v1'} if cache.get('input_sha256')==input_hash else {}
    programs=[];native_receipts=[];records=[];guard_controls=small_guard_controls()
    def checkpoint():
        report=dict(status='FINITE_EXACT_NAMED_FOLLOWUP' if len(records)==2 else 'PARTIAL_NAMED_FOLLOWUP',
            scope='Only two specified next endpoints after the actual repaired13835447 source. All parent physical names survive; outer root39, active mark4, old mark31 and spare6 stay distinct. No intermediate H/strongness invariant or arbitrary-depth closure.',
            input_sha256=input_hash,parent_H=source['H'],source_names=source['names'],source_rows=source['rows'],
            programs=programs,native_receipts=native_receipts,records=records,query_words=QUERY_WORDS,
            direct_guard_control_columns=['n','tournament_code','pair_a','pair_b','helper_out','first_path','second_path','pair_guard'],
            direct_guard_controls=guard_controls)
        if len(records)==2:
            pins=[inp,PACKET/'selector/named_repair.json',Path(__file__),OLD/'continuation/native_tree.py',OLD/'inverse/native.py',
                RESULTS/'adaptive_metagraph_reformulation_20261103/metagraph_probe.py',OLD/'interface/gap.py',OLD/'interface/check.py',
                OLD/'interface/phase.py',OLD/'inverse/reader.py',RESULTS/'adaptive_state_induction_20261103/two_port_envelope.json',
                RESULTS/'adaptive_metagraph_reformulation_20261103/METAGRAPH.md',HERE/'FOLLOWUP.md']
            report['source_sha256']={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}
        saved.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    for obj in targets:
        if obj['H'] in completed:program=completed[obj['H']]
        else:
            names=list(source['names']);rows=tuple(source['rows']);loc={v:i for i,v in enumerate(obj['names'])};births=[];changes=[]
            for v in obj['new_names']:
                out=sum(1<<i for i,u in enumerate(names) if obj['rows'][loc[v]]>>loc[u]&1)
                rows=birth(rows,out);names.append(v);births.append(dict(name=v,out_mask=out))
            need(names==obj['names'],'newborn order retains all actual parent names')
            for j,(a,b) in enumerate(obj['old_arc_differences']):
                ai,bi=names.index(a),names.index(b);U,V,W,p,q,case,inverse,search_order=direct_compile(rows,ai,bi)
                need(native.old.native(U,p)==V and native.old.native(V,q)==W,'literal fresh native guards')
                after=tuple(x&((1<<len(names))-1) for x in W[:-1])
                need(after==native.old.flip(rows,ai,bi),'single requested arc change')
                helper=41000+100*obj['offset']+j;local=names+[helper]
                changes.append(dict(changed_pair=[a,b],helper=helper,helper_out_mask=U[-1],
                    paths=[[local[v] for v in word] for word in (p,q)],case=case,inverse=inverse,
                    pair_guard_names=[names[v] for v in search_order],after_delete_sha256=digest(after)))
                rows=after;gc.collect()
            need(rows==tuple(obj['rows']),'compiled exact full follow-up graph')
            program=dict(parent_H=source['H'],target_H=obj['H'],branch='2222'+str(obj['offset']),
                guard_constructor='direct_pair_block_v1',permanent_births=births,changes=changes)
        receipt=replay(source,program,obj);programs.append(program);native_receipts.append(receipt);checkpoint()
        print('Compiled/replayed',obj['H'],'moves',receipt['native_moves'],'order',len(obj['names']),flush=True)
    for obj in targets:
        names=obj['names'];rows=tuple(obj['rows']);loc={v:i for i,v in enumerate(names)};bn=obj['component_names'];B=tuple(obj['component_rows']);w=bn.index(4)
        feedback=obj['structural_feedback_names'];fb=tuple(loc[v] for v in feedback);bf=tuple(bn.index(v) for v in feedback if v in bn)
        need(B==induced(rows,[loc[v] for v in bn]),'actual induced component')
        canonical=bn+obj['prefix']+[39]
        need(rows==induced(wrapper(B,w,len(obj['prefix'])),[canonical.index(v) for v in names]),'literal named wrapper')
        need(strong(B) and strong(rows),'strong component and full endpoint')
        need(names[:len(source['names'])]==source['names'] and obj['prefix']==[50,51] and obj['mark']==4,'retained source names and declared active frame')
        seed=obj['seed_reader'];seed_fields=reader.all_fields(tuple(seed['rows']))[seed['root']]
        need(list(seed_fields)==seed['fields'],'independent seed DP')
        nine=nine_fields(B,bf,w);h,s,q,F,e,P2,S2,E3,Q2=nine
        eight=[h,e,S2,F,Q2,q,P2,count(B,bf,3,starts=((0,w),))]
        need(eight==obj['component_eight_fields']==list(reader.run(seed_fields,seed['word'])),'independent component bins versus reader')
        lifted=lift(nine,len(obj['prefix']));need(lifted['H']==obj['H']==3*source['H']+obj['offset'],'exact actual follow-up H')
        pairs=[];charts={}
        for a,b in ((39,6),(31,6),(4,6),(50,51)):
            f=pair_fields(rows,fb,loc[a],loc[b]);obs=observer17(rows,fb,(loc[a],loc[b]))
            need(f[0]==obj['H']==obs[0] and f[1]==lifted['P2'],'full H/P2 and universal17 baseline')
            if a==39:need(f[2:7]==lifted['root_fields'],'independent actual root profile')
            pairs.append(dict(ports=[a,b],H=f[0],P2=f[1],port_data=[f[2:7],f[7:12]],ordered_W2=f[12:],observer17=obs))
            for v,p in ((a,f[2:7]),(b,f[7:12])):
                end,_,_,_,start2=p;need((start2-end)%2==0,'integer rooted golden')
                golden=[f[0],(start2-end)//2,(start2+end)//2];charts[str(v)]=dict(golden=golden,charts=phase(golden))
            print('Counted',obj['H'],'pair',a,b,flush=True)
        need(charts['39']['golden']==lifted['golden'],'root golden versus lift')
        records.append(dict(parent_H=source['H'],offset=obj['offset'],names=names,rows=list(rows),feedback_names=feedback,
            root=39,active_mark=4,retained_old_mark=31,spare=6,prefix=obj['prefix'],new_names=obj['new_names'],
            component_nine=nine,component_eight=eight,**lifted,named_port_charts=charts,joint_pairs=pairs))
        checkpoint()
    need(len(records)==len(programs)==2,'both follow-up graphs complete')
if __name__=='__main__':main()
