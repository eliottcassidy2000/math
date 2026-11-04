"""Independent gap counts and named joint observers on all14 actual endpoints."""
from pathlib import Path
import sys,json,hashlib,time
HERE=Path(__file__).resolve().parent;PACKET=HERE.parent
sys.path.insert(0,str(PACKET/'interface'))
from gap import count,nine_fields,pair_fields,induced,need
from check import phase,lift,wrapper
from phase import observer17,QUERY_WORDS
sys.path.insert(0,str(PACKET/'inverse'))
import reader

def main():
    tree=json.loads((HERE/'named_tree.json').read_text());records=[]
    sources={tree['source']['H']:tree['source']};sources.update({z['H']:z for z in tree['targets']})
    for obj in tree['targets']:
        started=time.monotonic();names=obj['names'];rows=tuple(obj['rows']);loc={v:i for i,v in enumerate(names)}
        fb=tuple(loc[v] for v in obj['structural_feedback_names']);bn=obj['component_names'];B=tuple(obj['component_rows']);bw=bn.index(31)
        bf=tuple(bn.index(v) for v in obj['structural_feedback_names'] if v in bn)
        need(B==induced(rows,[loc[v] for v in bn]),'actual induced marked component')
        canonical_names=bn+obj['prefix']+[39]
        need(rows==induced(wrapper(B,bw,len(obj['prefix'])),[canonical_names.index(v) for v in names]),'actual whole relation equals named wrapper grammar')
        need(reader.strong(B) and reader.strong(rows),'independent strong component and endpoint')
        nine=nine_fields(B,bf,bw);h,s,q,F,e,P2,S2,E3,Q2=nine
        eight=[h,e,S2,F,Q2,q,P2,count(B,bf,3,starts=((0,bw),))]
        seed=obj['seed_reader'];seed_fields=reader.all_fields(tuple(seed['rows']))[seed['root']]
        need(eight==list(reader.run(seed_fields,seed['word'])),'literal component bins versus seed reader')
        if 'component_eight_fields' in obj:need(eight==obj['component_eight_fields'],'literal component versus inverse receipt')
        lifted=lift(nine,len(obj['prefix']));need(lifted['H']==obj['H']==3*obj['parent_H']+obj['offset'],'actual ternary return')
        pairs=[];named_charts={}
        for a,b in ((39,6),(31,6)):
            r,s=loc[a],loc[b];f=pair_fields(rows,fb,r,s)
            need(f[0]==obj['H'] and f[1]==lifted['P2'],'independent full graph H and P2')
            if a==39:need(f[2:7]==lifted['root_fields'],'independent full root profile')
            obs=observer17(rows,fb,(r,s));need(obs[0]==obj['H'],'actual universal17 baseline')
            pairs.append(dict(ports=[a,b],H=f[0],P2=f[1],port_data=[f[2:7],f[7:12]],ordered_W2=f[12:],observer17=obs))
            for v,p in ((a,f[2:7]),(b,f[7:12])):
                end1,start1,deletion,E2,start2=p;need((start2-end1)%2==0,'integer golden observer')
                g=[f[0],(start2-end1)//2,(start2+end1)//2]
                named_charts[str(v)]=dict(golden=g,charts=phase(g))
        need(named_charts['39']['golden']==lifted['golden'],'direct golden agrees with marked lift')
        parent=sources[obj['parent_H']]
        need(names[:len(parent['names'])]==parent['names'],'all parent physical names preserved in order')
        need(obj['prefix'][:2]==[50,51] and obj['root']==39 and obj['mark']==31 and obj['spare']==6,'fixed physical interfaces')
        records.append(dict(parent_H=obj['parent_H'],depth=obj['depth'],branch=obj['branch'],order=len(names),prefix=obj['prefix'],structural_feedback_names=obj['structural_feedback_names'],component_nine=nine,component_eight=eight,**lifted,root_charts=phase(lifted['golden']),named_port_charts=named_charts,joint_pairs=pairs,rows_sha256=hashlib.sha256(json.dumps(list(rows),separators=(',',':')).encode()).hexdigest()))
        print('Counted',obj['H'],'order',len(names),'feedback',len(fb),'seconds',round(time.monotonic()-started,2),flush=True)
    need(len(records)==14,'All14 named endpoints checked')
    pins=[HERE/'named_tree.json',HERE/'verify.py',PACKET/'interface/gap.py',PACKET/'interface/check.py',PACKET/'interface/phase.py',PACKET/'inverse/reader.py']
    report=dict(status='FINITE_EXACT_GAP_RECOUNTS_AND_NAMED_JOINT_OBSERVERS',scope='All14 endpoints through depth3 only. Independent literal gap counts verify reader/bin formulas and both ternary values. Named root/mark/spare phase data and two universal17 joint observers retained; arbitrary-depth closure OPEN.',query_words=QUERY_WORDS,records=records,source_sha256={str(p.relative_to(PACKET)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    (HERE/'verify.json').write_text(json.dumps(report,separators=(',',':'))+'\n')

if __name__=='__main__':main()
