"""Independent endpoint counts and joint ports for the seed-promotion repair."""
from pathlib import Path
import hashlib, json, sys, time

HERE=Path(__file__).resolve().parent; RESULTS=HERE.parent.parent
OLD=RESULTS/'zero_ready_component_returns_20261103'
sys.path.insert(0,str(OLD/'interface'))
from gap import count, nine_fields, pair_fields, induced, need
from check import phase, lift, wrapper
from phase import observer17, QUERY_WORDS, strong
sys.path.insert(0,str(OLD/'inverse'))
import reader


def main():
    inp=HERE.parent/'selector/named_repair.json'; data=json.loads(inp.read_text())
    source,obj=data['source'],data['target'];names=obj['names'];rows=tuple(obj['rows'])
    loc={v:i for i,v in enumerate(names)};bn=obj['component_names'];B=tuple(obj['component_rows'])
    feedback_names=[4,7,14,39]; fb=tuple(loc[v] for v in feedback_names)
    bf=tuple(bn.index(v) for v in feedback_names if v in bn);w=bn.index(obj['mark'])
    need(B==induced(rows,[loc[v] for v in bn]),'Actual induced component')
    canonical=bn+obj['prefix']+[39]
    need(rows==induced(wrapper(B,w,len(obj['prefix'])),[canonical.index(v) for v in names]),'Named wrapper relation')
    need(strong(B) and strong(rows),'Strong component and whole tournament')
    need(names==source['names'] and obj['prefix']==source['prefix'],'Old names and prefix')
    seed=obj['seed_reader']; seed_fields=reader.all_fields(tuple(seed['rows']))[seed['root']]
    need(list(seed_fields)==seed['fields'],'Independent small-seed DP')
    nine=nine_fields(B,bf,w);h,s,q,F,e,P2,S2,E3,Q2=nine
    eight=[h,e,S2,F,Q2,q,P2,count(B,bf,3,starts=((0,w),))]
    need(eight==obj['component_eight_fields']==list(reader.run(seed_fields,seed['word'])),'Independent component bins versus recurrence')
    lifted=lift(nine,len(obj['prefix']))
    need(lifted['H']==obj['H']==3*source['H']+2,'Actual fourth +2 count')
    pairs=[];charts={}
    for a,b in ((39,6),(31,6),(4,6),(50,51)):
        f=pair_fields(rows,fb,loc[a],loc[b]); obs=observer17(rows,fb,(loc[a],loc[b]))
        need(f[0]==obj['H']==obs[0] and f[1]==lifted['P2'],'Direct H/P2 and 17-query baseline')
        if a==39:need(f[2:7]==lifted['root_fields'],'Root lift')
        pairs.append(dict(ports=[a,b],H=f[0],P2=f[1],port_data=[f[2:7],f[7:12]],ordered_W2=f[12:],observer17=obs))
        for name,p in ((a,f[2:7]),(b,f[7:12])):
            end,start,deletion,E2,start2=p
            need((start2-end)%2==0,'Integer golden phase')
            golden=[f[0],(start2-end)//2,(start2+end)//2]
            charts[str(name)]=dict(golden=golden,charts=phase(golden))
        print('Verified named pair',a,b,flush=True)
    need(charts['39']['golden']==lifted['golden'],'Root golden lift')
    pins=[inp,Path(__file__),OLD/'interface/gap.py',OLD/'interface/check.py',OLD/'interface/phase.py',OLD/'inverse/reader.py']
    report=dict(status='FINITE-EXACT',scope='One fourth +2 return, all old names retained; old mark31 and new mark4 are distinct physical ports.',
        H=obj['H'],parent_H=source['H'],names=names,feedback_names=feedback_names,
        component_nine=nine,component_eight=eight,**{k:v for k,v in lifted.items() if k!='H'},
        named_port_charts=charts,joint_pairs=pairs,query_words=QUERY_WORDS,
        source_sha256={str(p.relative_to(RESULTS)):hashlib.sha256(p.read_bytes()).hexdigest() for p in pins})
    (HERE/'verify.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':main()
