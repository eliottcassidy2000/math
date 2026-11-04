"""Independent exact named-interface counts for the parent's actual targets."""
from pathlib import Path
import hashlib,json,time
from gap import count,nine_fields,pair_fields,induced,need
from check import phase,lift,birth,wrapper
HERE=Path(__file__).resolve().parent
QUERY_WORDS=('', '1','2','3','01','11','21','31','02','12','22','32','201','301','021','221','321')

def strong(rows):
    for r in range(len(rows)):
        reached=1<<r
        while True:
            expanded=reached
            for v in range(len(rows)):
                if reached>>v&1:expanded|=rows[v]
            if expanded==reached:break
            reached=expanded
        if reached!=(1<<len(rows))-1:return False
    return True

def observer17(rows,feedback,pair):
    out=[];feedback=tuple(sorted(set(feedback)|set(pair)))
    for word in QUERY_WORDS:
        T=rows
        for letter in word:
            mask=int(letter);incoming=sum(1<<pair[i] for i in range(2) if mask>>i&1)
            T=birth(T,incoming)
        out.append(count(T,feedback))
    return out

def main():
    source=HERE.parent/'inverse/named_targets.json';data=json.loads(source.read_text());records=[]
    old_source=HERE.parent.parent/'triangle_return_170807_20261103/zero_ready.json'
    old=json.loads(old_source.read_text())[1]
    old_names=old['suffix_names']+old['prefix_names']+[old['root']]
    old_rows=wrapper(tuple(old['suffix_rows']),old['suffix_names'].index(old['inport']),len(old['prefix_names']))
    need(tuple(data['source_rows'])==induced(old_rows,[old_names.index(v) for v in data['source_names']]),
         'actual named 170807 source equals inherited marked-component wrapper')
    for obj in data['targets']:
        names=obj['names'];rows=tuple(obj['rows']);fb=tuple(names.index(v) for v in obj['structural_feedback_names'])
        bn=obj['component_names'];B=tuple(obj['component_rows']);bw=bn.index(31)
        need(B==induced(rows,[names.index(v) for v in bn]),'literal named component induced in target')
        wrapped_names=bn+obj['prefix']+[obj['root']]
        wrapped=wrapper(B,bw,len(obj['prefix']))
        need(rows==induced(wrapped,[wrapped_names.index(v) for v in names]),'exact full target zero-ready wrapper')
        need(strong(B) and strong(rows),'strong target component and full wrapper')
        bf=tuple(bn.index(v) for v in obj['structural_feedback_names'] if v in bn)
        nine=nine_fields(B,bf,bw);h,s,q,F,e,P2,S2,E3,Q2=nine
        eight=[h,e,S2,F,Q2,q,P2,count(B,bf,3,starts=((0,bw),))]
        lifted=lift(nine,2);need(lifted['H']==obj['H']==3*170807+obj['offset'],'actual target exact H')
        pairs=[]
        for a,b in ((39,6),(31,6),(50,51),(39,31)):
            r=names.index(a);s=names.index(b);fields=pair_fields(rows,fb,r,s)
            need(fields[0]==obj['H'] and fields[1]==lifted['P2'],'full named target independent gap H,P2')
            if a==39:need(fields[2:7]==lifted['root_fields'],'full root phase versus polynomial lift')
            x=observer17(rows,fb,(r,s));need(x[0]==fields[0],'observer17 first coordinate H')
            pairs.append(dict(ports=[a,b],H=fields[0],P2=fields[1],port_data=[fields[2:7],fields[7:12]],
                              ordered_W2=fields[12:],observer17=x))
        records.append(dict(offset=obj['offset'],names=names,rows_sha256=hashlib.sha256(json.dumps(list(rows)).encode()).hexdigest(),
            structural_feedback_names=obj['structural_feedback_names'],component_nine=nine,component_eight=eight,
            **lifted,charts=phase(lifted['golden']),joint_pairs=pairs))
    sources=[source,old_source,HERE/'phase.py',HERE/'gap.py',HERE/'check.py',
             HERE.parent.parent/'adaptive_state_induction_20261103/two_port_envelope.json']
    print(json.dumps(dict(scope='FINITE-EXACT named graphs; four joint pairs and universal17 continuation queries per branch',
        query_words=QUERY_WORDS,records=records,source_sha256={str(p.relative_to(HERE.parent.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}),indent=2,sort_keys=True))
if __name__=='__main__':main()
