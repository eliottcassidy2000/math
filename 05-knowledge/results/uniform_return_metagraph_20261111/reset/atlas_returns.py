"""Complete one-tail/two-sink/apex dictionary on the sixteen neutral4 rays."""
from pathlib import Path
from itertools import permutations
import json
from classify_modules import H,strong,need,inherited as a,two_birth as b,converse
HERE=Path(__file__).resolve().parent
ATLAS=[('00 01',9,-1),('00 10',11,1),('01 01',9,-3),('01 11',15,-1),
       ('10 10',12,1),('10 11',16,1),('01 00 10',25,5),('01 01 10',34,5),
       ('01 10 01',30,5),('01 10 10',30,7),('01 11 10',34,5),('10 00 01',22,5),
       ('10 01 01',25,5),('10 01 10',33,9),('10 10 01',33,7),('10 11 01',31,5)]


def construction(tail,drop,out,k):
    source=a.ready_source(tail,k)
    kept=[v for v in range(len(source)) if v!=k+drop]
    parent=a.induced(source,kept)
    raw=a.attach(a.attach(parent,0),0)
    # Expand the one prefix point to a homogeneous transitive k-module.
    lifted=((out>>1)<<k)|(((1<<k)-1) if out&1 else 0)
    raw=a.attach(raw,lifted)
    sinks=[v for v,row in enumerate(raw) if row==0]
    need(all(v>=len(parent) for v in sinks),'only newly introduced neutral padding removed')
    final=a.induced(raw,[v for v in range(len(raw)) if v not in sinks])
    return source,kept,parent,raw,final,sinks


def cycle(rows):
    for rest in permutations(range(1,len(rows))):
        word=(0,)+rest
        if all(rows[x]>>y&1 for x,y in zip(word,word[1:]+word[:1])):return word
    raise RuntimeError('missing Hamilton cycle')


def kernel(rows,root):
    ports=[v for v in range(len(rows)) if rows[v]>>root&1]
    ids=[v for v in range(len(rows)) if v not in(0,root)]
    W=a.induced(rows,ids);J=[ids.index(v) for v in ports]
    hp,en=a.dp(W);_,st=a.dp(converse(W));full=len(hp)-1
    E=sum(sum(en[S][r] for r in J)*hp[full^S] for S in range(1,full+1))
    C=sum(en[S][i]*st[full^S][j] for S in range(1,full) for i in J for j in J if i!=j)
    residual=a.induced(rows,[v for v in range(len(rows)) if v not in ports+[root]])
    return {'rows':W,'ports':J,'E':E,'q':hp[-1],'C':C,'residual_rows':residual,'residual_H':H(residual)}


def main():
    report=[];totalcuts=candidates=hits=0
    for word,A,B in ATLAS:
        tail=tuple({'00':0,'01':2,'10':1,'11':3}[s] for s in word.split())
        source=a.ready_source(tail,1);p=a.root_profiles(source)[-1]
        need(p['h']==2*A+B and p['gap']==4 and p['q']==p['s'],'inherited complete neutral4 atlas controls')
        row={'tail':word,'A':A,'B':B,'base_scalar_candidates':[],'uniform_operators':[]}
        for drop in range(len(tail)):
            kept=[v for v in range(len(source)) if v!=drop+1]
            parent=a.induced(source,kept);prepared=a.attach(a.attach(parent,0),0)
            vals=b.all_birth_counts(prepared);totalcuts+=len(vals)
            for out,h in enumerate(vals):
                offsets=[c for c in(2,4) if h==3*p['h']+c]
                if not offsets:continue
                offset=offsets[0];candidates+=1;raw=a.attach(prepared,out);co=b.module_runs(raw)
                need(H(raw)==h,'independent scalar candidate count')
                good=co[2]==3*A and not any(co[3:])
                row['base_scalar_candidates'].append({'delete_tail_index':drop,'outmask':out,
                    'offset':offset,'runs':co,'uniform':good})
                if not good:continue
                hits+=1;_,_,_,raw,final,sinks=construction(tail,drop,out,1)
                need(strong(final) and H(final)==h,'strong simplified target')
                root=len(parent)-1;K=kernel(final,root)
                need(K['E']==3*A and K['q']-K['C']==3*B+offset,'exact target kernel identity')
                controls=[]
                for k in range(1,7):
                    src,ids,W,U,V,removed=construction(tail,drop,out,k);r=len(W)-1
                    need(H(src)==A*(1<<k)+B and H(U)==3*H(src)+offset and H(V)==H(U),'literal all-k target')
                    need(a.induced(V,range(len(W)))==a.induced(src,ids),'every retained old arc and root/ports')
                    need(strong(V),'all literal outputs strong')
                    prof=a.root_profiles(V)[r]
                    need((prof['s'],prof['q'],prof['gap'],prof['delta'])==(K['q'],K['q'],K['C'],K['C']),
                         'literal target neutral profile')
                    need(V==b.inflate(final,0,k),'fixed pointed target template')
                    controls.append({'k':k,'kept_old_labels':ids,'removed_new_sinks':removed,
                                     'parent_H':H(W),'target_rows':V,'profile':prof})
                row['uniform_operators'].append({'delete_tail_index':drop,'outmask_k1_two_sink':out,
                    'offset':offset,'raw_strong':strong(raw),'removed_new_sink_labels_k1':sinks,
                    'base_rows':final,'base_hamilton_cycle':cycle(final),'kernel':K,'literal_controls':controls})
        report.append(row)
    need((totalcuts,candidates,hits)==(9216,77,6),'complete operator-universe totals')
    ready=json.loads((HERE/'ready425_inverse.json').read_text())
    need(ready['status']=='PASS' and ready['limit']>=199 and
         all(199 not in vals for vals in ready['lane_values'].values()),
         'complete inherited ready-image hostile at first+4 target199')
    data={'status':'PASS','scope':'All16 neutral4 tails, one old-tail deletion, two global sinks, all homogeneous-prefix apex cuts',
          'literal_final_cuts':totalcuts,'scalar_candidates':candidates,'uniform_operators':hits,'atlas':report,
          'first_plus4_target199_ready_lanes':[]}
    (HERE/'atlas_returns.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print('PASS: 9216 cuts,77 scalar candidates,6 uniform instructions across5 source families; actual+4 ray')


if __name__=='__main__':main()
