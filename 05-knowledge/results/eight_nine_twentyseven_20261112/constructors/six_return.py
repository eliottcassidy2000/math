"""Internal-root preparation in the incoming marked six-state resource carrier."""
from pathlib import Path
from math import gcd
import json
from escape import child,context,require,root_state,strong

HERE=Path(__file__).resolve().parent

def carrier(phase):
    r,s=phase;a=(0,1);b=(2,3)
    order=a[:r]+b[:s]+a[r:]+b[s:]
    module=[0]*4
    for i,u in enumerate(order):
        for v in order[i+1:]:
            module[u]|=1<<v
    return context(tuple(module),(2,0),(False,True)),order,r+s

def move(phase):
    r,s=phase
    return 2-s,r+s-1

def run(phase,k):
    t,_,_=carrier(phase);root=0
    for _ in range(k):
        t,root=child(t,root)
    t,root=child(t,root,True)
    pre_seal=root_state(t,root)
    t,root=child(t,root)
    return t,root,pre_seal

def main():
    phase=(2,1);phases=[]
    for _ in range(6):
        t,order,cut=carrier(phase);z=root_state(t,0)
        require(z[0]==17 and strong(t),'carrier lost its full count or strongness')
        require(root_state(t,5)==(17,1,1,17),'outside-root neutrality failed')
        phases.append({'phase':phase,'order':order,'cut':cut,'rows':t,
                       'physical_root':0,'state_h_s_q_F':z})
        phase=move(phase)
    require(phase==(2,1),'marked resource state did not return')
    require(phases[0]['state_h_s_q_F']==(17,8,9,18),'wrong initial state')
    require(phases[2]['phase']==(0,2) and phases[2]['state_h_s_q_F']==(17,2,9,44),
            'wrong prepared state')
    returns=[];unprepared=[]
    for k in (1,2,3):
        t,root,before=run((0,2),k);z=root_state(t,root)
        b=105*(1<<k)+7;q=35*(1<<k)-11
        require(z==(b,q,q,b) and gcd(b,q)==1 and gcd(*before)==1 and strong(t),
                'prepared return failed')
        returns.append({'k':k,'rows':t,'root':root,'order':len(t),
                        'state_h_s_q_F':z,'pre_last_R_state':before})
        old,oldroot,oldbefore=run((2,1),k);oz=root_state(old,oldroot)
        oq=9*((1<<k)+1)
        require(oz==(3*oq,oq,oq,3*oq) and gcd(*oldbefore)==9,
                'unprepared hostile failed')
        unprepared.append({'k':k,'state_h_s_q_F':oz,'existing_pair_content':9})
    boundary,root,_=run((0,2),0)
    require(root_state(boundary,root)==(105,17,17,105),'wrong k=0 boundary')
    receipt={'status':'Exact controls; all-k proof in SIX_RETURN.md',
             'marked_resource_cycle':phases,'prepared_returns':returns,
             'unprepared_returns':unprepared,'k0_prepared_boundary':[105,17,17,105]}
    (HERE/'six_return_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'passed':True,'full_target':17,'prepared_defect':20,
                      'return_counts':[r['state_h_s_q_F'][0] for r in returns]}))

if __name__=='__main__':
    main()
