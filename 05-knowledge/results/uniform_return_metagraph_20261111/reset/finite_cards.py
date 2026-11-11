"""Optimal initial one-deletion dip at141, and hostile prefix lifts."""
from pathlib import Path
from itertools import combinations
import json
from classify_modules import H,strong,need,inherited as a,two_birth as b,frontier
from two_clock import multiruns,exponential
HERE=Path(__file__).resolve().parent

def main():
    source=a.ready_source((2,3,1),2)
    cards=[{'deleted':v,'H':H(a.induced(source,[w for w in range(len(source)) if w!=v]))}
           for v in range(len(source))]
    need(max(x['H'] for x in cards)==73,'complete one-deletion maximum')
    kept=[0,2,3,4,5,6,7];parent=a.induced(source,kept)
    need(parent==a.ready_source((2,3,1),1),'prefix deletion is exact smaller source')
    witnesses=[]
    for out in (53,117):
        child=a.attach(parent,out)
        need(H(child)==425 and strong(child),'actual one-birth target')
        need(a.induced(child,range(7))==parent,'all seven old vertices and arcs retained')
        co=b.module_runs(child,0)
        need(co[3]==116 and not any(co[4:]),'higher growth channel retained')
        witnesses.append({'source_rows':source,'kept_old_labels':kept,'parent_rows':parent,
                          'new_outmask':out,'target_rows':child,'trace':[141,73,425],
                          'module_coefficients':co})
    prefix_probes=[]
    for k in range(2,9):
        parent=a.ready_source((2,3,1),k-1);old_h=34*(1<<k)+5;values=b.all_birth_counts(parent)
        hits={str(offset):[i for i,v in enumerate(values) if v==3*old_h+offset] for offset in (2,4)}
        for offset,outs in hits.items():
            for out in outs:need(H(a.attach(parent,out))==3*old_h+int(offset),'prefix probe literal hit')
        prefix_probes.append({'original_k':k,'parent_rows':parent,'cuts':len(values),'hits':hits})
    split=[];parent=a.ready_source((2,3,1),3)
    for out in (106,482):
        target=a.attach(parent,out);need(H(target)==1651,'positive finite+4 target')
        singles=[b.module_runs(target,v) for v in range(3)]
        need(not any(co[2]==816 and not any(co[3:]) for co in singles),'no one-block uniform+4 lift')
        pairs=[]
        for modules in combinations(range(3),2):
            exp=exponential(multiruns(target,modules))
            need(exp!={1:19,4:408},'no two-block uniform+4 lift')
            pairs.append({'modules':modules,'exponential':exp})
        need(pairs[-1]['exponential'].get(9,0)>0,'three-block lift contains faster nine-mode subfamily')
        split.append({'outmask':out,'rows':target,'single_runs':singles,'two_module_responses':pairs})
    collar=(240,5,9,3,238,206,142,14)
    need(H(collar)==141 and strong(collar),'minimal-defect collar source')
    vals=b.all_birth_counts(collar)
    need(425 not in vals and 427 not in vals,'collar one-birth boundary')
    collar_layers=frontier(collar,427)
    need([z['count'] for z in collar_layers]==[27,4,0],'complete collar strong frontier')
    need(not any(c['h'] in(425,427) for z in collar_layers for c in z['children']),'collar all-order target exclusion')
    collar_card_hits=[]
    for drop in range(len(collar)):
        ids=[v for v in range(len(collar)) if v!=drop];W=a.induced(collar,ids)
        cc=b.all_birth_counts(W)
        for out,h in enumerate(cc):
            if h not in(425,427):continue
            child=a.attach(W,out);need(H(child)==h,'collar card literal hit')
            collar_card_hits.append({'deleted':drop,'card_H':H(W),'kept_labels':ids,
                                     'outmask':out,'target_H':h,'target_rows':child})
    report={'status':'PASS','source_cards':cards,'optimal_card_H':73,'witnesses':witnesses,
            'prefix_cut_universes':prefix_probes,'split_lift_hostiles':split,
            'collar_rows':collar,'collar_frontier427':collar_layers,'collar_card_hits':collar_card_hits,
            'collar_scope':'Complete strong frontier under427, promoted by the low-ratio ear theorem and proper-divisor test.'}
    (HERE/'finite_cards.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print('PASS: optimal141->73->425 initial deletion, finite cut probes and split-clock hostiles')
if __name__=='__main__':main()
