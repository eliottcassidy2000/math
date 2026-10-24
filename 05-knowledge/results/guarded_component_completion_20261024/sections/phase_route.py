"""Literal endpoint admission, private-core swaps, and exact endpoint restoration."""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1];COUNT=0;U,V=1965,2234;FULL=(1<<21)-1
CHORDS=[(i,j)for i in range(8)for j in range(i+2,8)];REV=[CHORDS.index((7-j,7-i))for i,j in CHORDS]
def need(v,m):
    global COUNT
    COUNT+=1
    if not v:raise RuntimeError(m)
def add(a,b):return((((a>>4)+(b>>4))%128)<<4)|((a^b)&15)
def sub(a,b):return((((a>>4)-(b>>4))%128)<<4)|((a^b)&15)
def tile(m):
    r=sum(((m>>i)&1)<<j for i,j in enumerate(REV));return{m,r,FULL^m,FULL^r}
def main():
    paths=[RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json',
           RESULTS/'native_packet_decomposition_20261023/native/domains.json',HERE/'native_2125_receipt.json',
           HERE/'phase_sections_receipt.json',HERE.parent/'native/xor_tree_receipt.json',
           RESULTS/'metagraph_order8_electrical_20261023/receipt.json',HERE/'admission_receipt.json']
    supply,old,new,phase,host,electric,admission=[json.loads(p.read_text())for p in paths];old=old['sources'];coreids=sorted(z['source']for z in old)
    E=electric['order8'];X=E['higher_layers'][-1]['lifts'][0];Y=E['radical_lifts_mod4'];potentials=[X]+Y[:5];moduli=[256]+[4]*5
    records={z['source']:{r['target']:r for r in z['exact_all_shared_clocks_target_domain']}for z in old}
    for z in supply['helpers']:records[z['source']]={U:z['selected']['before'],V:z['selected']['after']}
    records[2125]={r['target']:r for r in new['records']};helpers=sorted(set(records)-set(coreids));ids=sorted(records)
    need(len(coreids)==28 and len(helpers)==54 and len(ids)==82,'exact whole native supply')
    labels={s:((X[s]%128)<<4)|sum((Y[j][s]%2)<<(j-1)for j in range(1,5))for s in ids}
    costs={s:sum(a!=b for a,b in zip(records[s][U]['path'],records[s][V]['path']))for s in ids}
    need(max(costs[s]for s in coreids)==4,'each old core replacement costs at most four clocks')
    def charge(ss):
        q=0
        for s in ss:q=add(q,labels[s])
        return q
    free=[s for s in helpers if s not in(368,776,3139)];opposite={}
    for ss in combinations(free,2):opposite.setdefault(charge([776,3139]+list(ss)),set(ss)|{776,3139})
    need(len(opposite)==972,'opposite forced phase has genuine missing boundary values')
    A={176,1394};initial={776,1596,2125,3139};total=add(charge(A),charge(initial))
    B=next(set(ss)for ss in combinations(coreids,2)if set(ss).isdisjoint(A)and sub(total,charge(ss))in opposite)
    final=opposite[sub(total,charge(B))]
    def hub(C):return{368}|set(phase['fixed_phase_section'][str(sub(total,charge(C)))])
    nearest={tuple(r['source']):set(r['target'])for r in admission['all_minimum_witnesses']}
    need(admission['worst_minimum_clock_cost']==49 and len(nearest)==1275,'exact complete source-specific admission bank')
    states=[(set(A),set(initial)),(set(A),nearest[tuple(sorted(initial))])];C=set(A)
    for a,b in zip(sorted(A-B),sorted(B-A)):
        C=(C-{a})|{b};states.append((C,nearest[tuple(sorted(final))]if C==B else hub(C)))
    states.append((set(B),set(final)))
    need(len(states)==5 and states[0][1]==initial and states[-1][1]==final,'two admissions and two private core exchanges')
    baseline=None;rows=[];actions=[]
    for C,H in states:
        S=C|H;need(len(C)==2 and len(H)==4 and add(charge(C),charge(H))==total,'exact core/helper margins and full charge')
        bits=[int(s in H)for s in(368,776,3139)]
        actual=next(t for t in host['actual_guard_truth_table']if t['bits']==bits)
        need(actual['allowed']and len(set(actual['actual_owners']))==1,'literal shared owner equations at arrival')
        words=set();current=[defaultdict(int)for _ in potentials]
        for s in ids:
            t=V if s in S else U;r=records[s][t];ww=tile(r['mask'])
            need(words.isdisjoint(ww)and len(ww)==4,'all source inverse converse resources injective at arrival');words|=ww
            for p,f in enumerate(potentials):
                d=2*(f[s]-f[t]);current[p][s]+=d;current[p][t]-=d
        need(len(words)==328,'exact actual82-unit arrival')
        if baseline is None:baseline=current
        for j,m in enumerate(moduli):need(all((current[j][v]-baseline[j][v])%m==0 for v in ids+[U,V]),'all physical generator rows including target degree terms')
        rows.append(dict(core_V=sorted(C),helper_V=sorted(H),shared_owner=actual['actual_owners'][0],
                         total_charge=total,helper_charge=charge(H),closed_word_count=len(words)))
        actions.append(current[0])
    events=[]
    for j,((a,h),(b,k))in enumerate(zip(states,states[1:])):
        changed=(a|h)^(b|k);support=len(changed);cost=sum(costs[s]for s in changed)
        need(support<=8 and cost<=51,'strong global coupled event bound')
        debt={v:actions[j+1][v]-actions[j][v]for v in ids+[U,V]if actions[j+1][v]!=actions[j][v]}
        need(all(d%256==0 for d in debt.values()),'actual integer-X debt paid modulo256')
        events.append(dict(changed=sorted(changed),support=support,clock_cost=cost,X_physical_row_delta=debt))
    out=dict(status='FINITE-EXACT literal global hub-route control',source_count=82,states=rows,events=events,
             bounds=dict(changed_slots=8,clock_occurrences=51,core_swap_clock=51,admission_clock=49),checks=COUNT,
             inputs={str(p.relative_to(RESULTS)):sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()for p in paths})
    out=json.loads(json.dumps(out));p=HERE/'phase_route_receipt.json'
    if '--write'in sys.argv:p.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',newline='\n')
    else:need(out==json.loads(p.read_text()),'frozen literal phase route')
    print(f'PASS {COUNT}; global8/51 bounds; literal event costs '+str([e['clock_cost']for e in events]))
if __name__=='__main__':main()
