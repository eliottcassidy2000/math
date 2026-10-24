#!/usr/bin/env python3
"""Independent enumeration and literal current-row audit of the charged helper route."""
from pathlib import Path
from itertools import combinations,permutations
from collections import defaultdict,deque
import json,subprocess,csv,gzip
from array import array
import sys
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def raw(p):return subprocess.check_output(['git','show','0f9dec9c3e1d:05-knowledge/results/'+p],cwd=ROOT)
def main():
    cert=json.loads((HERE/'route18.normal.json').read_text());H=cert['selected_helpers'];sources=cert['old_sources']+H
    need(H==[1638] and len(sources)==19,'nineteen literal sources')
    rows={r['source']:r for r in json.loads(raw('session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json'))['records']}
    pot=json.loads(raw('metagraph_order8_electrical_20261023/receipt.json'))['order8'];X=pot['higher_layers'][-1]['lifts'][0];Y=pot['radical_lifts_mod4']
    chars=[X]+[[64*x for x in y] for y in Y[:4]]
    globalindex={s:i for i,s in enumerate(cert['sources'])};old=cert['route'][0]
    wanted=[sum(chi[s] for s in sources if old>>globalindex[s]&1)%128 for chi in chars]
    states=[]
    # Complete 19 choose9 enumeration, with no translated-layer table or producer graph.
    for V in combinations(sources,9):
        if all(sum(chi[s] for s in V)%128==v for chi,v in zip(chars,wanted)):
            states.append(sum(1<<globalindex[s] for s in V))
    need(len(states)==42,'complete nineteen-source fibre')
    for b in (4,6):
        unseen=set(states);components=[]
        while unseen:
            todo=[unseen.pop()];seen=set(todo)
            for m in todo:
                neighbors=[v for v in unseen if (m^v).bit_count()<=b]
                unseen.difference_update(neighbors);seen.update(neighbors);todo.extend(neighbors)
            components.append(seen)
        row=next(r for r in cert['components'] if r['helpers']==H and r['budget']==b)
        need(sorted(map(len,components))==row['component_sizes'],'independent complete component census')
    nodes=list(csv.DictReader(raw('metagraph_order8_20261023/merged_nodes.tsv').decode().splitlines(),delimiter='\t'))
    own=array('H');own.frombytes(gzip.decompress(raw('metagraph_order8_20261023/owner_map.u16.gz')))
    if sys.byteorder!='little':own.byteswap()
    merged={int(n[k]):int(n['merged_id']) for n in nodes for k in ('class_a','class_b')}
    chords=[(i,j) for i,j in combinations(range(8),2) if j>i+1];index={e:i for i,e in enumerate(chords)};full=(1<<21)-1;options={};reservoir=set()
    for s in sources:
        code=int(nodes[s]['canon_a']);A=[[0]*8 for _ in range(8)]
        for k,(i,j) in enumerate(combinations(range(8),2)):A[i][j]=(code>>(27-k))&1;A[j][i]=1-A[i][j]
        p=rows[s]['before']['path'];q=rows[s]['after']['path'];fixed=[i for i in range(8) if p[i]==q[i]];eligible=[]
        for path in permutations(range(8)):
            if not all(path[i]==p[i] for i in fixed) or not all(A[i][j] for i,j in zip(path,path[1:])):continue
            w=sum(A[path[i]][path[j]]<<k for k,(i,j) in enumerate(chords));target=merged[own[full^w]]
            if target not in (1965,2234):
                need(s==1923 and target==1559,'only raw extra path, excluded by exact zero target count');continue
            need(merged[own[w]]==s,'actual source owner')
            rev=sum(((w>>k)&1)<<index[7-j,7-i] for k,(i,j) in enumerate(chords));ports={w,full^w,rev,full^rev}
            need(len(ports)==4 and not reservoir&ports,'all selected before/after normalized resources distinct');reservoir|=ports
            eligible.append(target);options[s,target]=ports
        need(sorted(eligible)==[1965,2234],'complete pinned U/V alternatives')
    route=cert['route'];need(set(route)<=set(states),'every arrival belongs to independently enumerated fibre')
    def packet_current(mask,chi):
        out=defaultdict(int)
        for s in sources:
            t=2234 if mask>>globalindex[s]&1 else 1965
            d=2*(chi[s]-chi[t]);out[s]+=d;out[t]-=d
        return out
    for m in route:
        used=set()
        for s in sources:
            t=2234 if m>>globalindex[s]&1 else 1965
            need(not used&options[s,t],'actual injection at every arrival');used|=options[s,t]
        for chi in chars:
            a=packet_current(old,chi);c=packet_current(m,chi)
            need(all((c[v]-a[v])%256==0 for v in set(a)|set(c)),'every original physical current row, diagonal included')
    supports=[(a^b).bit_count() for a,b in zip(route,route[1:])]
    need(max(supports)==6 and len(route)==5,'four six-slot events')
    need(cert['total_owner_clock_work']==57,'retained owner-clock work')
    return {'status':'PASS same-agent independent direct census and physical-row audit','sources':19,'states':42,
      'atomic_supports':supports,'events':4,'owner_clock_work':57,'reservoir_closed_ports':len(reservoir),'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
