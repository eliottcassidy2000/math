#!/usr/bin/env python3
"""Original order-eight electrical phases under literal contiguous path edits.

These are declared observer-block moves, not asserted old metagraph edges.
All nonempty endpoint/phase fibres of the four frozen sources are included.
"""
from array import array
from collections import Counter, defaultdict, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import csv, gzip, json, sys

ROOT=Path(__file__).resolve().parents[3]
RESULTS=ROOT/'05-knowledge/results'
CHECKS=0


def need(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(label)


def components(adj):
    unseen=set(adj);out=[]
    while unseen:
        todo=[min(unseen)];unseen.remove(todo[0]);part=[]
        while todo:
            u=todo.pop();part.append(u)
            for v in adj[u]:
                if v in unseen:unseen.remove(v);todo.append(v)
        out.append(sorted(part))
    return out


def block_chain(paths,b):
    """Retain actual owners outside each internal b-clock block."""
    N=len(paths);windows=7-b;P=defaultdict(F)
    for start in range(1,8-b):
        groups=defaultdict(list)
        for i,p in enumerate(paths):groups[p[:start]+p[start+b:]].append(i)
        for group in groups.values():
            for i in group:
                for j in group:P[i,j]+=F(1,windows*len(group))
    adj={i:[] for i in range(N)}
    for (i,j),value in P.items():
        if i!=j and value:adj[i].append(j)
    return P,adj


def full_word_matrix(canon):
    A=[[0]*8 for _ in range(8)]
    for bit,(i,j) in enumerate(combinations(range(8),2)):
        A[i][j]=(canon>>(27-bit))&1;A[j][i]=1-A[i][j]
    return A


def route_certificate(paths,P,adj):
    N=len(paths);load=Counter();maxlen=0;allroutes={}
    for x in range(N):
        previous={x:None};todo=deque([x])
        while todo:
            a=todo.popleft()
            for b in sorted(adj[a]):
                if b not in previous:previous[b]=a;todo.append(b)
        need(len(previous)==N,'complete legal-move connectivity')
        for y in range(N):
            route=[];u=y
            while u is not None:route.append(u);u=previous[u]
            route.reverse();length=len(route)-1;maxlen=max(maxlen,length)
            allroutes[x,y]=route
            for a,b in zip(route,route[1:]):load[a,b]+=length
    congestion=max((F(c,N)/P[e] for e,c in load.items()),default=F(0))
    return congestion,maxlen,allroutes


def main():
    parent=json.loads((RESULTS/'native_guard_moves_farey_20261023/guards/receipt.json').read_text())['original_eight']
    potfile=RESULTS/'metagraph_order8_electrical_20261023/receipt.json'
    potraw=potfile.read_bytes()
    normalized=potraw.replace(b'\r\n',b'\n')
    need(sha256(normalized).hexdigest()=='417e4cc3af2688df39bbb71701789a6967677ae82ea8d0f0badb8ff50d739b27',
         'committed LF potential bytes')
    legacy=sha256(normalized.replace(b'\n',b'\r\n')).hexdigest()
    need(legacy=='d8ef442adbe0261bf8ae9892db925e220bd35bd67546d2933d1ab7b2f1c8f0d4',
         'legacy receipt pins exactly the CRLF view of the same bytes')
    need(parent['original_potential_sha256'] in (sha256(normalized).hexdigest(),legacy),
         'inherited metadata is either the recorded LF or CRLF input view')
    x=json.loads(potraw)['order8']['higher_layers'][-1]['lifts'][0]
    folder=RESULTS/'metagraph_order8_20261023'
    nodes=list(csv.DictReader((folder/'merged_nodes.tsv').open(),delimiter='\t'))
    mid=[None]*6880
    for row in nodes:
        for field in ('class_a','class_b'):mid[int(row[field])]=int(row['merged_id'])
    raw=gzip.decompress((folder/'owner_map.u16.gz').read_bytes())
    need(sha256(raw).hexdigest()==parent['owner_bytes_sha256'],'original owner-map gauge bytes')
    owner=array('H');owner.frombytes(raw)
    if sys.byteorder!='little':owner.byteswap()
    chords=[(i,j) for i in range(8) for j in range(i+2,8)];full=(1<<21)-1
    reports=[];witnesses=[]
    for source,expected in zip((1238,2242,2246,2315),(409,407,523,527)):
        ci=int(nodes[source]['class_a']);A=full_word_matrix(int(nodes[source]['canon_a']))
        paths=[p for p in permutations(range(8)) if all(A[u][v] for u,v in zip(p,p[1:]))]
        need(len(paths)==expected,'all actual source Hamiltonian paths')
        values={};targets={}
        for p in paths:
            mask=sum(A[p[i]][p[j]]<<b for b,(i,j) in enumerate(chords))
            need(owner[mask]==ci,'actual fixed source with original forward-bit gauge')
            target=mid[owner[full^mask]];targets[p]=target;values[p]=(x[source]-x[target])%8
            reverse=p[::-1]
            cmask=sum(A[reverse[j]][reverse[i]]<<b for b,(i,j) in enumerate(chords))
            need(mid[owner[cmask]]==source and mid[owner[full^cmask]]==target,'converse transports both merged owners')
        for keep_target in (False,True):
            fibres=defaultdict(list)
            for p in paths:
                key=(p[0],p[-1],values[p])+( (targets[p],) if keep_target else () )
                fibres[key].append(p)
            histogram=Counter();nontrivial=0;worst=F(0);worstlen=0;chosen=[];small_witness=None
            for key,pp in sorted(fibres.items()):
                if len(pp)==1:histogram[0]+=1;continue
                nontrivial+=1
                for b in range(2,7):
                    P,adj=block_chain(pp,b);cc=components(adj)
                    if len(cc)==1:break
                histogram[b]+=1
                c,d,rr=route_certificate(pp,P,adj);worst=max(worst,c);worstlen=max(worstlen,d)
                for i in adj:
                    for j in adj[i]:
                        g=dict(zip(pp[i],pp[j]))
                        automorphism=all(A[u][v]==A[g[u]][g[v]] for u,v in combinations(range(8),2))
                        if automorphism:continue
                        width=max(k for k in range(8) if pp[i][k]!=pp[j][k])-min(k for k in range(8) if pp[i][k]!=pp[j][k])+1
                        if small_witness is None or width<small_witness['block_width']:
                            small_witness=dict(source=source,keep_target=keep_target,nonautomorphism_move=True,
                                               block_width=width,actual_paths=[pp[i],pp[j]],
                                               target_owners=[targets[pp[i]],targets[pp[j]]],
                                               original_phase=values[pp[i]],
                                               fixed_owner_clocks={0:pp[i][0],7:pp[i][-1]})
                # Keep one strongest endpoint certificate for each source/filter.
                if not chosen or (b,len(pp))>(chosen[0]['minimum_block'],chosen[0]['states']):
                    best=max(rr,key=lambda pair:len(rr[pair]));route=rr[best]
                    chosen=[dict(source=source,keep_target=keep_target,endpoint_guard=key,
                                 minimum_block=b,states=len(pp),diameter=d,congestion=str(c),
                                 actual_paths=[pp[i] for i in route],
                                 target_owners=[targets[pp[i]] for i in route],
                                 original_phases=[values[pp[i]] for i in route],
                                 below_budget_components=components(block_chain(pp,b-1)[1]) if b>2 else [])]
                # Record a genuine same-phase, same-ports change of target owner.
                if not keep_target and len(set(targets[p] for p in pp))>1:
                    candidates=[(i,j) for i in adj for j in adj[i] if targets[pp[i]]!=targets[pp[j]]]
                    if candidates and not any(v['source']==source and v.get('cross_target') for v in witnesses):
                        i,j=min(candidates,key=lambda pair:sum(a!=b for a,b in zip(pp[pair[0]],pp[pair[1]])))
                        witnesses.append(dict(source=source,cross_target=True,minimum_fibre_block=b,
                                              actual_paths=[pp[i],pp[j]],target_owners=[targets[pp[i]],targets[pp[j]]],
                                              original_phase=values[pp[i]],fixed_owner_clocks={0:pp[i][0],7:pp[i][-1]}))
            witnesses.extend(chosen)
            if small_witness:witnesses.append(small_witness)
            reports.append(dict(source=source,keep_target=keep_target,nonempty_fibres=len(fibres),
                                nontrivial_fibres=nontrivial,minimum_block_histogram=dict(sorted(histogram.items())),
                                worst_canonical_congestion=str(worst),maximum_shortest_route=worstlen))
    return dict(status='FINITE-EXACT PASS; declared observer-block moves only',reports=reports,
                witnesses=witnesses,checks=CHECKS,
                source_pins={'potential_LF_sha256':sha256(normalized).hexdigest(),
                             'legacy_CRLF_sha256':legacy,
                             'owner_bytes_sha256':sha256(raw).hexdigest()})


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
