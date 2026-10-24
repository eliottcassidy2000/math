#!/usr/bin/env python3
"""Audit the richer native circuit from original matrices, full pin domains and Lz."""
from pathlib import Path
from itertools import combinations,permutations
from array import array
from collections import Counter
import csv,gzip,hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];REV='0f9dec9c3e1d'
CHECKS=0;PINS={}
def need(v,m):
    global CHECKS
    CHECKS+=1
    if not v:raise RuntimeError(m)
def read(path):
    raw=subprocess.check_output(['git','show',REV+':05-knowledge/results/'+path],cwd=ROOT)
    PINS[path]=hashlib.sha256(raw).hexdigest();return raw
def main():
    inventory=json.loads(read('session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json'))
    options={r['source']:r for r in inventory['records']}
    nodes=list(csv.DictReader(read('metagraph_order8_20261023/merged_nodes.tsv').decode().splitlines(),delimiter='\t'))
    owner=array('H');owner.frombytes(gzip.decompress(read('metagraph_order8_20261023/owner_map.u16.gz')))
    if sys.byteorder!='little':owner.byteswap()
    merged={int(n[k]):int(n['merged_id']) for n in nodes for k in ('class_a','class_b')}
    pot=json.loads(read('metagraph_order8_electrical_20261023/receipt.json'))['order8']
    Y=pot['radical_lifts_mod4'];X=pot['higher_layers'][-1]['lifts'][0]
    chars=[[(64*Y[j][i])%256 for i in range(3528)] for j in range(4)]+[X]
    plus=[176,1394,1457,1477,1532,1923,2087,2176,3369]
    minus=[1425,1549,1841,1921,2185,2553,2835,2905,3000]
    ss=plus+minus;sign=[1]*9+[-1]*9;U,V=1965,2234
    need(len(set(ss))==18 and U not in ss and V not in ss,'distinct sources outside targets')
    pairs=list(combinations(range(8),2));chords=[(i,j) for i,j in pairs if j>=i+2];idx={e:i for i,e in enumerate(chords)};full=(1<<21)-1
    def mirror(w):return sum(((w>>i)&1)<<idx[7-v,7-u] for i,(u,v) in enumerate(chords))
    def orbit(w):return {w,full^w,mirror(w),full^mirror(w)}
    def partner(B,p):
        E={frozenset(e) for e in zip(p,p[1:])}
        return [[0 if u==v else B[u][v] if frozenset((u,v)) in E else B[v][u] for v in range(8)] for u in range(8)]
    native=[];closed=set();selected=[];cost=0
    for s in ss:
        rec=options[s];code=int(nodes[s]['canon_a']);B=[[0]*8 for _ in range(8)]
        for bit,(u,v) in enumerate(pairs):B[u][v]=(code>>(27-bit))&1;B[v][u]=1-B[u][v]
        p,q=map(tuple,(rec['before']['path'],rec['after']['path']));fixed=[j for j in range(8) if p[j]==q[j]]
        shared=rec['pins'];start,end,j,c,d=shared
        need(p[0]==q[0]==start and p[-1]==q[-1]==end and p[j:j+2]==q[j:j+2]==(c,d),'actual shared pins and indexed arc')
        roles={start:'first',end:'last',c:'arc0',d:'arc1'}
        frame={i:roles.get(i,f'private:{s}:{i}') for i in range(8)}
        need(len(set(frame.values()))==8,'no owner identification within a slot')
        aut=[];paths=[]
        for pp in permutations(range(8)):
            if all(B[u][v]==B[pp[u]][pp[v]] for u,v in pairs):aut.append(pp)
            if all(pp[t]==p[t] for t in fixed) and all(B[u][v] for u,v in zip(pp,pp[1:])):paths.append(pp)
        need(aut==[tuple(range(8))],'full source automorphism audit; no dropped lifts')
        domain=[]
        for pp in paths:
            w=sum(B[pp[u]][pp[v]]<<bit for bit,(u,v) in enumerate(chords));ports=orbit(w)
            need(len(ports)==4 and not closed&ports,'distinct actual normalized resources')
            closed.update(ports);target=merged[owner[full^w]]
            need(merged[owner[w]]==s,'literal source owner')
            K=partner(B,pp);C=list(map(list,zip(*B)));rev=pp[::-1]
            need(partner(K,pp)==B,'literal complement involution')
            need(all(K[u][v] for u,v in zip(pp,pp[1:])),'target retains full path')
            need(partner(C,rev)==list(map(list,zip(*K))),'actual reflected converse square')
            need(sum(C[rev[u]][rev[v]]<<bit for bit,(u,v) in enumerate(chords))==mirror(w),'converse word gauge')
            need(all(pp[t]==p[t] for t in fixed),'all fixed owner clocks')
            domain.append({'path':pp,'target':target,'word':w,'closed_ports':sorted(ports)})
        eligible=[z for z in domain if z['target'] in (U,V)]
        need({tuple(z['path']) for z in eligible}=={p,q} and Counter(z['target'] for z in eligible)==Counter({U:1,V:1}),'full pinned domain after nonnegative target margins')
        need(len(domain)==(3 if s==1923 else 2),'raw pin-domain cardinality')
        for key,path in [('before',p),('after',q)]:
            row=next(z for z in eligible if tuple(z['path'])==path)
            need(row['word']==rec[key]['mask'] and row['target']==rec[key]['target'],'search witness literal authentication')
        need((X[s]-X[U])%8==(X[s]-X[V])%8,'frozen per-slot phase')
        label=[Y[t][s]%2 for t in range(4)]+[X[s]%128]
        need(label==rec['source_label'],'complete marked charge coordinates')
        cost+=8-len(fixed);selected.extend(eligible)
        native.append({'source':s,'canon':code,'fixed_clock_owners':[(t,p[t]) for t in fixed],
          'shared_roles':shared,'owner_frame':frame,'source_label':label,'domain':domain})
    need(cost==45,'exact changed-owner clock cost')
    need(len(closed)==148 and len(set(w for z in selected for w in z['closed_ports']))==144,'full and margin-compatible port counts')
    edges_raw=read('metagraph_order8_20261023/merged_edges.tsv.gz')
    edges=[tuple(int(row[k]) for k in ('source','target','blue','black')) for row in csv.DictReader(gzip.decompress(edges_raw).decode().splitlines(),delimiter='\t')]
    def lap(x):
        out=[0]*3528
        for i,j,b,k in edges:
            d=(b+k)*(x[i]-x[j]);out[i]+=d;out[j]-=d
        return out
    for chi in chars:
        need(all(z%256==0 for z in lap(chi)),'genuine full original harmonic separator')
        need(2*(chi[U]-chi[V])%256==0,'paired target charge equality')
    witnesses=json.loads((HERE/'image_witnesses.json').read_text())
    need(witnesses['input_sha256']==hashlib.sha256(edges_raw).hexdigest() and witnesses['modulus']==256,'exact image witness input')
    expected=[{s:2*t for s,t in zip(ss,sign)},{U:2,V:-2}]
    for sol,want in zip(witnesses['solutions'],expected):
        need(dict(sol['charge'])==want,'requested source and target charges')
        im=lap(sol['potential']);cor=dict(sol['integer_correction_nonzero'])
        need(all(z==want.get(i,0)+256*cor.get(i,0) for i,z in enumerate(im)),'every full integer Lz equation')
    need(len(witnesses['solutions'])==2,'complete image witness pair')
    count=[0]*(1<<18);mass=[0]*(1<<18);parity=[0]*(1<<18);balanced=0;low=[];accepted=[];reject=Counter()
    for bits in range(1,1<<18):
        bit=bits&-bits;j=bit.bit_length()-1;prev=bits^bit
        count[bits]=count[prev]+sign[j];mass[bits]=mass[prev]+sign[j]*X[ss[j]]
        parity[bits]=parity[prev]^sum((Y[t][ss[j]]%2)<<t for t in range(4))
    for bits in range(1<<18):
        if count[bits]:continue
        balanced+=1
        if parity[bits]==0 and mass[bits]%4==0:low.append(bits)
        if parity[bits]==0 and mass[bits]%128==0:accepted.append(bits)
        else:
            t=(parity[bits]&-parity[bits]).bit_length()-1 if parity[bits] else 4
            charge=2*sum(sign[j]*chars[t][ss[j]] for j in range(18) if bits>>j&1)
            need(charge%256!=0,'proper balanced subset has a nonzero target-owned current')
            reject[t]+=1
    need(accepted==[0,(1<<18)-1],'complete all-dyadic menu has two states')
    return {'status':'FINITE-EXACT PASS; native 18-slot primitive packet','slots':18,'plus':plus,'minus':minus,
      'integer_atomic_events':1,'owner_clock_cost':cost,'joint_UV_choices':1<<18,
      'target_balanced_choices':balanced,'H8_states':len(low),'H256_states':accepted,
      'raw_pinned_paths':sum(len(x['domain']) for x in native),'margin_compatible_paths':len(selected),
      'closed_endpoint_ports':72,'closed_union_ports':144,'raw_pinned_closed_ports':len(closed),
      'separators':dict(reject),'integer_X_debt':2*sum(t*X[s] for s,t in zip(ss,sign)),
      'literal_sources':native,'input_revision':REV,'input_sha256':PINS,'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
