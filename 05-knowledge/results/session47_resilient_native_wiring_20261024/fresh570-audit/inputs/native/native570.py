"""Fresh common-clock sources reach570 without changing any old28 pin."""
from common import *
from itertools import permutations,combinations
from collections import defaultdict,Counter
from shelter import partner
CHECKS=0
def check(v,m):
    global CHECKS
    CHECKS+=1;need(v,m)
def main():
    nodes,owner,merged,X,Y=data();old=inventory();S=[1130,1673,1754,2678];A={1130,1673};B={1754,2678};T=(570,2234)
    candidates=json.loads((HERE/'contexts570.json').read_text())['pairs'];domain=[];reservoir=set();raw_counts=[]
    for r in old:
        for key in ('before','after'):
            O=set(orbit(r[key]['mask']));check(len(O)==4 and not O&reservoir,'old alternative reservoir');reservoir.update(O)
    old_reservoir=set(reservoir)
    for s in S:
        options=[r for r in candidates if r['source']==s and r['target']==2234 and r['fixed']]
        r=min(options,key=lambda r:(r['cost'],r['before_path'],r['after_path']));p=tuple(r['before_path']);q=tuple(r['after_path']);F=r['fixed']
        M=matrix(nodes[s]['canon_a']);raw=[];caps=Counter();aut=[]
        for h in permutations(range(8)):
            if all(M[i][j]==M[h[i]][h[j]] for i,j in PAIRS):aut.append(h)
            if not all(M[i][j] for i,j in zip(h,h[1:])):continue
            w=word(M,h);t=merged[owner[FULL^w]];caps[t]+=1
            if all(h[i]==p[i] for i in F):raw.append((h,w,t))
        check(aut==[tuple(range(8))],'fresh literal source has no hidden automorphism multiplicity')
        chosen=[r for r in raw if r[2] in T];check({h for h,w,t in chosen}=={p,q} and len(chosen)==2,'complete fresh binary domain after zero outside-target margins')
        alternatives=[]
        for h,w,t in chosen:
            check(merged[owner[w]]==s,'literal source lookup')
            O=set(orbit(w));check(len(O)==4 and not O&reservoir,'fresh resources disjoint from full old and fresh selected alternatives');reservoir.update(O)
            K=partner(M,h);C=list(map(list,zip(*M)));rev=h[::-1]
            check(partner(K,h)==M,'actual complement involution')
            check(partner(C,rev)==list(map(list,zip(*K))),'actual converse square')
            check(all(K[i][j] for i,j in zip(h,h[1:])),'target retains the Hamiltonian path')
            check(word(C,rev)==mirror(w),'reflected clock gauge')
            alternatives.append({'path':h,'word':w,'target':t,'closed_ports':sorted(O),'target_matrix':K})
        anchor=F[0];anchor_owner=p[anchor]
        frame={i:('new_shared_anchor' if i==anchor_owner else f'private:{s}:{i}') for i in range(8)}
        check(len(set(frame.values()))==8 and p[anchor]==q[anchor],'injective frame, genuinely unchanged shared anchor')
        check(all((X[s]-X[t])%8==(X[s]-X[T[0]])%8 for t in T),'original ordinary selected X phase')
        domain.append({'source':s,'matrix':M,'label':label(X,Y,s),'X':X[s],
                       'fixed_clock_owners':[(i,p[i]) for i in F],'shared_anchor_clock':anchor,'owner_frame':frame,
                       'cost':r['cost'],'raw_domain_size':len(raw),'raw_targets':dict(Counter(t for h,w,t in raw)),
                       'complete_target_path_capacities':{str(t):caps[t] for t in T},'domain':alternatives})
        raw_counts.append(len(raw))
    # Check true original black capacities and authentic full-bank generators.
    edges=list(csv.DictReader(gzip.decompress(read('metagraph_order8_20261023/merged_edges.tsv.gz')).decode().splitlines(),delimiter='\t'))
    E=[tuple(int(e[k]) for k in ('source','target','blue','black')) for e in edges]
    edge={(i,j):(b,k) for i,j,b,k in E}
    for d in domain:
        s=d['source']
        for t in T:check(edge[tuple(sorted((s,t)))]==(0,2*d['complete_target_path_capacities'][str(t)]),'original genuine paired-black edge capacity')
    chars=[X]+[[64*y for y in row] for row in Y[:4]]
    for chi in chars:
        L=[0]*len(nodes)
        for i,j,b,k in E:z=(b+k)*(chi[i]-chi[j]);L[i]+=z;L[j]-=z
        check(all(z%256==0 for z in L),'actual original harmonic generator')
        check(2*(chi[T[0]]-chi[T[1]])%256==0,'full paired target charge equality')
    def current(bits,chi):
        row=defaultdict(int)
        for j,s in enumerate(S):
            t=T[1] if bits>>j&1 else T[0];v=2*(chi[s]-chi[t]);row[s]+=v;row[t]-=v
        return dict(row)
    start=12;finish=3;baseline=[current(start,chi) for chi in chars];accepted=[];separators=[]
    for bits in range(16):
        if bits.bit_count()!=2:continue
        rows=[current(bits,chi) for chi in chars];diff=[{u:row.get(u,0)-base.get(u,0) for u in set(row)|set(base)} for row,base in zip(rows,baseline)]
        bad=[(j,u,z%256) for j,row in enumerate(diff) for u,z in row.items() if z%256]
        if bad:separators.append({'mask':bits,'character':bad[0][0],'physical_owner':bad[0][1],'residual_mod256':bad[0][2]})
        else:accepted.append(bits)
    check(accepted==[finish,start],'complete fresh four-source margin2 fibre is a primitive two-state return')
    debt={u:current(finish,X).get(u,0)-current(start,X).get(u,0) for u in S+list(T)}
    check(debt[T[0]]==256 and debt[T[1]]==-256 and all(debt[s]==0 for s in S),'exact physical integer X debt')
    check(sum(d['cost'] for d in domain)==23,'actual owner-clock work')
    check(len(reservoir)==256 and len(old_reservoir)==224,'64 reserved resources for old28 plus fresh4 selected alternatives')
    return {'status':'FINITE-EXACT actual order8 primitive four-source570/V return; old28 pins unchanged',
            'positive_570_to_V':sorted(A),'negative_V_to_570':sorted(B),'target_pair':T,'start_mask':start,'finish_mask':finish,
            'accepted_masks':accepted,'proper_balanced_separators':separators,'owner_clock_work':23,'raw_domain_sizes':raw_counts,
            'literal_sources':domain,'integer_X_current_delta':debt,'fresh_reserved_resources':8,'fresh_closed_ports':32,
            'combined_reserved_resources':64,'combined_closed_ports':256,'combined_arrival_resources':32,'combined_arrival_closed_ports':128,
            'scope':'old28 frozen at any chosen native state; fresh sources retain their own common clocks and one new shared owner, not fixed first/last roles; disjoint compatible remainder required',
            'input_revision':REV,'inputs':INPUTS,'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
