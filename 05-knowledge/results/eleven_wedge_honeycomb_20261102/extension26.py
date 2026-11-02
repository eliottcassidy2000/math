#!/usr/bin/env python3
"""A native extension hostile to the fixed-G25 complement carrier."""
import json
from itertools import combinations
from math import isqrt
def need(b,s):
    if not b:raise RuntimeError(s)
def bank(n):
    E=[(i,j) for i,j in combinations(range(1,n+1),2) if isqrt(i+j)**2==i+j]
    A={i:[] for i in range(1,n+1)}
    for i,j in E:A[i].append(j);A[j].append(i)
    need([v for v in A if len(A[v])==1]==[18],'forced endpoint')
    paths=[]
    def dfs(p,seen):
        if len(p)==n:paths.append(p);return
        for v in A[p[-1]]:
            if v not in seen:dfs(p+(v,),seen|{v})
    dfs((18,),{18})
    return E,sorted(paths+[p[::-1] for p in paths])
def encode(p,E):
    mask=val=0
    for i,j in zip(p,p[1:]):
        k=E.index(tuple(sorted((i,j))));mask|=1<<k
        if i<j:val|=1<<k
    return mask,val
def main():
    E,P=bank(25);EE,PP=bank(26);new=[e for e in EE if e not in E]
    need(len(P)==20 and len(PP)==24 and new==[(10,26),(23,26)],'source banks')
    old=[encode(p,E) for p in P];ext=[encode(p,E+new) for p in PP];full=(1<<len(E))-1
    def legal(x,patterns):return tuple(i for i,(m,v) in enumerate(patterns) if x&m==v)
    def state(x):
        L=legal(x,old)
        return L,tuple(legal(x^(full^old[i][0]),old) for i in L)
    active=set()
    for m,v in old:
        free=[i for i in range(len(E)) if not m>>i&1]
        for bits in range(1<<len(free)):
            active.add(v|sum(1<<k for j,k in enumerate(free) if bits>>j&1))
    states={x:state(x) for x in active};need(len(active)==4864 and len(set(states.values()))==32,'inherited active carrier')
    need(len({(states[x],tuple(states[x^(full^old[i][0])] for i in states[x][0])) for x in active})==32,'old carrier transition closure')
    groups={};scores={}
    for x in sorted(active):
        for a in range(4):
            h=len(legal(x|(a<<32),ext));scores[x,a]=h
            groups.setdefault((states[x],a),{}).setdefault(h,x)
    split=[(s,a,d) for (s,a),d in groups.items() if len(d)>1]
    # Minimal active-source observer for arbitrary old commands followed by
    # one terminal extension26 query, whose two new edge bits are supplied.
    ordered=sorted(active)
    outputs={x:(states[x],tuple(scores[x,a] for a in range(4))) for x in ordered}
    def compress(signatures):
        names={s:i for i,s in enumerate(sorted(set(signatures.values())))}
        return {x:names[s] for x,s in signatures.items()}
    classes=compress(outputs);rounds=[len(set(classes.values()))]
    for _ in range(len(active)):
        sig={x:(classes[x],tuple(classes[x^(full^old[i][0])] for i in states[x][0])) for x in ordered}
        nxt=compress(sig);rounds.append(len(set(nxt.values())))
        if rounds[-1]==rounds[-2]:break
        classes=nxt
    else:raise RuntimeError('refinement did not stabilize')
    # Start with just current legal names and the terminal outputs instead;
    # this must recover the same equivalence without assuming the old quotient.
    alt=compress({x:(states[x][0],tuple(scores[x,a] for a in range(4))) for x in ordered})
    for _ in range(len(active)):
        nxt=compress({x:(alt[x],tuple(alt[x^(full^old[i][0])] for i in states[x][0])) for x in ordered})
        if len(set(nxt.values()))==len(set(alt.values())):break
        alt=nxt
    need(len(set(zip(classes.values(),alt.values())))==len(set(classes.values()))==len(set(alt.values())),'two observer initializations')
    need(rounds==[36,38,38] and len(split)==4 and len({s for s,a,d in split})==4,'declared native refinement census')
    # Reproduce the incoming full G26 carrier, then preserve that entire
    # target behaviour across the native birth of26, not merely its count.
    active26=set();full26=(1<<34)-1
    for m,v in ext:
        free=[i for i in range(34) if not m>>i&1]
        for bits in range(1<<len(free)):
            active26.add(v|sum(1<<k for j,k in enumerate(free) if bits>>j&1))
    L26={x:legal(x,ext) for x in sorted(active26)}
    c26=compress(L26);r26=[len(set(c26.values()))]
    for _ in range(len(active26)):
        nc=compress({x:(c26[x],tuple(c26[x^(full26^ext[i][0])] for i in L26[x])) for x in L26})
        r26.append(len(set(nc.values())))
        if r26[-1]==r26[-2]:break
        c26=nc
    need(len(active26)==11232 and r26==[36,66,70,70],'incoming full target carrier')
    targets={x:tuple(c26.get(x|(a<<32),-1) for a in range(4)) for x in ordered}
    birth=compress({x:(states[x],targets[x]) for x in ordered});rb=[len(set(birth.values()))]
    for _ in range(len(active)):
        nb=compress({x:(birth[x],tuple(birth[x^(full^old[i][0])] for i in states[x][0])) for x in ordered})
        rb.append(len(set(nb.values())))
        if rb[-1]==rb[-2]:break
        birth=nb
    need(len(set(zip(birth.values(),classes.values())))==len(set(birth.values())),'birth carrier refines count-only carrier')
    need(rb==[38,38] and len(set(zip(birth.values(),classes.values())))==38,'full birth and terminal-count observers coincide')
    birth_table=[];birth_mass={}
    for label in sorted(set(birth.values())):
        members=[x for x in ordered if birth[x]==label];x=members[0]
        transitions={str(i):birth[x^(full^old[i][0])] for i in states[x][0]}
        need(all(targets[y]==targets[x] and {str(i):birth[y^(full^old[i][0])] for i in states[y][0]}==transitions for y in members),'actual prebirth state table')
        birth_table.append({'class':label,'source_weight':len(members),'representative_bits':x,
                            'old_legal_names':states[x][0],'old_command_targets':transitions,
                            'target26_classes_by_new_bits':targets[x]})
        for target in targets[x]:birth_mass[target]=birth_mass.get(target,0)+len(members)
    need(sum(r['source_weight'] for r in birth_table)==4864 and sum(birth_mass.values())==19456,'source-weighted birth mass')
    weights={r['class']:r['source_weight'] for r in birth_table}
    need(all(weights[t]==r['source_weight'] for r in birth_table for t in r['old_command_targets'].values()),'partial involution preserves source weights')
    # Complete ALL2^32 old sources without enumerating the huge dead fibre.
    # Only old-active sources or projections of new-active sources can ever
    # answer a guard positively. Everything else is one weighted dead state.
    relevant=active|{z&full for z in active26};domain=sorted(relevant)+[-1]
    def old_names(x):return states[x][0] if x in states else ()
    full_targets={x:tuple(c26.get(x|(a<<32),-1) for a in range(4)) if x!=-1 else (-1,)*4 for x in domain}
    joint=compress({x:(old_names(x),full_targets[x]) for x in domain});rjoint=[len(set(joint.values()))]
    for _ in range(len(domain)):
        nj=compress({x:(joint[x],tuple(joint[x^(full^old[i][0])] for i in old_names(x))) for x in domain})
        rjoint.append(len(set(nj.values())))
        if rjoint[-1]==rjoint[-2]:break
        joint=nj
    need(rjoint==[89,97,97] and len(relevant)==14048,'complete native birth census')
    need({x^full for x in relevant}==relevant,'relevant-source converse closure')
    converse26={-1:-1}
    for x in L26:
        label,target=c26[x],c26[x^full26]
        need(label not in converse26 or converse26[label]==target,'target converse descends')
        converse26[label]=target
    for x in domain:
        y=x^full if x!=-1 else -1
        need(all(full_targets[y][a^3]==converse26[full_targets[x][a]] for a in range(4)),'native birth converse square')
    symbolic_weight=2**32-len(relevant);need(symbolic_weight>0,'residual dead source mass')
    joint_table=[];joint_mass={}
    for label in sorted(set(joint.values())):
        members=[x for x in domain if joint[x]==label];x=members[0]
        weight=sum(symbolic_weight if y==-1 else 1 for y in members)
        transitions={str(i):joint[x^(full^old[i][0])] for i in old_names(x)}
        reverse=joint[x^full if x!=-1 else -1]
        need(all(full_targets[y]==full_targets[x] and {str(i):joint[y^(full^old[i][0])] for i in old_names(y)}==transitions for y in members),'all-source transition closure')
        need(all(joint[y^full if y!=-1 else -1]==reverse for y in members),'prebirth converse descends')
        joint_table.append({'class':label,'source_weight':weight,'representative_bits':x,
                            'old_legal_names':old_names(x),'old_command_targets':transitions,
                            'converse_class':reverse,'target26_classes_by_new_bits':full_targets[x]})
        for target in full_targets[x]:joint_mass[target]=joint_mass.get(target,0)+weight
    expected_mass={-1:2**34-len(active26)}
    for label in c26.values():expected_mass[label]=expected_mass.get(label,0)+1
    need(sum(r['source_weight'] for r in joint_table)==2**32 and joint_mass==expected_mass,'exact full source/target mass conservation')
    need(sum(bool(r['old_legal_names']) for r in joint_table)==38 and len(relevant-active)==9184,'38 active and59 old-dead classes')
    full_weights={r['class']:r['source_weight'] for r in joint_table}
    need(all(full_weights[t]==r['source_weight'] for r in joint_table for t in r['old_command_targets'].values()),'full weighted involutions')
    need(all(full_weights[r['converse_class']]==r['source_weight'] for r in joint_table),'converse source weights')
    old_active_new_active=sum((z&full) in active for z in active26)
    need(old_active_new_active==2048 and sum(m for c,m in birth_mass.items() if c!=-1)==2048,'active/active incidence')
    collision=None;seen={}
    for x in ordered:
        key=classes[x]
        if key in seen and birth[seen[key]]!=birth[x]:collision=(seen[key],x);break
        seen[key]=x
    witness=None
    for x in sorted(active):
        for k in range(32):
            y=x^(1<<k)
            if y<=x or y not in states or states[y]!=states[x]:continue
            for a in range(4):
                if scores[x,a]!=scores[y,a]:witness=(x,y,k,a);break
            if witness:break
        if witness:break
    need(witness is not None,'one-edge native hostile')
    x,y,k,a=witness
    # Independent literal orientation/path comparison for the two witnesses.
    def arcs(z):return {(i,j) if z>>k&1 else (j,i) for k,(i,j) in enumerate(E+new)}
    sources=[]
    for z in [x,y]:
        w=z|(a<<32);A=arcs(w)
        good=[p for p in PP if all(e in A for e in zip(p,p[1:]))]
        need(len(good)==scores[z,a],'literal consequence cross-check')
        sources.append({'guard_bits':z,'H26':len(good),'legal_paths26':good})
    return {'status':'passed: complete weighted native birth interface','scope':'All2^32 G25 orientations via14048 explicit relevant sources and one weighted residual state; all G26 target behaviours via11232 active sources and one dead state. Also retains the separate active25-only38-state result.',
            'edge_order25':E,'new_edges':new,'path_counts_directed':[len(P),len(PP)],
            'paths26_with_new_endpoint':sum(p[0]==26 or p[-1]==26 for p in PP),
            'paths26_with_new_interior':sum(p[0]!=26 and p[-1]!=26 for p in PP),
            'active_carrier_states':32,'fixed_attachment_groups':len(groups),'split_groups':len(split),
            'split_old_states':len({s for s,a,d in split}),
            'terminal_extension_refinement_counts':rounds,'minimal_active_terminal_extension_states':rounds[-1],
            'full_target26_refinement_counts':r26,'active_target26_sources':len(active26),
            'full_target26_states_with_dead':r26[-1]+1,
            'birth_refinement_counts':rb,'minimal_active_prebirth_states':rb[-1],
            'prebirth_controller_table':birth_table,'target26_mass_from_this_birth':birth_mass,
            'target26_mass_scope':'Pushforward of active25 source weights and all four attachments; do not replace with the full26 uniform fibre weights.',
            'full_prebirth_refinement_counts':rjoint,'full_prebirth_controller_table':joint_table,
            'explicit_relevant_old_sources':len(relevant),'symbolic_residual_source_weight':symbolic_weight,
            'new_active_from_old_active':old_active_new_active,'new_active_from_old_dead':len(active26)-old_active_new_active,
            'count_only_collision_for_full_birth':None if collision is None else [{'guard_bits':x,'four_target26_classes':targets[x],'count_only_class':classes[x],'birth_class':birth[x]} for x in collision],
            'witness':{'flipped_edge':E[k],'new_bits':a,'new_arcs':[e if a>>i&1 else e[::-1] for i,e in enumerate(new)],
                       'same_old_state':states[x],'old_legal_paths':[P[i] for i in states[x][0]],'sources':sources}}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
