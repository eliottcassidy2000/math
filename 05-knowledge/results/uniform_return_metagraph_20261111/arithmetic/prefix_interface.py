"""Complete seven-coordinate prefix response and physical S5a/S5b switch."""
from itertools import permutations
from pathlib import Path
import json
import audit_defect as audit
import dyadic_defect as dyadic

HERE=Path(__file__).resolve().parent
S5A=(2,4,9,19,7)
FLIPPED=(6,4,8,19,7)  # Reverse the inherited arc2->0, keeping every name.
S5B=(2,4,25,3,11)


def profile(rows,root=None):
    n=len(rows);root=n-1 if root is None else root;full=(1<<n)-1
    starts=audit.starts(rows)
    ends=audit.starts([sum(((rows[w]>>v)&1)<<w for w in range(n)) for v in range(n)])
    H=[1]+[sum(x) for x in starts[1:]];deletion=full^(1<<root)
    U=sum(starts[S][root]*H[full^S] for S in range(1,full+1) if S>>root&1)
    V=sum(ends[S][root]*H[full^S] for S in range(1,full+1) if S>>root&1)
    K=sum(H[S]*H[deletion^S] for S in dyadic.submasks(deletion))
    Q=sum(H[S]*H[full^S] for S in range(full+1))
    return {'H':H[full],'S':starts[full][root],'E':ends[full][root],
            'd':H[deletion],'U':U,'V':V,'K':K},Q


def predicted(P,V):
    p,Q=profile(P);kernel,_=profile(V)
    h=p['H'];q0=kernel['d'];QW=kernel['K']
    C=kernel['V']-kernel['H']-q0+kernel['S'];a0=kernel['V']-q0
    H=a0*Q+(q0-C)*h
    return {'H':H,'S':h*q0,'E':h*kernel['E'],'d':h*q0,
            'U':Q*QW+h*(kernel['U']-QW),'V':H+h*C,'K':Q*QW}


def amalgam(A,B):
    # Both input roots are last. Keep all A names and append B's nonroots.
    a=len(A);b=len(B);root=a-1
    rows=list(A)+[0]*(b-1)
    for v in range(a-1):rows[v]|=((1<<(b-1))-1)<<a
    for v in range(b-1):
        for w in range(b-1):
            if B[v]>>w&1:rows[a+v]|=1<<(a+w)
        if B[v]>>(b-1)&1:rows[a+v]|=1<<root
        else:rows[root]|=1<<(a+v)
    return rows,root


def main():
    old,qold=profile(S5A);new,qnew=profile(FLIPPED)
    audit.check((old['H'],qold,new['H'],qnew)==(9,84,9,84),'equal prefix H,Q2')
    changes=[(v,w) for v in range(5) for w in range(5) if (S5A[v]^FLIPPED[v])>>w&1]
    audit.check(changes==[(0,2),(2,0)],'one physical named arc')
    maps=[]
    for p in permutations(range(5)):
        mapped=tuple(sum(((FLIPPED[p[v]]>>p[w])&1)<<w for w in range(5)) for v in range(5))
        if mapped==S5B:maps.append(p)
    audit.check(bool(maps),'flipped core belongs to S5b ordinary class')
    controls=[]
    for left,right in ((0,0),(1,0),(0,1),(1,1)):
        for tail in ((),(1,)):
            V=audit.graph(tail);pair=[]
            for core in (S5A,FLIPPED):
                P=dyadic.ordinal([(0,)*left,core,(0,)*right])
                rows=dyadic.prefix_insert(P,V)
                actual,_=profile(rows)
                audit.check(actual==predicted(P,V),'literal seven-coordinate formula')
                pair.append(actual)
            audit.check(pair[0]==pair[1],'full interface-neutral physical flip')
            controls.append({'left':left,'right':right,'tail':tail,'common_seven':pair[0]})
    # Same total gap, six distinct ordered-core/gap sources, identical state.
    rearrangements=[]
    for cores in (((2,4,1),(2,4,9,3)),((2,4,9,3),(2,4,1))):
        for gap in range(3):
            parts=[]
            for j in range(3):
                if j==gap:parts.append((0,))
                if j<2:parts.append(cores[j])
            P=dyadic.ordinal(parts);V=audit.graph(())
            state,_=profile(dyadic.prefix_insert(P,V))
            audit.check(state==predicted(P,V),'permutation/gap response')
            rearrangements.append({'core_orders':list(map(len,cores)),'gap':gap,'seven':state})
    audit.check(all(r['seven']==rearrangements[0]['seven'] for r in rearrangements),'same interface across six sources')
    context=(2,4,1);grafts=[]
    for side in ('left','right'):
        states=[]
        for core in (S5A,FLIPPED):
            T=dyadic.prefix_insert(core,audit.graph(()))
            rows,root=amalgam(context,T) if side=='left' else amalgam(T,context)
            state,_=profile(rows,root);states.append(state)
        audit.check(states[0]==states[1],'literal left/right root amalgamation')
        grafts.append({'context_rows':context,'context_side':side,'common_seven':states[0]})
    report={'status':'PASS literal full-seven controls; all-order preservation proved in DYADIC_DEFECT section6',
            'coordinate_order':['H','S','E','d','U','V','K'],
            'inherited_named_core_flip':{'source':S5A,'target':FLIPPED,'old_arc':[2,0],
                                         'target_to_S5b_relabellings':maps,'common_H':9,'common_Q2':84},
            'physical_context_flip_pairs':controls,'same_interface_rearrangements':rearrangements,
            'literal_amalgam_context_pairs':grafts,
            'scope':'A one-arc replacement inside an ordinary prefix preserves all seven original-root coordinates. Prefix core permutation is observer-equivalent, but no neutral physical permutation path is claimed.'}
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'prefix_interface.json').write_text(text,encoding='utf-8')
    print(json.dumps({'status':report['status'],'physical_pairs':len(controls),'rearrangements':len(rearrangements),'amalgam_pairs':len(grafts)}))


if __name__=='__main__':main()
