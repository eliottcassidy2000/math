"""Independent deletion-lift and graph census; uses independently rebuilt DFA."""
from independent_audit import automaton,graph_counts,slot_counts,BASE
import itertools,json,collections,time,hashlib
start=time.time()
_,_,_,initial,T,O=automaton()
# Bitset, rather than frozenset, tracks possible once-deleted residual states.
queue=[(initial,0,0)];index={queue[0]:0};words=[()];bad=[]
for i,(s,deleted,length) in enumerate(queue):
    outputs={O[d] for d in range(len(T)) if (deleted>>d)&1}
    if O[s] in (1,2) and O[s] not in outputs:
        bad.append((O[s],words[i]));assert length <= (2 if O[s]==1 else 3)
    for a in range(4):
        nd=1<<s
        for d in range(len(T)):
            if (deleted>>d)&1:nd|=1<<T[d][a]
        new=(T[s][a],nd,min(4,length+1))
        if new not in index:index[new]=len(queue);queue.append(new);words.append(words[i]+(a,))
assert len(queue)==1066
expected=[(1,(2,)),(1,(3,0)),(2,(1,2,0)),(2,(2,1,2)),(2,(2,3,0)),(2,(3,0,2)),(2,(3,1,0)),(2,(3,2,0)),(2,(3,2,1))]
assert bad==expected
# Direct graph validation of irreducibility includes each possible base deletion.
for c,w in expected:
    assert graph_counts(w)[1]==c
    assert all(graph_counts(w[:i]+w[i+1:])[1]!=c for i in range(len(w)))
minima={};inventory={};checks=0;values=set()
for m in range(1,8):
    inv=collections.Counter()
    for w in itertools.product(range(4),repeat=m):
        s=initial
        for a in w:s=T[s][a]
        c=O[s]
        if c not in (1,2):continue
        h,gc=graph_counts(w);assert gc==c
        assert slot_counts(w)==(h,c)
        checks+=1;inv[c]+=1
        key=str((m,c))
        if key not in minima or h<minima[key][0]:minima[key]=[h,list(w)]
        if m<=6:values.add(h)
    inventory[str(m)]={str(c):v for c,v in inv.items()}
assert checks==1817 and 67 not in values
saved=json.loads((BASE/'deletion-source/lowcycle_minima.json').read_text())
assert saved['inventory']==inventory and saved['minima']==minima
out={'deletion_lift_states':len(queue),'irreducible_words':bad,'graph_checks':checks,'inventory':inventory,'minima':minima,'67_absent_through_base_6':True,'elapsed_seconds':time.time()-start}
(BASE/'deletion_results.json').write_text(json.dumps(out,indent=2)+'\n')
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((BASE/'deletion-source').iterdir()) if p.is_file()}
(BASE/'deletion_hashes.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('PASS deletion graph states',len(queue),'irreducible',bad,'independent graph checks',checks,'elapsed',time.time()-start)
