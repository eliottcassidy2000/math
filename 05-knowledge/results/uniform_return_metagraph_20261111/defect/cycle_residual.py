"""Complete d=2,b=6 sector with residual TT_left -> C3 -> TT_right.

The two-inport positive automaton is extended by one exact C3 block.
Core contact words are quotiented only by its actual cyclic automorphisms.
Cost pruning uses the independently proved tagged-capacity monotonicity,
and the factor3 from disjoint C3 odd-cycle packings, never a fitted rule.
"""
from collections import Counter, deque
from itertools import permutations, product
from pathlib import Path
import json
import two_inport as two


CORE=(2,4,1)
CONTACTS=tuple(sorted({min(t,t[1:]+t[:1],t[2:]+t[:2])
                       for t in product(range(4),repeat=3)}))


def core_assignments():
    result=[]
    for mask in range(8):
        first=[i for i in range(3) if mask>>i&1]
        second=[i for i in range(3) if not(mask>>i&1)]
        for a in permutations(first):
            if any(not(CORE[u]>>v&1) for u,v in zip(a,a[1:])):
                continue
            for b in permutations(second):
                if all(CORE[u]>>v&1 for u,v in zip(b,b[1:])):
                    result.append((a,b))
    two.require(len(result)==12,'C3 ordered two-bin count')
    return result


ASSIGNMENTS=core_assignments()


def block_table(contact):
    table=[]
    for runs,_ in two.SKELETONS:
        origin_rows=[]
        for st in two.STATES:
            dest=Counter()
            for paths in ASSIGNMENTS:
                nxt=list(st)
                legal=True
                for j,(path,(left,right)) in enumerate(zip(paths,runs)):
                    if not path:
                        continue
                    if st[j]==0 and left>=0 and contact[path[0]]>>left&1:
                        legal=False
                        break
                    nxt[j]=1 if right<0 or contact[path[-1]]>>right&1 else 2
                if legal:
                    dest[two.STATE_ID[tuple(nxt)]]+=1
            origin_rows.append(tuple(sorted(dest.items())))
        table.append(tuple(origin_rows))
    return tuple(table)


TABLES={c:block_table(c) for c in CONTACTS}


def apply_block(weights,contact,cap=None):
    answer=[]
    for sk,old in enumerate(weights):
        new=[0]*9
        for i,count in enumerate(old):
            for j,multiplicity in TABLES[contact][sk][i]:
                new[j]+=count*multiplicity
        answer.append(tuple(min(cap,x) for x in new) if cap else tuple(new))
    return tuple(answer)


def literal(left,contact,right):
    word=left+contact+right
    rows=list(two.rows_for(word))
    start=3+len(left)
    for i in range(3):
        rows[start+i]&=~(7<<start)
        rows[start+i]|=CORE[i]<<start
    return tuple(rows)


def literal_defect(rows):
    hp,ends=two.path_tables(rows)
    full=(1<<len(rows))-1
    f=sum(ends[m][0]*hp[full^m] for m in range(1,full+1) if m&1)
    return f-hp[full],hp[full]


def controls():
    queries=0
    for ell in range(3):
        for rr in range(3-ell):
            for left in product(range(4),repeat=ell):
                for right in product(range(4),repeat=rr):
                    for contact in CONTACTS:
                        w=two.START
                        for x in left:
                            w=two.advance(w,x)
                        w=apply_block(w,contact)
                        for x in right:
                            w=two.advance(w,x)
                        target,_=literal_defect(literal(left,contact,right))
                        two.require(target==two.value(w),'C3 block/literal mismatch')
                        queries+=1
    return queries


def classify():
    # state: phase0 empty, phase1 nonempty left, phase2 core already inserted
    start=(0,two.START,False)
    states=[start]
    index={start:0}
    edges=[]
    for phase,weights,last in states:
        row=[]
        if phase<2:
            options=[c for c in CONTACTS if phase or c==(3,3,3)]
            for c in options:
                w=apply_block(weights,c,7)
                if two.value(w)>6:
                    continue
                nxt=(2,w,any(c))
                if nxt not in index:
                    index[nxt]=len(states);states.append(nxt)
                row.append(('C'+','.join(map(str,c)),3,index[nxt]))
        for letter in range(4):
            if phase==0 and letter!=3:
                continue
            w=two.advance(weights,letter,7)
            bound=2 if phase<2 else 6
            if two.value(w)>bound:
                continue
            nxt=(max(1,phase),w,bool(letter) if phase==2 else False)
            if nxt not in index:
                index[nxt]=len(states);states.append(nxt)
            row.append((str(letter),1,index[nxt]))
        edges.append(row)
    accepted={i for i,(phase,w,last) in enumerate(states) if phase==2 and last and two.value(w)==6}
    reverse=[[] for _ in states]
    for i,row in enumerate(edges):
        for _,_,j in row:
            reverse[j].append(i)
    useful=set(accepted);queue=deque(accepted)
    while queue:
        for i in reverse[queue.popleft()]:
            if i not in useful:
                useful.add(i);queue.append(i)
    loops=[(i,label,degree) for i in sorted(useful) for label,degree,j in edges[i] if i==j]
    indegree={i:0 for i in useful}
    for i in useful:
        for _,_,j in edges[i]:
            if j in useful and i!=j:
                indegree[j]+=1
    queue=deque(i for i,d in indegree.items() if not d);topo=[]
    while queue:
        i=queue.popleft();topo.append(i)
        for _,_,j in edges[i]:
            if j in useful and i!=j:
                indegree[j]-=1
                if not indegree[j]:queue.append(j)
    patterns=[]
    if len(topo)==len(useful):
        def walk(i,parts):
            own=[label for node,label,_ in loops if node==i]
            two.require(len(own)<=1,'multiple self-loop alternatives')
            if own:parts=parts+[own[0]+'*']
            if i in accepted:patterns.append(parts)
            for label,_,j in edges[i]:
                if j in useful and i!=j:walk(j,parts+[label])
        if 0 in useful:walk(0,[])
    # Weight letters byone, the unique C3 marker bythree, root+ports bythree.
    by_degree=[Counter() for _ in range(41)]
    by_degree[3][0]=1
    counts=[0]*41
    for degree,dist in enumerate(by_degree):
        counts[degree]=sum(c for i,c in dist.items() if i in accepted)
        for i,c in dist.items():
            for _,cost,j in edges[i]:
                if j in useful and degree+cost<len(by_degree):
                    by_degree[degree+cost][j]+=c
    return {'status':'COMPLETE under proved d2,b6,C3-residual filters',
            'contact_orbits':len(CONTACTS),'reachable_states':len(states),'productive_states':len(useful),
            'self_loops':loops,'acyclic_after_self_loops':len(topo)==len(useful),
            'exact_disjoint_path_regex':patterns,'counts_by_total_order':counts,
            'productive_edges':[(i,label,cost,j) for i in sorted(useful) for label,cost,j in edges[i] if j in useful],
            'accepting':sorted(accepted),'initial':0}


def main():
    out=classify()
    out['literal_DP_controls']=controls()
    Path(__file__).with_name('cycle_residual.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='productive_edges'},indent=2))


if __name__=='__main__':main()
