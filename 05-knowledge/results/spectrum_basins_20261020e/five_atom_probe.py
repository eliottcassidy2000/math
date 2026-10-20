"""Discovery: even prefix-preserving black reductions of strong order-five atoms."""
from inventory_probe import *
from scc_probe import canon

def strong(t):
    reach=[row|(1<<i) for i,row in enumerate(t)]
    for v in range(len(t)):
        for u in range(len(t)):
            if reach[u]>>v&1:reach[u]|=reach[v]
    return all(row==(1<<len(t))-1 for row in reach)

def atoms():
    out={}
    pairs=[(i,j) for i in range(5) for j in range(i+2,5)]
    for mask in range(1<<len(pairs)):
        t=list(word_tournament((1,)*5))
        for bit,(i,j) in enumerate(pairs):
            if mask>>bit&1:t[i]^=1<<j;t[j]^=1<<i
        t=tuple(t)
        if strong(t):out.setdefault(canon(t),t)
    return list(out.values())

def main():
    rows=[]; target=word_tournament((1,1,1,1,3)); edge=frozenset((0,1))
    targets=[]
    for p in permutations(range(2,7)):
        t=relabel(target,(0,1)+p)
        if t not in targets:targets.append(t)
    back={}
    for t in targets:
        for state,route in ball2(t,edge).items():back.setdefault(state,(t,route))
    for b in atoms():
        source=join((2,0),b)
        forward=ball2(source,edge)
        common=forward.keys()&back.keys()
        if common:
            mid=min(common); t,rev=back[mid]
            route=forward[mid]+tuple(tuple(reversed(q)) for q in reversed(rev))
            state=source
            for p in route:
                need(not blue(state,p),'blue discovery route')
                state=reverse_path(state,p)
            need(state==t,'target mismatch')
            row=dict(atom=b,source=source,target=t,paths=route)
        else:row=dict(atom=b,found=False,forward_states=len(forward))
        rows.append(row);print(json.dumps(row),flush=True)
    HERE.joinpath('five_atom_probe.json').write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':main()
