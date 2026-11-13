"""Independent neutral-component and selector controls; no selector imports."""
from itertools import combinations
from math import gcd
from pathlib import Path
import json
from root_phase import state, require


def H(rows):
    n=len(rows); full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,full+1):
        for v,c in enumerate(dp[mask]):
            if not c:continue
            todo=rows[v] & (full^mask)
            while todo:
                bit=todo&-todo;todo-=bit
                dp[mask|bit][bit.bit_length()-1]+=c
    return sum(dp[-1])


def flip(rows,a,b):
    out=list(rows);out[a]^=1<<b;out[b]^=1<<a
    return tuple(out)


def twins(rows,a,b):
    return not ((rows[a]^rows[b]) & ~((1<<a)|(1<<b)))


def sccs(rows):
    n=len(rows); reach=[row|(1<<v) for v,row in enumerate(rows)]
    for k in range(n):
        for v in range(n):
            if reach[v]>>k&1:reach[v]|=reach[k]
    remaining=set(range(n)); parts=[]
    while remaining:
        v=min(remaining)
        part={u for u in remaining if reach[v]>>u&1 and reach[u]>>v&1}
        remaining-=part;parts.append(tuple(sorted(part)))
    return sorted(parts,key=lambda part:-sum(bool(rows[part[0]]>>other[0]&1)
                                             for other in parts if other!=part))


def skeleton(rows):
    answer=[];run=[]
    for part in sccs(rows):
        if len(part)==1:run.extend(part)
        else:
            if run:answer.append(('run',tuple(sorted(run))));run=[]
            answer.append(('strong',part))
    if run:answer.append(('run',tuple(sorted(run))))
    return answer


def relabel(rows,order):
    return tuple(sum(((rows[old]>>other)&1)<<j for j,other in enumerate(order))
                 for old in order)


def marked(rows,root):
    return state(relabel(rows,[root]+[v for v in range(len(rows)) if v!=root]))


def main():
    counts={'labelled_tournaments':0,'cross_scc_flips':0,'strict_cross_flips':0,
            'adjacent_singleton_neutral':0,'all_neutral_edges':0}
    for n in range(1,6):
        pairs=list(combinations(range(n),2)); all_rows=[]
        for code in range(1<<len(pairs)):
            rows=[0]*n
            for bit,(a,b) in enumerate(pairs):
                if code>>bit&1:rows[a]|=1<<b
                else:rows[b]|=1<<a
            all_rows.append(tuple(rows))
        values={rows:H(rows) for rows in all_rows}
        for rows in all_rows:
            parts=sccs(rows); partof={v:i for i,p in enumerate(parts) for v in p}
            for a,b in pairs:
                other=flip(rows,a,b)
                if partof[a]!=partof[b]:
                    counts['cross_scc_flips']+=1
                    i,j=sorted((partof[a],partof[b]))
                    exception=j==i+1 and len(parts[i])==len(parts[j])==1
                    require((values[other]==values[rows])==exception,'cross-SCC equality')
                    if exception:counts['adjacent_singleton_neutral']+=1
                    else:
                        require(values[other]>values[rows],'cross-SCC strict increase')
                        counts['strict_cross_flips']+=1
                if values[other]==values[rows]:
                    counts['all_neutral_edges']+=1
                    require(skeleton(rows)==skeleton(other),'neutral SCC skeleton')
            counts['labelled_tournaments']+=1
    receipt=json.loads((Path(__file__).parent.parent/'selector/receipt.json').read_text())
    sources=[]
    for index,src in enumerate(receipt['sources']):
        rows=tuple(src['rows']);root=src['root'];h=H(rows)
        require(h==881,'source H')
        neutral=[(a,b) for a,b in combinations(range(len(rows)),2) if H(flip(rows,a,b))==h]
        require(all(twins(rows,a,b) for a,b in neutral),'only twin source neighbors')
        expected=[(4,5),(7,8)] if index==0 else [(5,6)]
        require(neutral==expected,'expected neutral edges')
        seen={rows};todo=[rows]
        while todo:
            here=todo.pop()
            for a,b in combinations(range(len(rows)),2):
                nxt=flip(here,a,b)
                if nxt not in seen and H(nxt)==h:
                    seen.add(nxt);todo.append(nxt)
        require(len(seen)==(4 if index==0 else 2),'full neutral component')
        original=marked(rows,root)
        require(all(marked(t,root)==original for t in seen),'root profile constant on component')
        for v,p in enumerate(src['root_profiles']):
            require(marked(rows,v)==tuple(p),'independent full root atlas')
        order=list(range(len(rows)));order[root],order[3]=order[3],order[root]
        require(relabel(rows,order)==tuple(src['prepared_rows']),'literal role swap')
        prepared=marked(src['prepared_rows'],root)
        require(prepared==tuple(src['prepared_profile']),'prepared state')
        for ret in src['literal_returns']:
            require(marked(ret['rows'],ret['root'])==tuple(ret['profile']),'actual return state')
        h,s,q,F=prepared;d=F-q;D=F+s-h-q
        for k in range(1,201):
            t=d*2**k+q-D;b=3*t+2*D
            good=gcd(b,t)==1 and t%3==0 and b%3!=0
            require(good==(k%2==1 if index==0 else k%10!=7),'explicit complete sieve controls')
        sources.append({'root':root,'neutral_arcs':neutral,'component_size':len(seen),
                        'root_state_constant':original,'prepared_state':prepared,
                        'literal_returns_checked':len(src['literal_returns'])})
    result={'status':'FINITE-EXACT','optimized':not __debug__,'counts':counts,
            'sources':sources,'scope':'All labelled tournaments through n=5; full labelled H881 neutral components at fixed orders; literal return graphs k=1,2,3.'}
    name='neutral_scc.optimized.json' if not __debug__ else 'neutral_scc.json'
    Path(__file__).with_name(name).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
