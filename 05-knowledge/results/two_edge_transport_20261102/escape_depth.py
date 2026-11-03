"""Exact shortest distance and least maximum-error escape, retaining every step."""
from collections import Counter, deque
from heapq import heappop, heappush
from itertools import combinations
from pathlib import Path
import json
from scout import rows,h_count
HERE=Path(__file__).resolve().parent

def main():
    counts=[h_count(rows(6,x)) for x in range(1<<15)]
    traps=[r['source'] for r in json.loads((HERE/'escape_scout.json').read_text())['records']]
    goals=[x for x,h in enumerate(counts) if h==17]
    distance=[-1]*len(counts);nxt=[None]*len(counts);q=deque(goals)
    for x in goals:distance[x]=0
    while q:
        x=q.popleft()
        for e in range(15):
            y=x^(1<<e)
            if distance[y]<0:distance[y]=distance[x]+1;nxt[y]=x;q.append(y)
    # First minimize the maximum error. Then BFS inside that sublevel set.
    # A lexicographic (peak,length) Dijkstra is not generally isotone under max.
    peaks=[10**9]*len(counts);heap=[]
    for x in goals:peaks[x]=0;heappush(heap,(0,x))
    while heap:
        peak,x=heappop(heap)
        if peak!=peaks[x]:continue
        for e in range(15):
            y=x^(1<<e);new=max(peak,(counts[y]-17)**2)
            if new<peaks[y]:peaks[y]=new;heappush(heap,(new,y))
    constrained={}
    for peak in sorted(set(peaks[x] for x in traps)|{16}):
        dist=[-1]*len(counts);links=[None]*len(counts);q=deque(goals)
        for x in goals:dist[x]=0
        while q:
            x=q.popleft()
            for e in range(15):
                y=x^(1<<e)
                if dist[y]<0 and (counts[y]-17)**2<=peak:
                    dist[y]=dist[x]+1;links[y]=x;q.append(y)
        constrained[peak]=(dist,links)
    cost=[(peaks[x],constrained[peaks[x]][0][x]) if x in traps else None for x in range(len(counts))]
    def walk(x,links):
        ans=[x]
        while links[x] is not None:x=links[x];ans.append(x)
        return {'masks':ans,'counts':[counts[y] for y in ans],
                'flipped_edge_indices':[(a^b).bit_length()-1 for a,b in zip(ans,ans[1:])]}
    summary=Counter((distance[x],*cost[x]) for x in traps)
    reps=[]
    for key in sorted(summary):
        x=next(x for x in traps if (distance[x],*cost[x])==key)
        shell3=[x^sum(1<<e for e in es) for es in combinations(range(15),3)]
        reps.append({'key':key,'source':x,'count':summary[key],
                     'shortest':walk(x,nxt),'minimax':walk(x,constrained[peaks[x]][1]),
                     'three_flip_targets':sum(counts[y]==17 for y in shell3)})
    result={'all_source_distance_histogram':dict(sorted(Counter(distance).items())),
            'trap_types':reps,'goal_count':len(goals),
            'traps_reachable_at_peak16':sum(constrained[16][0][x]>=0 for x in traps)}
    (HERE/'escape_depth.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
