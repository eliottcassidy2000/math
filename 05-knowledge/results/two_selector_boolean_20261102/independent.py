"""Independent complete subset-DP and bounded-word observer certificate."""
from collections import Counter
from pathlib import Path
import json
from experiment import rows, ham, need, R

HERE=Path(__file__).resolve().parent


def main():
    n=7;N=1<<15;full=N-1;graphs=[rows(n,x,True) for x in range(N)]
    counts=[ham(g) for g in graphs]
    claimed=json.loads((HERE/'experiment.json').read_text())['carrier']
    need({str(k):v for k,v in Counter(counts).items()}==claimed['H_histogram'],
         'all 32768 independent subset counts histogram')
    from experiment import path_cube_counts
    need(counts==path_cube_counts(n,True),'every path-subcube count versus endpoint DP')
    free=[(a,b) for a in range(n) for b in range(a+2,n)]
    def transition(order):
        result=[]
        for g in graphs:
            if any(not(g[order[i]]>>order[i+1]&1) for i in range(n-1)):
                result.append(-1)
            else:
                target=0
                for k,(a,b) in enumerate(free):
                    if g[order[a]]>>order[b]&1:target|=1<<k
                result.append(target)
        return result
    C=[x^full for x in range(N)];r=transition(R);s=transition(claimed['second_selector'])
    def word_partition(transitions,depth):
        ids={};labels=[]
        for x in range(N):
            frontier=[x];signature=[counts[x]]
            for _ in range(depth):
                frontier=[t[y] if y>=0 else -1 for y in frontier for t in transitions]
                signature.extend(counts[y] if y>=0 else -1 for y in frontier)
            key=tuple(signature)
            if key not in ids:ids[key]=len(ids)
            labels.append(ids[key])
        # Every split has a bounded distinguishing word. Congruence gives
        # sufficiency for all future words, independently of refinement rounds.
        destinations={}
        for x,cls in enumerate(labels):
            value=tuple(labels[t[x]] if t[x]>=0 else -1 for t in transitions)
            if cls in destinations:need(destinations[cls]==value,'bounded-word partition is closed')
            else:destinations[cls]=value
        return len(ids)
    one=word_partition([C,r],3);two=word_partition([C,r,s],4)
    need(one==10455 and two==15874,'independent bounded-word minimal observers')
    result={'status':'PASS','all_count_states':N,'one_selector_word_bound':3,
            'two_selector_word_bound':4,'one_selector_classes':one,'two_selector_classes':two,
            'method':'endpoint subset DP; all command words through the bounds; exact congruence check'}
    (HERE/'independent.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8')
    print(json.dumps(result))


if __name__=='__main__':main()
