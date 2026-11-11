"""Finite-feedback positive defect compiler and capped three-inport strata.

Raw path-run automata are exactly lumped by weighted bisimulation before
the capped-vector reachability search. Numerical agreement is not used to
merge states. Only the transitive residual chart is enumerated here.
"""
from collections import Counter, deque
from itertools import permutations, product
from pathlib import Path
import argparse
import json


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def compile_machine(port_rows):
    d = len(port_rows)
    skeletons = []
    for p in permutations(range(d)):
        # A path starting at an inport: no backbone before first port.
        skeletons.append(tuple(zip(p,p[1:]))+((p[-1],-1),))
        for split in range(1,d):
            left,right = p[:split],p[split:]
            # First bin ends at its last port, second begins at its first.
            skeletons.append(((-1,left[0]),)+tuple(zip(left,left[1:]))
                             +tuple(zip(right,right[1:]))+((right[-1],-1),))
    run_states = tuple(product(range(3), repeat=d))
    ids = {s:i for i,s in enumerate(run_states)}
    width = len(run_states)
    raw_trans, raw_final, initial = [], [], Counter()
    for sk,runs in enumerate(skeletons):
        initial[sk*width] += 1
        for st in run_states:
            raw_final.append(int(all(v!=2 and (v!=0 or l<0 or r<0 or port_rows[l]>>r&1)
                                     for v,(l,r) in zip(st,runs))))
            transitions = []
            for letter in range(1<<d):
                dest = Counter()
                for j,(left,right) in enumerate(runs):
                    if st[j]==0 and left>=0 and letter>>left&1:
                        continue
                    nxt = list(st)
                    nxt[j] = 1 if right<0 or letter>>right&1 else 2
                    dest[sk*width+ids[tuple(nxt)]] += 1
                transitions.append(dest)
            raw_trans.append(transitions)
    # Exact right-language bisimulation, including all transition weights.
    classes = raw_final[:]
    while True:
        signatures = {}
        updated = []
        for i, row in enumerate(raw_trans):
            images = []
            for dest in row:
                grouped = Counter()
                for j,multiplicity in dest.items():
                    grouped[classes[j]] += multiplicity
                images.append(tuple(sorted(grouped.items())))
            sig = (raw_final[i],tuple(images))
            if sig not in signatures:
                signatures[sig] = len(signatures)
            updated.append(signatures[sig])
        if updated==classes:
            break
        classes = updated
    size = max(classes)+1
    trans, final = [None]*size, [None]*size
    for i,row in enumerate(raw_trans):
        c = classes[i]
        if trans[c] is not None:
            continue
        images=[]
        for dest in row:
            grouped=Counter()
            for j,multiplicity in dest.items():
                grouped[classes[j]] += multiplicity
            images.append(tuple(sorted(grouped.items())))
        trans[c],final[c] = tuple(images),raw_final[i]
    start = [0]*size
    for i,count in initial.items():
        start[classes[i]] += count
    return {'d':d,'raw_states':len(raw_trans),'states':size,'trans':trans,
            'final':tuple(final),'start':tuple(start),'port_rows':tuple(port_rows)}


def advance(machine, weights, letter, cap=None):
    answer=[0]*len(weights)
    for i,count in enumerate(weights):
        if not count:
            continue
        for j,multiplicity in machine['trans'][i][letter]:
            answer[j] += count*multiplicity
    if cap is not None:
        answer=[min(cap,x) for x in answer]
    return tuple(answer)


def evaluate(machine,word):
    weights=machine['start']
    for letter in word:
        weights=advance(machine,weights,letter)
    return sum(x*y for x,y in zip(weights,machine['final']))


def literal(word,port_rows):
    d=len(port_rows)
    n=d+len(word)+1
    rows=[sum(1<<j for j in range(d+1,n))]
    rows += [(row<<1)|1 for row in port_rows]
    rows += [0]*len(word)
    for j,letter in enumerate(word):
        v=d+1+j
        rows[v] |= sum(1<<w for w in range(v+1,n))
        for p in range(d):
            if letter>>p&1:
                rows[v] |= 1<<(p+1)
            else:
                rows[p+1] |= 1<<v
    return tuple(rows)


def direct(rows):
    n=len(rows)
    full=(1<<n)-1
    dp=[[0]*n for _ in range(full+1)]
    for v in range(n):
        dp[1<<v][v]=1
    for mask in range(1,full+1):
        for v,count in enumerate(dp[mask]):
            if not count:
                continue
            todo=rows[v]&~mask
            while todo:
                bit=todo&-todo
                todo-=bit
                dp[mask|bit][bit.bit_length()-1]+=count
    hps=[1]+[sum(v) for v in dp[1:]]
    f=sum(dp[mask][0]*hps[full^mask] for mask in range(1,full+1) if mask&1)
    return f-hps[full],hps[full]


def closure(port_rows,subset):
    while True:
        grown=subset
        for v,row in enumerate(port_rows):
            if subset>>v&1:
                grown|=row
        if grown==subset:
            return subset
        subset=grown


def classify(machine,budget,limit,monotone6=False):
    d=machine['d']
    full=(1<<d)-1
    closures=[closure(machine['port_rows'],s) for s in range(1<<d)]
    initial=(tuple(min(budget+1,x) for x in machine['start']),0,False)
    states=[initial]
    index={initial:0}
    edges=[]
    for weights,seen,last in states:
        row=[]
        for letter in range(1<<d):
            if monotone6 and (not letter or (not seen and letter!=full)
                              or (seen and letter & ~last)):
                row.append(-1)
                continue
            nxt=(advance(machine,weights,letter,budget+1),closures[seen|letter],
                 letter if monotone6 else bool(letter))
            if nxt not in index:
                if len(states)>=limit:
                    return {'status':'INCOMPLETE state limit','states':len(states),'budget':budget}
                index[nxt]=len(states)
                states.append(nxt)
            row.append(index[nxt])
        edges.append(row)
    accepted=[i for i,(w,seen,last) in enumerate(states)
              if seen==full and last and sum(x*y for x,y in zip(w,machine['final']))==budget]
    reverse=[[] for _ in states]
    for i,row in enumerate(edges):
        for j in row:
            if j>=0:
                reverse[j].append(i)
    useful=set(accepted)
    todo=deque(accepted)
    while todo:
        for i in reverse[todo.popleft()]:
            if i not in useful:
                useful.add(i)
                todo.append(i)
    # Preserve the productive machine: an independent evaluator can audit it.
    acc=set(accepted)
    hist=[]
    dist={0:1}
    for n in range(31):
        hist.append(sum(c for i,c in dist.items() if i in acc))
        new=Counter()
        for i,c in dist.items():
            for j in edges[i]:
                if j in useful:
                    new[j]+=c
        dist=new
    return {'status':'COMPLETE capped DFA','budget':budget,'reachable_states':len(states),
            'productive_states':len(useful),'counts_by_backbone_length':hist,
            'productive_edges':[(i,a,j) for i in sorted(useful) for a,j in enumerate(edges[i]) if j in useful],
            'accepting':accepted,'initial':0}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--budget',type=int,default=6)
    p.add_argument('--limit',type=int,default=100000)
    p.add_argument('--controls',type=int,default=3)
    args=p.parse_args()
    reports=[]
    for name,ports in [('TT3',(6,4,0)),('C3',(2,4,1))]:
        machine=compile_machine(ports)
        controls=0
        for length in range(args.controls+1):
            for word in product(range(8),repeat=length):
                need(evaluate(machine,word)==direct(literal(word,ports))[0],f'direct mismatch {name} {word}')
                controls+=1
        print(json.dumps({'ports':name,'raw_states':machine['raw_states'],
                          'lumped_states':machine['states'],'controls':controls}),flush=True)
        answer=classify(machine,args.budget,args.limit,args.budget==6)
        answer.update(ports=name,raw_states=machine['raw_states'],lumped_states=machine['states'],controls=controls)
        reports.append(answer)
        print(json.dumps({k:v for k,v in answer.items() if k!='productive_edges'}),flush=True)
    Path(__file__).with_name('multi_inport.json').write_text(json.dumps(reports,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    main()
