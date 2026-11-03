"""Concrete guard charts and black-line compilation; integer registers stay external."""
from itertools import product,combinations
from math import lcm
from pathlib import Path
import json
from normal_forms import FRACTIONS,machine,valuations,need
from bank import VALUES,subset_h,converse
HERE=Path(__file__).resolve().parent

def reverse(A,a,b):
    B=list(A);B[a]^=1<<b;B[b]^=1<<a;return tuple(B)

def native(A,P):
    need(sorted(P)==list(range(len(A))) and all(A[a]>>b&1 for a,b in zip(P,P[1:])), 'native path guard')
    kept={tuple(sorted(e)) for e in zip(P,P[1:])};B=A
    for a,b in combinations(range(len(A)),2):
        if (a,b) not in kept:B=reverse(B,a,b)
    return B

def guard(n):return tuple(int(n%b==0) for a,b in FRACTIONS)

def marked_chart(bits):
    n=6;A=[0]*n;chords={(i,i+2):i for i in range(4)}
    for i,j in combinations(range(n),2):
        forward=1 if j==i+1 else bits[chords[(i,j)]] if (i,j) in chords else bits[0]
        if forward:A[i]|=1<<j
        else:A[j]|=1<<i
    return tuple(A)

def normalized_score(A):
    scores=sorted(row.bit_count() for row in A)
    return min(tuple(scores),tuple(sorted(len(A)-1-v for v in scores)))

def main():
    feasible={}
    for bits in product((0,1),repeat=4):
        least=lcm(*(b for bit,(a,b) in zip(bits,FRACTIONS) if bit))
        if guard(least)==bits:feasible[bits]=least
    compatible=sorted(bits for bits in feasible if tuple(1-b for b in bits) in feasible)
    need(len(feasible)==13 and len(compatible)==10,'divisibility-compatible native charts')
    representatives=sorted({min(b,tuple(1-x for x in b)) for b in compatible})
    quotient_edges=set();lift_edges=[]
    for a,b in combinations(compatible,2):
        if sum(x!=y for x,y in zip(a,b))!=1:continue
        lift_edges.append((a,b))
        x=min(a,tuple(1-t for t in a));y=min(b,tuple(1-t for t in b))
        quotient_edges.add(tuple(sorted((representatives.index(x),representatives.index(y)))))
    degrees=[sum(i in e for e in quotient_edges) for i in range(5)]
    need(sorted(degrees)==[2,2,2,3,3] and len(quotient_edges)==6 and len(lift_edges)==12,'K2,3 quotient')
    need(all(sum(representatives[i])%2!=sum(representatives[j])%2 for i,j in quotient_edges),'quotient bipartition')
    for bits in product((0,1),repeat=4):
        need(native(marked_chart(bits),tuple(range(6)))==marked_chart(tuple(1-b for b in bits)),'literal chart complement')
    need(guard(25)==guard(625) and guard(27)!=guard(675),'snapshot loses future state')
    forward={(0,1),(0,3),(0,4),(0,5),(1,2),(2,3),(2,6),(3,4),(3,6),(4,5),(4,6),(5,6)}
    W=[0]*7
    for a,b in combinations(range(7),2):
        if (a,b) in forward:W[a]|=1<<b
        else:W[b]|=1<<a
    W=tuple(W);cycle=[(i,(i+1)%7) for i in range(7)]
    paths=[tuple((g+1+i)%7 for i in range(7)) for g in range(7)]
    leaves=[reverse(W,a,b) for a,b in cycle];Z=native(leaves[0],paths[0])
    need(len({normalized_score(A) for A in [Z]+leaves})==8,'eight distinct merged classes')
    routes=[]
    for g,h in product(range(7),repeat=2):
        need(native(leaves[g],paths[g])==Z and native(Z,paths[h])==leaves[h],'actual two-line route')
        routes.append({'from':g,'to':h,'paths':[paths[g],paths[h]]})
    # A spoke could be blue only if reversal of its marked path were an anti-automorphism.
    for g,P in enumerate(paths):
        anti={P[i]:P[-1-i] for i in range(7)}
        is_anti=all(not (Z[i]>>j&1) or (Z[anti[j]]>>anti[i]&1) for i,j in product(range(7),repeat=2) if i!=j)
        need(not is_anti,'black spoke')
    executions=[]
    for n in (1,25,625,49,21,35,3**4*5**7*7**9,101*3**2*5**4*7**3):
        end,trace=machine(n);state=n;g=6;clock=0;history=[]
        for rule in trace:
            selected=next(i for i,(a,b) in enumerate(FRACTIONS) if state%b==0)
            need(selected==rule-1,'priority')
            dest=selected;mid=native(leaves[g],paths[g]);out=native(mid,paths[dest])
            need(out==leaves[dest],'native compiler endpoint')
            a,b=FRACTIONS[selected];new=state//b*a
            history.append({'integer_before':state,'integer_after':new,'leaf_before':g,'leaf_after':dest,
                            'rule':rule,'paths':[paths[g],paths[dest]],'native_clock':[clock,clock+1,clock+2]})
            state=new;g=dest;clock+=2
        need(state==end and clock==2*len(trace),'exact execution clock')
        executions.append({'input':n,'output':end,'steps':history,'native_steps':clock})
    bank=json.loads((HERE/'bank.json').read_text());table=bank['transition']
    need(len(table)==7,'seven primitive logical states')
    for g,row in enumerate(table):
        for h in row:need(native(native(leaves[g],paths[g]),paths[h])==leaves[h],'every primitive selector transition compiled')
    out={'status':'PROVED_SCOPED','feasible_guards':[{'bits':b,'least_integer':n} for b,n in sorted(feasible.items())],
         'complement_compatible':compatible,'quotient_representatives':representatives,
         'quotient_edges':sorted(quotient_edges),'lift_edges':lift_edges,
         'snapshot_hostile':[{'input':n,'guard':guard(n),'successor':n//25*27,'next_guard':guard(n//25*27)} for n in (25,625)],
         'hub':{'W':W,'center':Z,'leaves':leaves,'paths':paths,'scores':[normalized_score(A) for A in [Z]+leaves],
                'H':[subset_h(A)[-1] for A in [Z]+leaves]},'all_routes':routes,'executions':executions,
         'primitive_selector_compiled_transitions':70,'primitive_alphabet':VALUES}
    (HERE/'controllers.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: 13 feasible / 10 complement-compatible guards; K2,3; 49 native routes; 70 selector transitions')
if __name__=='__main__':main()
