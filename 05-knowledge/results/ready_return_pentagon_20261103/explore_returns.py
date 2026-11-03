"""Constructive return probes and exact bounded joint-response languages."""
from pathlib import Path
import sys,json
from collections import Counter,deque
from itertools import product

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'two_entry_funnels_20261103'))
from research import *
from verify import rooted_data,start_paths_at_most_two,joint_pair_word
HERE=Path(__file__).resolve().parent

ALPHABET=((0,0),(0,1),(1,0),(1,1))


def transition(state,letter,cap):
    g1,d1,g2,d2,ap,bp=state;a,b=letter
    v=(a*(g1+d1)+g1+(1-b)*ap,
       (1-a)*(g1+d1)+d1+(1-b)*(1-ap),
       b*(g2+d2)+g2+(1-a)*bp,
       (1-b)*(g2+d2)+d2+(1-a)*(1-bp))
    return tuple(min(cap,x) for x in v)+(a,b)


def machine(target):
    initial=(0,0,0,0,1,1);states=[initial];index={initial:0};rows=[]
    for state in states:
        row=[]
        for letter in ALPHABET:
            w=transition(state,letter,target+1)
            if w not in index:index[w]=len(states);states.append(w)
            row.append(index[w])
        rows.append(row)
    final={i for i,(g1,d1,g2,d2,a,b) in enumerate(states) if (a or b) and a+b+g1+g2==target}
    reverse=[[] for _ in states]
    for i,row in enumerate(rows):
        for j in row:reverse[j].append(i)
    good=set(final);queue=deque(final)
    while queue:
        for i in reverse[queue.popleft()]:
            if i not in good:good.add(i);queue.append(i)
    order=sorted(good);renumber={i:j for j,i in enumerate(order)}
    return [states[i] for i in order],[[renumber.get(j,-1) for j in rows[i]] for i in order],{renumber[i] for i in final},renumber.get(0,-1),len(states)


def finite_tails(target):
    states,rows,final,start,total=machine(target)
    require(start>=0,'productive initial')
    tails=[];stack=[(start,())]
    while stack:
        i,w=stack.pop()
        if i in final:tails.append(w)
        for letter,j in zip(ALPHABET,rows[i]):
            if j<0 or (i==start and letter==(1,1)):continue
            require(len(w)<2*(target.bit_length()-1)+4,'proved positive-energy tail bound')
            stack.append((j,w+(letter,)))
    return dict(target=target,reachable_states=total,productive_states=len(states),
                tails=sorted(tails,key=lambda w:(len(w),w)),max_tail=max(map(len,tails)),
                state_rows=rows,states=states,final=sorted(final),start=start)


def energy(word):
    """Independent closed count by first/last assignments to the two paths."""
    out=sum(word[-1])
    for (a,b),(c,d) in zip(word,word[1:]):out+=a*(1-d)+b*(1-c)
    for j,(a,b) in enumerate(word):
        for i in range(j+1,len(word)):
            c,d=word[i]
            out+=2**(i-j-1)*((1-b)*c+(1-a)*d)
    return out


def coefficients(tail,target):
    T,r=incoming_pair_source(tail)
    W=induced(T,list(range(r)))
    bins=[bin_values(W,port,2) for port in (0,1)]
    a=sum(z[1][2] for z in bins);b=tables(W)[0][-1]-target
    tables.cache_clear()
    return a,b


def main():
    result=[]
    for target in (2,4,6,8,10,12):
        r=finite_tails(target)
        print('target',target,'states',r['productive_states'],'tails',len(r['tails']),'length',r['max_tail'],flush=True)
        if target<=6:
            families=[]
            for tail in r['tails']:
                vals=[];a,b=coefficients(tail,target)
                for k in (1,2,3):
                    T,root=incoming_pair_source(((1,1),)*k+tail)
                    h,s,q,F=rooted_data(T,root)
                    require(F-h==target and s==q and strong(T),'actual neutral defect source')
                    require(h==a*2**k+b,'finite-core all-size Hamilton formula')
                    vals.append(h);tables.cache_clear()
                families.append(dict(tail=tail,A=a,B=b,H_k1_to3=vals))
            r['families']=families
        result.append(r)
    (HERE/'languages.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('C4',result[1]['families'])
    energy_checks=direct_checks=0
    accepted={r['target']:set(r['tails']) for r in result}
    for m in range(1,10):
        for tail in product(ALPHABET,repeat=m-1):
            word=((1,1),)+tail;c=energy(word)
            require(c==joint_pair_word(word),'closed energy vs recurrence')
            k=0
            while k<len(word) and word[k]==(1,1):k+=1
            if word[-1]!=(0,0) and c in accepted:
                require(word[k:] in accepted[c],'independent word belongs to complete atlas')
            if m<=4:
                T,r=incoming_pair_source(word);W=induced(T,list(range(r)))
                require(c==joint_values(W,{0,1},1),'energy vs subset paths')
                tables.cache_clear();direct_checks+=1
            energy_checks+=1
    print('PASS energy',energy_checks,'independent subset',direct_checks)


if __name__=='__main__':main()
