"""Exact interface and wheel controls; standard library, explicit gates."""
from itertools import product
from math import gcd, isqrt
from collections import Counter
import json

DEAD=-2
FREE=-1


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def classify(mask):
    if not mask:
        return DEAD
    if mask & (mask-1):
        return FREE
    return mask.bit_length()-1


def step(state, mask):
    if state == DEAD:
        return DEAD
    if state >= 0:
        mask &= ~(1 << state)
    return classify(mask)


def path_machine(lists):
    state=FREE
    for mask in lists:
        state=step(state,mask)
    return state != DEAD


def cycle_machine(lists,q):
    for a in range(q):
        if not (lists[0]>>a & 1):
            continue
        state=a
        for mask in lists[1:]:
            state=step(state,mask)
        if state==FREE or (state>=0 and state!=a):
            return True
    return False


def proper_words(q,n,cycle=False):
    return [word for word in product(range(q),repeat=n)
            if all(a!=b for a,b in zip(word,word[1:])) and
            (not cycle or word[-1]!=word[0])]


def literal(lists,words):
    return any(all(mask>>colour&1 for mask,colour in zip(lists,word))
               for word in words)


def obstruction(lists):
    if 0 in lists:
        return True
    for start,left in enumerate(lists):
        if left.bit_count()!=1:
            continue
        forced=left.bit_length()-1
        for right in range(start+1,len(lists)):
            mask=lists[right]
            if mask.bit_count()==1:
                if mask == 1<<forced:
                    return True
                break
            if mask.bit_count()!=2 or not (mask>>forced&1):
                break
            forced=(mask^(1<<forced)).bit_length()-1
    return False


def prime(n):
    return n>=2 and all(n%d for d in range(2,isqrt(n)+1))


def main():
    transition_count=0; distinguish=0
    for q in range(2,9):
        for previous in range(1<<q):
            for allowed in range(1<<q):
                actual=0
                for b in range(q):
                    if (allowed>>b&1) and any((previous>>a&1) and a!=b for a in range(q)):
                        actual |= 1<<b
                need(step(classify(previous),allowed)==classify(actual),
                     ('transition mismatch',q,previous,allowed))
                transition_count+=1
        states=[DEAD,FREE]+list(range(q))
        reached={step(FREE,mask) for mask in range(1<<q)}
        need(reached==set(states),('unreachable state',q))
        for i,a in enumerate(states):
            for b in states[i+1:]:
                if DEAD in (a,b):
                    need((a!=DEAD)!=(b!=DEAD),'empty suffix distinction')
                else:
                    forced=a if a>=0 else b
                    need((step(a,1<<forced)!=DEAD)!=(step(b,1<<forced)!=DEAD),
                         ('singleton suffix distinction',q,a,b))
                distinguish+=1
    path_count=0
    for q,maximum in ((2,6),(3,4)):
        for n in range(1,maximum+1):
            words=proper_words(q,n)
            for lists in product(range(1<<q),repeat=n):
                actual=literal(lists,words)
                need(path_machine(lists)==actual,('path',q,lists))
                need(obstruction(lists)==(not actual),('chain certificate',q,lists))
                path_count+=1
    cycles=Counter()
    for q,n in ((3,6),(4,4)):
        words=proper_words(q,n,True)
        masks=[m for m in range(1,1<<q) if m.bit_count()<=2]
        for lists in product(masks,repeat=n):
            actual=literal(lists,words)
            need(cycle_machine(lists,q)==actual,('cycle',q,lists))
            cycles[(q,n,actual)]+=1
    bad=(1,3,6,5,3,3)
    good=(1,3,6,5,6,3)
    need(not cycle_machine(bad,5) and cycle_machine(good,5),'anchor hostile')
    need([m.bit_count() for m in bad]==[m.bit_count() for m in good],
         'same-size hostile')
    for lists in (bad,good):
        forbidden=[31 & ~(mask|16) for mask in lists]
        need(sum(m.bit_count() for m in forbidden)==13,'twenty-vertex carrier size')
        need(all((31 & ~(f|16))==m for f,m in zip(forbidden,lists)),
             'physical list reconstruction')
    first=[sum(all(m>>c&1 for m,c in zip((3,),word)) for word in proper_words(3,1)
               if word[-1]==j) for j in range(3)]
    second=[sum(all(m>>c&1 for m,c in zip((7,3),word)) for word in proper_words(3,2)
                if word[-1]==j) for j in range(3)]
    need(first==[1,1,0] and second==[2,2,0],'weighted hostile')
    need(sum(first)-first[0]==1 and sum(second)-second[0]==2,'singleton count distinction')
    p=30030
    need(p==2*3*5*7*11*13,'primorial')
    need(all(gcd(p,d)>1 for d in range(2,17)),'first wheel gap')
    named=[30013,30047,29741,30319,30011,30029]
    need(all(prime(n) for n in named),'named endpoint primality')
    need(30049==151*199 and 30007==37*811 and 30053==41*733,'endpoint factors')
    need(30031==59*509,'future-query bridge factor')
    offsets={str(sign):next(d for d in range(2,290) if prime(p+sign*d)) for sign in (-1,1)}
    need(set(offsets.values())=={17},'first prime endpoint')
    prime_endpoints=0
    for d in range(2,289):
        for sign in (-1,1):
            if prime(p+sign*d):
                prime_endpoints+=1
                need(gcd(p,d)==1 and prime(d),'short-window theorem')
    need(not prime(289) and prime(p-289) and prime(p+289),'sharp square boundary')
    a=(30013,30047); b=(30029,30031)
    need(sum(a)==sum(b)==60060,'equal owner degree')
    need(all(sum(gcd(c,d) for c in a)==sum(gcd(c,d) for c in b)==2
             for d in range(1,59,2)), 'bounded cycle-query equivalence')
    need(sum(gcd(c,59) for c in a)==2 and sum(gcd(c,59) for c in b)==60,
         'first split at59')
    print(json.dumps({'status':'PASS','transition_identities':transition_count,
        'distinguishing_state_pairs':distinguish,'path_list_words':path_count,
        'cycle_list_words':[{'q':q,'length':n,'colourable':yes,'count':count}
                            for (q,n,yes),count in sorted(cycles.items())],
        'states_at_q4_q5_q6':[6,7,8],'physical_hostile_vertices':20,
        'prime_endpoints_in_open_square_window':prime_endpoints,
        'first_offsets':offsets,'bounded_query_bridge':{'degree':60060,'last_equal_odd_query':57,
                                                     'first_different_query':59,'responses':[2,60]},
        'scope':'exact finite controls; no HN or Fortune lower-bound/conjecture closure'},indent=2))


if __name__=='__main__':
    main()
