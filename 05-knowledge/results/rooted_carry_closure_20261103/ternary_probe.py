"""Complete bounded six-command descendant search from the two named H497 sources."""
from collections import Counter
import json
from observer import HERE,old,observer,load,need
live=load('inherited_live','05-knowledge/results/adaptive_state_induction_20261103/live_port.py')

def vector(T,ports,states):
    out=Counter()
    for (sk,m),c in observer(T,ports).items():
        state=tuple(sorted(sk+((-1,-1),)*m))
        if sum(b==-1 for a,b in state)<=1:out[state]+=c
    need(set(out)<=set(states),'exact cap-one projection to supplied-source carrier')
    return tuple(out[s] for s in states)

def main():
    states,matrices,_,obs=live.compile_carrier();w,u=old.chart(3,4);records=[]
    for forward in (True,False):
        T=old.source(w,u,forward);ports=(7,8);v=vector(T,ports,states)
        words={v:()};queue=[v];tested=0;counts=set();max_depth=0
        for v in queue:
            word=words[v];h=live.output(v,obs);counts.add(h);max_depth=max(max_depth,len(word))
            U=T;pp=list(ports)
            for letter in word:
                r,drop=divmod(letter,3);x=len(U);U=live.birth(U,pp[r])
                if drop<2:pp[drop]=x
            need(old.hp(U)==h,'literal representative Hamilton count')
            for letter,M in enumerate(matrices):
                W=live.advance(v,M);hh=live.output(W,obs);tested+=1
                need(hh>h,'strict actual count growth')
                if hh<=1495 and W not in words:words[W]=word+(letter,);queue.append(W)
        need(not counts&{1493,1495},'declared six-command language misses both H497 ternary offsets')
        records.append({'orientation':'plus' if forward else 'minus','initial_rows':T,'ports':ports,'cutoff':1495,
                        'distinct_observers':len(words),'transitions':tested,'maximum_depth':max_depth,'attained_H':sorted(counts),
                        'missing_targets':[1493,1495],'words':[list(word) for word in words.values()]})
    (HERE/'ternary_probe.json').write_text(json.dumps({'status':'FINITE_EXACT_COMPLETE_BOUNDED_LANGUAGE','controls':records},indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: exact named H497 six-command cutoff1495; both ternary targets require an additional operation')
if __name__=='__main__':main()
