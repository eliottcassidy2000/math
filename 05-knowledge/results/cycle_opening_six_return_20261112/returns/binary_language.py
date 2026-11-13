"""Recognize the inherited-ancestry L/K binary and L/T Fibonacci sublanguages."""
from itertools import product
from collections import Counter
import json
import constructor_language as c
K=((2,2,0),(0,0,0),(0,0,1))

def build(word):
    rows=(2,4,1);root=0
    for letter in word:
        rows=c.birth(rows,root)
        if letter=='L':root=len(rows)-1
    return rows,root

def decode(rows,marked,base_order=3):
    current=tuple(rows);ids=list(range(len(rows)));records=[]
    while len(current)>base_order:
        n=len(current);full=(1<<n)-1
        inc=[full^row^(1<<i) for i,row in enumerate(current)]
        low=[i for i,mask in enumerate(inc) if mask.bit_count()==1]
        c.need(len(low)==2,'two indegree-one interface vertices')
        born=[i for i in low if inc[i].bit_length()-1 in low]
        c.need(len(born)==1,'canonical newest birth')
        v=born[0];port=inc[v].bit_length()-1
        records.append((ids[v],ids[port]))
        current=c.a.induced(current,[i for i in range(n) if i!=v]);ids.pop(v)
    records.reverse();out=[]
    for i,(v,port) in enumerate(records):
        chosen=records[i+1][1] if i+1<len(records) else marked
        c.need(chosen in (v,port),'choice must be one of the retained indegree-one pair')
        out.append('L' if chosen==v else 'K')
    return ''.join(out)

def words(cost):
    if cost==0:yield ''
    if cost>=1:
        for w in words(cost-1):yield w+'L'
    if cost>=2:
        for w in words(cost-2):yield w+'T'

def hist(counter):return [dict(H=h,classes=count) for h,count in sorted(counter.items())]

def main():
    c.need(c.mm(c.L,K)==c.T,'T is chronological K then L')
    records=[];checks=0;literal=0
    seed=(2,0,1)
    for cost in range(1,11):
        ordinary={};states=Counter()
        for letters in product('LK',repeat=cost):
            word=''.join(letters);graph,r=build(word)
            c.need(decode(graph,r)==word,'binary rooted history recovered intrinsically')
            reverse=list(reversed(range(len(graph))));relabeled=c.a.induced(graph,reverse)
            c.need(decode(relabeled,len(graph)-1-r)==word,'binary decoding survives relabeling')
            v=seed
            for letter in word:v=c.mv(c.L if letter=='L' else K,v)
            states[v]+=1
            if cost<=5:
                c.need(c.state(graph,r)==v,'literal full Hamiltonian interface');literal+=1
            key=word[:-1]
            if key in ordinary:c.need(ordinary[key][0]==graph,'last marking choice leaves same graph')
            ordinary[key]=(graph,v[0]+v[2]);checks+=1
        c.need(len(ordinary)==2**(cost-1),'binary ordinary family census')
        hcounts=Counter(value[1] for value in ordinary.values())
        rooted_H=Counter()
        for v,mult in states.items():rooted_H[v[0]+v[2]]+=mult
        c.need(all(rooted_H[h]==2*number for h,number in hcounts.items()),'marked-to-ordinary exact factor two')
        fibonacci_H=Counter();fibonacci_words=[]
        for word in words(cost):
            expanded=word.replace('T','KL');graph,r=build(expanded)
            c.need(decode(graph,r)==expanded and expanded.endswith('L') and 'KK' not in expanded,'Fibonacci sublanguage code')
            v=c.mv(c.word_matrix(word),seed);fibonacci_H[v[0]+v[2]]+=1;fibonacci_words.append(word)
        f0,f1=0,1
        for _ in range(cost+1):f0,f1=f1,f0+f1
        c.need(len(fibonacci_words)==f0,'Fibonacci cost census')
        records.append(dict(added_vertices=cost,rooted_binary_words=2**cost,ordinary_binary_classes=len(ordinary),binary_H_distribution=hist(hcounts),binary_interface_states=len(states),Fibonacci_classes=len(fibonacci_words),Fibonacci_H_distribution=hist(fibonacci_H)))
    fermat_checks=[];value=seed
    for births in range(17):
        c.need(value==(2**(births+1),0,1),'closed all-K branch')
        graph,r=build('K'*births)
        order=[2,0]+list(reversed(range(3,births+3)))+[1]
        book=c.a.induced(graph,order);N=births+3
        expected=(2,sum(1<<j for j in range(2,N)))+tuple(1+sum(1<<j for j in range(i+1,N)) for i in range(2,N))
        c.need(book==expected,'explicit C3[1,1,TT] book isomorphism')
        fermat_checks.append(dict(births=births,H=value[0]+value[2],module_order=births+1))
        value=c.mv(K,value)
    core_checks=[]
    for word in ('LL','KL'):
        graph,r=build(word)
        matches=[name for name in ('S5a','S5b') if c.a.iso(graph,c.a.CORES[name]) is not None]
        c.need(len(matches)==1,'two strong H9 five-core representatives')
        core_checks.append(dict(word=word,rows=graph,root=r,core=matches[0],state=c.state(graph,r)))
    c.need(len({x['core'] for x in core_checks})==2,'both distinct H9 cores occur')
    out=dict(status='VERIFIED restricted-language corollary of inherited adaptive ancestry; exact coordinate and graph controls',inheritance='adaptive_metagraph_reformulation_20261103/ANCESTRY.md A1-A3,A5',max_added_vertices=10,binary_rooted_checks=checks,literal_interface_checks=literal,K=K,all_K_book_controls=fermat_checks,H9_cores=core_checks,records=records)
    (c.HERE/'binary_language.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('PASS',checks,'binary rooted histories,',literal,'literal interface checks; Fibonacci restrictions and both H9 cores')
if __name__=='__main__':main()
