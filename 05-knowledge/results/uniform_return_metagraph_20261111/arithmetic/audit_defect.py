"""Independent literal-path audit of the root lane's positive defect grammar.

Own vertex labels and first-vertex DP. Direct run assignments check the four
positive skeletons without importing their transition table.
"""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import importlib.util,json,random

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'defect/two_inport.py'
RECEIPT=HERE.parent/'defect/two_inport.json'
spec=importlib.util.spec_from_file_location('producer',SOURCE)
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)


def check(ok,message):
    if not ok:raise RuntimeError(message)


def graph(word):
    m=len(word);a=m;b=m+1;r=m+2
    rows=[sum(1<<j for j in range(i+1,m)) for i in range(m)]+[(1<<b)|(1<<r),1<<r,(1<<m)-1]
    for i,x in enumerate(word):
        for v,bit in ((a,1),(b,2)):
            if x&bit:rows[i]|=1<<v
            else:rows[v]|=1<<i
    return rows


def starts(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for mask in range(1,1<<n):
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                dp[mask][v]=1 if not rest else sum(dp[rest][w] for w in range(n) if rows[v]>>w&1)
    return dp


def literal(word):
    rows=graph(word);n=len(rows);full=(1<<n)-1;r=n-1
    forward=starts(rows)
    backward=starts([sum(((rows[w]>>v)&1)<<w for w in range(n)) for v in range(n)])
    H=[1]+[sum(v) for v in forward[1:]]
    f=sum(backward[mask][r]*H[full^mask] for mask in range(1,full+1) if mask>>r&1)
    q=H[full^(1<<r)];s=forward[full][r]
    # No reliance on producer's two-run states: build literal path words.
    a=n-3;b=n-2;m=len(word)
    legal=lambda seq:all(rows[u]>>v&1 for u,v in zip(seq,seq[1:]))
    start_ports=0;two_bins=0
    for mask in range(1<<m):
        left=[i for i in range(m) if mask>>i&1]
        right=[i for i in range(m) if not(mask>>i&1)]
        start_ports+=legal([a]+left+[b]+right)+legal([b]+left+[a]+right)
        two_bins+=(legal(left+[a]) and legal([b]+right))+(legal(left+[b]) and legal([a]+right))
    check(start_ports==q-s,'q-s positive paths')
    check(start_ports+two_bins==f-H[full],'positive defect identity')
    return f-H[full],start_ports,two_bins


def main():
    tested=0
    hist={}
    for m in range(5):
        counts=Counter()
        for w in product(range(4),repeat=m):
            b,_,_=literal(w)
            check(b==p.automaton_value(w),f'positive automaton mismatch{w}')
            tested+=1
            if m and w[-1] and any(x&1 for x in w):counts[b]+=1
        hist[m]=dict(counts)
    rng=random.Random(711141425)
    long_checks=0
    for m in range(5,10):
        for _ in range(10):
            w=tuple(rng.randrange(4) for j in range(m))
            b,_,_=literal(w)
            check(b==p.automaton_value(w),'long direct path assignment')
            long_checks+=1
    data=json.loads(RECEIPT.read_text())
    grammar={}
    for aut in data['automata']:
        budget=aut['budget']
        check(aut['acyclic_after_removing_self_loops'],'productive cycle certificate')
        coeff=[0]*41
        for regex in aut['exact_disjoint_path_regex']:
            stars=[x for x in regex if x.endswith('*')]
            check(len(stars)<=1,'one-star series scope')
            length=len(regex)-len(stars)
            if stars:
                for n in range(length,41):coeff[n]+=1
            elif length<=40:coeff[length]+=1
        check(coeff==aut['counts_by_backbone_length_0_through_40'],'regex coefficient agreement')
        if budget==2:expect=[0,1,3]+[5]*38
        elif budget==3:expect=[0]*41
        elif budget==4:expect=[0,1,1,6]+[16]*37
        elif budget==5:expect=[0,0,2]+[0]*38
        else:expect=coeff
        check(coeff==expect,f'claimed all-length series{budget}')
        grammar[budget]={'regex_paths':len(aut['exact_disjoint_path_regex']),
                         'star_paths':sum(any(s.endswith('*') for s in r) for r in aut['exact_disjoint_path_regex']),
                         'first10_coefficients':coeff[:10]}
    check(p.automaton_value((1,))==4,'nonneutral b4 singleton')
    check(p.automaton_value((2,1))==4,'nonneutral b4 length2')
    check(p.automaton_value((1,2))==p.automaton_value((2,3))==5,'nonneutral b5 words')
    report={'status':'PASS','complete_literal_words_lengths_0_through4':tested,
        'additional_direct_words_lengths5_through9':long_checks,
        'scope':'Independent literal graphs, first-vertex DP and direct run assignments; proof review covers regularity/Burnside and all-length regex interpretation, not a second full DFA discovery.',
        'reviewed_sha256':{str(x.name):sha256(x.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for x in [SOURCE,RECEIPT,HERE.parent/'defect/PROOF.md']},
        'grammar':grammar}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (HERE/'defect_audit.json').write_text(out,encoding='utf-8')
    print(out)


if __name__=='__main__':main()
