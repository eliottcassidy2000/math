#!/usr/bin/env python3
"""Error-refined SCC series check on an explicit bounded actual-core alphabet.
This is a formula hostile control, not an all-tournament census.
"""
from pathlib import Path
from collections import Counter
import argparse,json
import check

def need(ok,label):
    if not ok:raise RuntimeError(label)
def multiply(P,Q,N):
    out=Counter()
    for (n,e),v in P.items():
        for (m,f),w in Q.items():
            if n+m<=N:out[n+m,e+f]+=v*w
    return +out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='series-normal.json');args=ap.parse_args();N=16
    data=json.loads((check.HERE/'normal.json').read_text())
    asymmetric=tuple(next(r for r in data['representatives'] if r['source']=='strong_five' and not r['self_converse'])['rows'])
    BAD7=(74,84,97,22,37,11,56)
    cores=[('one',(0,)),('three',check.cyc(3)),('bad_seven',BAD7),('ten',asymmetric),('ten_dual',check.reverse(asymmetric))]
    alphabet=[]
    for name,T in cores:
        need(check.strong(T),'every alphabet symbol is an actual strong tournament')
        e=check.error(T)[0];dual=[i for i,(_,U) in enumerate(cores) if check.mappings(check.reverse(T),U,True)]
        need(len(dual)==1,'alphabet contains exactly one converse partner per core')
        alphabet.append({'name':name,'order':len(T),'error':e,'dual':dual[0]})
    need(all(alphabet[a['dual']]['dual']==i for i,a in enumerate(alphabet)),'core converse involution')
    S=Counter((a['order'],a['error']) for a in alphabet)
    Fs=Counter((a['order'],a['error']) for i,a in enumerate(alphabet) if a['dual']==i)
    A=Counter({(0,0):1});power=Counter({(0,0):1})
    for _ in range(N):power=multiply(power,S,N);A.update(power)
    directA=Counter();directF=Counter()
    def visit(word,n,e):
        directA[n,e]+=1
        if tuple(alphabet[i]['dual'] for i in word[::-1])==word:directF[n,e]+=1
        for i,a in enumerate(alphabet):
            if n+a['order']<=N:visit(word+(i,),n+a['order'],e+a['error'])
    visit((),0,0)
    need(A==directA,'SCC ordinary series equals literal ordered-word enumeration')
    double=Counter({(2*n,2*e):c for (n,e),c in A.items() if 2*n<=N})
    center=Fs+Counter({(0,0):1});F=multiply(center,double,N)
    need(F==directF,'self-converse error series equals literal reversed-dual words')
    # Differentiate by the error marker at t=1 and compare the full formula.
    countsA=Counter();countsF=Counter();errorsA=Counter();errorsF=Counter();errorsS=Counter();errorsFs=Counter()
    for (n,e),c in A.items():countsA[n]+=c;errorsA[n]+=e*c
    for (n,e),c in F.items():countsF[n]+=c;errorsF[n]+=e*c
    for (n,e),c in S.items():errorsS[n]+=e*c
    for (n,e),c in Fs.items():errorsFs[n]+=e*c
    rhs=Counter()
    for a,v in countsA.items():
        for b,w in errorsFs.items():
            if 2*a+b<=N:rhs[2*a+b]+=v*w
        for b,w in countsF.items():
            for c,x in errorsS.items():
                if 2*a+b+2*c<=N:rhs[2*a+b+2*c]+=2*v*w*x
    need(+rhs==+errorsF,'differentiated self-converse error series with factor two')
    incorrect=Counter({(2*n,e):c for (n,e),c in A.items() if 2*n<=N})
    wrong=multiply(center,incorrect,N)
    first_wrong=min(n for n in range(N+1) if any(wrong[n,e]!=F[n,e] for e in {x[1] for x in wrong}|{x[1] for x in F}))
    need(first_wrong==14,'dropping doubled error first fails at paired bad seven')
    # Omitting lower-order centers also fails before the paired-core threshold.
    first_nonstrong_positive=min(n for (n,e),c in (F-center).items() if c and e)
    need(first_nonstrong_positive==9,'odd central positive core already pads at nine')
    out={'status':'PASS','scope':'actual finite core alphabet, all SCC words of total order at most sixteen','alphabet':alphabet,'word_count':sum(directA.values()),'self_converse_word_count':sum(directF.values()),'first_failure_if_t_squared_omitted':first_wrong,'first_nonstrong_positive_selfconverse_order_in_alphabet':first_nonstrong_positive,'selfconverse_positive_error_coefficients':[[n,e,c] for (n,e),c in sorted(F.items()) if e],'selfconverse_weighted_error_coefficients':sorted((n,c) for n,c in errorsF.items() if c)}
    Path(args.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,indent=2,sort_keys=True))
if __name__=='__main__':main()
