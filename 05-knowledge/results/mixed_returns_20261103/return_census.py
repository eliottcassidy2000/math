"""Exhaustive short returns, retaining full words and bulk projective data."""
from fractions import Fraction as F
from collections import Counter
import json
from support import *

def interval(q):
    vals=[F(q[0][1]+r*q[1][1],q[0][0]+r*q[1][0]) for r in (0,3)]
    return min(vals),max(vals)

def main():
    states=[(T,())];returns=[];census=[]
    inputs=[vector(w) for w in SEEDS]+list(REPAIRS)
    ratios=[F(bulk(v)[1],bulk(v)[0]) for v in inputs]
    for depth in range(1,7):
        nxt=[]
        for v,w in states:
            for l,p in OLD_MATS.items():
                z=step(v,p)
                if min(z)>=0 or max(z)<=0:continue
                word=w+(l,);nxt.append((z,word))
                if z[2]==0 and z[0]==-2*z[1]:
                    mat=matrix(word);q=quotient(mat);lo,hi=interval(q)
                    returns.append({'word':word,'gain':-z[1],'cost':sum(COSTS[l] for l in word),
                        'matrix':mat,'bulk_matrix':q,'interval':[str(lo),str(hi)],
                        'candidate_seeds':[i for i,r in enumerate(ratios) if lo<=r<=hi]})
        states=nxt
        census.append({'length':depth,'mixed_words':len(states),'mixed_rows':len({z for z,w in states}),
            'return_words':sum(len(r['word'])==depth for r in returns)})
    candidates=[r for r in returns if r['candidate_seeds']]
    result={'status':'PASS','alphabet':'fixed 65; all words of one through six letters; only irreversible one-sided error rows pruned',
        'census':census,'returns':returns,'seed_states':inputs,'seed_ratios':[str(r) for r in ratios],
        'candidate_returns':candidates,'all_gains':sorted({r['gain'] for r in returns})}
    (HERE/'return_census.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'returns':len(returns),'census':census,'candidate_count':len(candidates),
        'first_candidates':[(r['word'],r['candidate_seeds']) for r in candidates[:12]],'gains':result['all_gains']}))
if __name__=='__main__':main()
