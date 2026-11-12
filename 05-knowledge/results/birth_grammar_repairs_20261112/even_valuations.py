"""Complete eventual valuation classification of the478 even higher axes.

Odd modes are1,3,5. Bounded branches use exact residues modulo8 after
removing common content. Unbounded branches have a simple9-adic lifting
law, including the one mixed3/5 example. Even modes vanish at high precision.
"""
from common import *
from collections import Counter

def run():
    out=[];hist=Counter();controls=0
    for r in read('affine_births.json.gz',True)['records']:
        u=r['Stirling_coefficients']
        if not any(u[2:]) or u[0]%2:continue
        e=expcoeff(u);odd={i+1:x for i,x in enumerate(e) if i%2==0 and x};even={i+1:x for i,x in enumerate(e) if i%2==1 and x}
        need(set(odd)<={1,3,5},'declared entire odd-mode alphabet')
        g=min(v2(x) for x in odd.values());f={i:c//2**g for i,c in odd.items()}
        vals=[];roots=[]
        if sum(f.values())%2:
            typ='bounded_constant';vals=[g,g]
        else:
            # In this finite bank, each parity's derivative after a=2t+p
            # has exact content8: odd coefficient of3 dominates any5 term.
            need(f.get(3,0)%2==1 and f.get(5,0)%2==0,'simple lifting derivative')
            for p in [0,1]:
                z=sum(c*b**p for b,c in f.items())
                if z%8:vals.append(g+v2(z))
                else:
                    vals.append(None);roots.append(p)
            typ='unbounded' if roots else 'bounded_parity'
        rec={**r,'exponential_coefficients':e,'odd_content_valuation':g,'normalized_odd_modes':f,'type':typ,'eventual_valuations_by_parity':vals,'unbounded_parities':roots}
        if not roots:
            J=max(vals);start=max([1]+[(J+1-v2(c)+v2(b)-1)//v2(b) for b,c in even.items()])
            rec['valuation_start']=start;rec['finite_prefix']=[{'a':a,'U':value(e,a),'v2':v2(value(e,a))} for a in range(1,start)]
            for a in range(start,start+128):need(v2(value(e,a))==vals[a%2],'bounded valuation literal controls');controls+=1
        else:
            lifts=[]
            for p in roots:
                t=0
                for q in range(3,33):
                    mod=2**q
                    need(sum(c*pow(b,2*t+p,mod) for b,c in f.items())%mod==0,'compatible odd-mode root')
                    # Pick a large positive representative so all even modes vanish.
                    period=2**(q-2);a=2*t+p
                    M=g+q
                    while a<max(M,1):a+=period
                    need(sum(c*pow(i+1,a,2**M) for i,c in enumerate(e))%2**M==0,'actual slope high-valuation witness')
                    lifts.append({'precision':M,'parity':p,'t_residue':t,'t_modulus':2**(q-3),'actual_a':a})
                    cand=[t,t+2**(q-3)]
                    good=[z for z in cand if sum(c*pow(b,2*z+p,2**(q+1)) for b,c in f.items())%2**(q+1)==0]
                    need(len(good)==1,'unique next Hensel digit');t=good[0]
            rec['lifting_controls']=lifts
        hist[typ]+=1;out.append(rec)
    need(len(out)==478,'complete even higher-mode input')
    result={'scope':'All478 even higher-mode single birth axes through10. Eventual bounded valuations are proved by fixed odd-mode residues plus an explicit even-mode cutoff. Unbounded valuations follow the simple9/25 lifting criterion verified for every classified root branch. Finite controls are not used in place of that symbolic criterion.','classification':dict(hist),'bounded_literal_controls':controls,'records':out}
    save('even_valuations.json.gz',result);print(json.dumps({'classification':dict(hist),'bounded_literal_controls':controls,'largest_bounded_cutoff':max(r.get('valuation_start',0) for r in out)}),flush=True)
if __name__=='__main__':run()
