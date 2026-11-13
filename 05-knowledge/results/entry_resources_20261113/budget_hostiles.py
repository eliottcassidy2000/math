"""Higher-order budget obstructions, and local-density controls."""
from common import *
from itertools import combinations
from math import comb
import random
from independent_checks import signed_cut

def run():
    paley=[]
    for q in [3,7,11,19]:
        squares={x*x%q for x in range(1,q)}
        rr=[sum(1<<j for j in range(q) if (j-i)%q in squares) for i in range(q)]
        h,u,w,marks=tables(np.array(rr,dtype=np.int64));weights={int(w[i,j]) for i in range(q) for j in range(i+1,q)}
        need(len(weights)==1,'unordered-pair transitivity produces constant weights')
        weight=next(iter(weights));r=next(k for k in range(2,q+1) if comb(k,2)*weight>int(h)-1)
        paley.append({'q':q,'H':int(h),'pair_weight':weight,'minimal_nonface_size':r,'unit_markings':len(marks)})
    rng=random.Random(131313);hostile=None;tested=0
    for n in [9,10]:
        for attempt in range(2000):
            rr=labelled(n,rng.randrange(1<<(n*(n-1)//2)));h,u,w,_=tables(np.array(rr,dtype=np.int64));tested+=1
            for subset in combinations(range(n),4):
                cost=sum(int(w[a,b]) for a,b in combinations(subset,2))
                if cost<=int(h)-1:continue
                if all(sum(int(w[a,b]) for a,b in combinations(t,2))<=int(h)-1 for t in combinations(subset,3)):
                    mask=sum(1<<v for v in subset);checks=[]
                    for s in [mask]+[mask^(1<<v) for v in subset]:
                        d=int(signed_cut(np.array(rr,dtype=np.int64),s));checks.append([s,d])
                    need(checks[0][1]<1 and all(d>=1 for s,d in checks[1:]),'independent minimal four-way obstruction')
                    hostile={'n':n,'attempt':attempt,'rows':rr,'H':int(h),'selected':mask,'subset':subset,'signed_cut_checks':checks};break
            if hostile:break
        if hostile:break
    residues=[]
    for r in range(1,13):
        mod=1<<r;values={pow(3,a,mod)*pow(5,b,mod)%mod for a in [0,1] for b in range(1<<max(0,r-2))}
        need(values==set(range(1,mod,2)),'3 and5 generate every odd residue')
        residues.append({'r':r,'odd_residues':len(values)})
    save('budget_hostiles',{'scope':'Constant pair weights on Paley3,7,11,19; fixed-seed search for a minimal four-vertex budget nonface after the complete through8 census. Dyadic residue coverage is an arithmetic control, not global Hamiltonian coverage.','paley':paley,'sampled_tournaments':tested,'four_way_hostile':hostile,'dyadic_controls':residues})
    print('Paley',paley,'hostile',hostile,'sampled',tested,flush=True)
if __name__=='__main__':run()
