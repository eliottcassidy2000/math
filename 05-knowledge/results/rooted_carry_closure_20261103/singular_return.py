"""The incoming moving guard is NR: actual rank-one reset completes the residue monoid."""
from collections import Counter
from itertools import product
import json,hashlib
from observer import HERE,junction,need,old,observer,root_profile
from returns import cycle
from field_clock import M3,M5,ID,mul,mv,act,pair_point
from golden_transport import sealed,virtual_cycle,apply as cover_step,obs,phase
Q=(3,0,1,0)

def nr(T,root,ports,b,t,literal=False):
    pp=list(ports);slot=pp.index(root);fixed=pp[1-slot];r=len(T)
    U=junction.ear(T,1<<root);pp[slot]=r
    steps=[{'root':root,'newborn':r,'ports':tuple(pp),'discard':'old guard'}]
    last=len(U);U=junction.ear(U,1<<r)
    steps.append({'root':r,'newborn':last,'ports':tuple(pp),'discard':'newborn'})
    pp=tuple(pp);need(fixed in pp,'fixed other observer port survives the singular reset')
    if literal:need(junction.profile(U,r)==(3*b,b,b,3*b),'literal rank-one balanced reset')
    need(junction.strong(U),'strong singular reset')
    return U,r,pp,3*b,b,{'command':'N','steps':steps,'rows':U,'root':r,'ports':pp,'pair':[3*b,b]}

def main():
    generators={'3':tuple(x%7 for x in M3),'5':tuple(x%7 for x in M5),'N':Q}
    words={ID:''};queue=[ID]
    for T in queue:
        for letter,G in generators.items():
            U=mul(G,T,7)
            if U not in words:words[U]=words[T]+letter;queue.append(U)
    need(set(words)==set(product(range(7),repeat=4)),'full M2(F7), including all singular matrices')
    ranks=Counter(0 if not any(T) else 2 if (T[0]*T[3]-T[1]*T[2])%7 else 1 for T in words)
    need(ranks=={0:1,1:384,2:2016},'exact rank strata')
    receipts=json.loads((HERE/'returns.json').read_text());seals={x['k']:x for x in receipts['seals']}
    v=tuple(seals[2]['pair']);w=tuple(seals[3]['pair']);point=pair_point(v,w)
    witnessed={}
    for A,word in words.items():witnessed.setdefault(act(A,point),{'word':word,'matrix':A})
    need(set(witnessed)==set(range(51)),'all fifty-one classes from actual paired sources')
    all_actions=tuple(words)
    coarse={};fine={}
    for p in range(51):
        c=tuple(act(A,p) in (49,50) for A in all_actions)
        f=tuple((act(A,p) in (49,50),act(A,p) in (7,50)) for A in all_actions)
        coarse.setdefault(c,[]).append(p);fine.setdefault(f,[]).append(p)
    need(sorted(map(len,coarse.values()))==[1]*9+[42] and len(fine)==51,'10/51 quotients persist under singular actions')
    zero_word=words[(0,0,0,0)];literal=[]
    wa,wb=old.chart(3,4);base=old.source(wa,wb,False);base_cover=observer(base,(7,8))
    for k in (2,3):
        source=seals[k];T=tuple(source['rows']);root=source['root'];ports=tuple(source['ports']);b,t=source['pair'];events=[]
        CT,CP,CV,slot=sealed(base,base_cover,k)
        need(CT==T and CP==ports and CP[slot]==root,'full observer supplied source alignment')
        # An independently counted first singular step; start over for the shortest zero word.
        nr(T,root,ports,b,t,literal=True)
        for letter in zero_word:
            if letter=='N':
                CT,CP,CV=cover_step(CT,CP,CV,(0,1<<slot,slot))
                CT,CP,CV=cover_step(CT,CP,CV,(0,1<<slot,2))
            else:CT,CP,CV,slot,_,_=virtual_cycle(CT,CP,CV,slot,int(letter))
            if letter=='N':out=nr(T,root,ports,b,t)
            else:out=cycle(T,root,ports,b,t,int(letter))
            T,root,ports,b,t,receipt=out
            need(CT==T and CP==ports and CP[slot]==root and root_profile(CV,slot)==(b,t,t,b),'full-cover and actual native singular-word transport')
            events.append({'letter':letter,'order':len(T),'root':root,'ports':ports,'pair':[b,t],
                           'receipt_sha256':hashlib.sha256(json.dumps(receipt,sort_keys=True).encode()).hexdigest()})
        need(b%7==t%7==0,'actual zero-phase source with positive integer counts')
        literal.append({'initial_k':k,'word':zero_word,'events':events,'final_rows':T,'root':root,'ports':ports,'pair':[b,t],
                        'full_golden_observer':obs(CV,slot),'full_golden_phase':phase(obs(CV,slot))})
    out={'status':'PROVED_ALL_WORDS_WITH_FINITE_EXACT_MONOID_CERTIFICATE','singular_map':Q,'matrix_monoid_size':2401,
         'rank_counts':dict(ranks),'zero_word':zero_word,'source_pair':list(map(list,(v,w))),
         'projective_witnesses':[{'point':p,**witnessed[p]} for p in range(51)],'coarse_classes':list(coarse.values()),
         'fine_classes':len(fine),'literal_zero_phase_sources':literal,
         'scope':'raw actual balanced-source macros; primitive arithmetic guards intentionally not asserted after singular reset'}
    (HERE/'singular_return.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: actual NR reset, full 2401-matrix monoid, all51 strata, and literal zero-phase sources; zero word',zero_word)
if __name__=='__main__':main()
