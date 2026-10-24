"""Exact sparse completion of the six-slot alphabet, with binary contexts retained."""
from common import *
from collections import defaultdict,Counter
from itertools import combinations
def components(S,budget):
    unseen=set(S);out=[]
    while unseen:
        q=[min(unseen)];unseen.remove(q[0])
        for a in q:
            hits=[b for b in unseen if (a^b).bit_count()<=budget]
            unseen.difference_update(hits);q.extend(hits)
        out.append(sorted(q))
    return out
class DSU:
    def __init__(self,n):self.p=list(range(n))
    def root(self,a):
        while self.p[a]!=a:self.p[a]=self.p[self.p[a]];a=self.p[a]
        return a
    def join(self,a,b):
        a=self.root(a);b=self.root(b)
        if a==b:return False
        self.p[b]=a;return True
def main():
    rec=inventory();ids=[r['source'] for r in rec]
    lab=[(sum(b<<j for j,b in enumerate(r['source_label'][:4])),r['source_label'][4]) for r in rec]
    costs=[r['cost'] for r in rec]
    def work(m):return sum(c for i,c in enumerate(costs) if m>>i&1)
    def charge(m):
        p=x=0
        for j in range(28):
            if m>>j&1:p^=lab[j][0];x+=lab[j][1]
        return p,x%128
    bank={};comps={};where={};counts={}
    for k in (4,5):
        D=defaultdict(list)
        for I in combinations(range(28),k):
            m=sum(1<<i for i in I);D[charge(m)].append(m)
        bank[k]=D;comps[k]={g:components(S,6) for g,S in D.items()};where[k]={m:(g,j) for g,C in comps[k].items() for j,S in enumerate(C) for m in S}
        counts[k]={'states':sum(map(len,D.values())),'bad_fibres':sum(len(C)>1 for C in comps[k].values()),'extra_components':sum(len(C)-1 for C in comps[k].values())}
    # Every four-margin additional edge has disjoint support8. Group by its
    # intrinsic four-subset charge, not by a higher-margin context.
    candidates=[];cuts=[]
    for g,C in sorted(comps[5].items()):
        if len(C)>1:
            need(len(C)==2,'each margin5 hostile is a two-component cut')
            cuts.append((g,set(C[0]),set(C[1])))
    for g,C in sorted(comps[4].items()):
        for i,j in combinations(range(len(C)),2):
            for a in C[i]:
                for b in C[j]:
                    need(not a&b,'six-disconnected four-subsets are disjoint')
                    covers=[]
                    for h,(cg,L,R) in enumerate(cuts):
                        for slot in range(28):
                            bit=1<<slot
                            if (a|b)&bit:continue
                            if ((a|bit) in L and (b|bit) in R) or ((b|bit) in L and (a|bit) in R):covers.append(h);break
                    candidates.append((a,b,g,i,j,tuple(covers)))
    all_candidate_count=len(candidates)
    candidates=[c for c in candidates if work(c[0]|c[1])<=22]
    ds={g:DSU(len(C)) for g,C in comps[4].items()};chosen=[];uncovered=set(range(len(cuts)))
    while uncovered:
        possible=[c for c in candidates if ds[c[2]].root(c[3])!=ds[c[2]].root(c[4]) and uncovered&set(c[5])]
        need(possible,'margin5 cuts admit a forest-compatible four-margin bridge')
        c=min(possible,key=lambda c:(-len(uncovered&set(c[5])),work(c[0]|c[1]),c[0],c[1]));chosen.append(c);ds[c[2]].join(c[3],c[4]);uncovered.difference_update(c[5])
    for c in sorted(candidates,key=lambda c:(work(c[0]|c[1]),c[0],c[1])):
        if ds[c[2]].join(c[3],c[4]):chosen.append(c)
    need(len(chosen)==counts[4]['extra_components']==220,'sharp count of extra eight-slot move types')
    clock_witness=sorted(bank[4][6,123])
    need(clock_witness==[12864,69272576] and work(clock_witness[0]^clock_witness[1])==22,'complete two-state fibre proves sharp22 owner-clock floor')
    moves=[(a,b) for a,b,*_ in chosen]
    # Validate on every exceptional fibre, including contexts and reverse moves.
    for k in (4,5):
        for g,C in comps[k].items():
            if len(C)<2:continue
            d=DSU(len(C))
            for a,b in moves:
                for m in bank[k][g]:
                    if m&a==a and m&b==0:
                        n=(m^a)|b;need(n in where[k] and where[k][n][0]==g,'actual feasible same-fibre arrival')
                        d.join(where[k][m][1],where[k][n][1])
                    if m&b==b and m&a==0:
                        n=(m^b)|a;need(n in where[k] and where[k][n][0]==g,'reverse actual arrival')
                        d.join(where[k][m][1],where[k][n][1])
            need(len({d.root(i) for i in range(len(C))})==1,'sparse alphabet connects every small boundary fibre')
    return {'status':'FINITE-EXACT sparse optimal eight-slot completion of original28 binary supply',
            'sources':ids,'counts':counts,'eight_slot_candidate_bridges':all_candidate_count,'candidates_cost_at_most22':len(candidates),'additional_move_types':len(moves),
            'max_chosen_clock_cost':max(work(a|b) for a,b in moves),'sharp_clock_fibre':{'charge':[6,123],'states':clock_witness,'cost':22},
            'moves':moves,'margin5_cut_witnesses':[(g,sorted(L),sorted(R)) for g,L,R in cuts],
            'margin4_bad_components':[(g,C) for g,C in sorted(comps[4].items()) if len(C)>1],
            'input_revision':REV,'inputs':INPUTS}
if __name__=='__main__':
    o=main();(HERE/'bridges.json').write_text(json.dumps(o,indent=2,sort_keys=True)+'\n');print(json.dumps({k:o[k] for k in ('counts','eight_slot_candidate_bridges','additional_move_types')},indent=2))
