"""Literal X_n reset-and-complement quotient, independently refined through seven."""
from audit import canon,pres,partner,need,relabel
from collections import Counter,defaultdict
import json
def sigma(t):return relabel(partner(t,()),tuple(reversed(range(len(t)))))
def join(a,b):return tuple(x|(((1<<len(b))-1)<<len(a)) for x in a)+tuple(x<<len(a) for x in b)
def closure(obs,groups):
    n=len(obs);cells=[];ids={}
    for o in obs:cells.append(ids.setdefault(o,len(ids)))
    fibres=defaultdict(list)
    for i,g in enumerate(groups):fibres[g].append(i)
    rounds=0
    while True:
        rows={g:tuple(sorted(Counter(cells[i] for i in f).items())) for g,f in fibres.items()}
        ids={};new=[]
        for i in range(n):
            sig=(cells[i],cells[i^(n-1)],rows[groups[i]])
            new.append(ids.setdefault(sig,len(ids)))
        if len(ids)==len(set(cells)):return cells,rounds
        cells=new;rounds+=1
def main():
    out=[]
    for n in range(1,8):
        ts=[pres(n,m) for m in range(1<<((n-1)*(n-2)//2))]
        f=[canon(t) for t in ts];g=[canon(t,(0,n-1)) for t in ts];b=[sigma(t)==t for t in ts];N=len(ts)
        row={'n':n,'presentations':N,'protocols':{}}
        for rooted in (False,True):
            for colour in (False,True):
                q=g if rooted else f;obs=[(q[i],b[i]) if colour else q[i] for i in range(N)]
                for fixed in (False,True):
                    if fixed and not rooted:continue
                    cells,steps=closure(obs,g if fixed else f)
                    expect=[(obs[i],obs[i^(N-1)]) for i in range(N)]
                    ec=defaultdict(set);ce=defaultdict(set)
                    for e,c in zip(expect,cells):ec[e].add(c);ce[c].add(e)
                    need(all(len(v)==1 for v in ec.values()) and all(len(v)==1 for v in ce.values()),'reset-pair theorem exact partition')
                    key=('joint' if rooted else 'free')+('/colour' if colour else '')+('/fixed_reset' if fixed else '/free_reset')
                    row['protocols'][key]={'cells':len(set(cells)),'strict_refinements':steps}
        out.append(row)
    x=(2,12,41,17,39,75,31);y=(2,20,9,51,37,71,31)
    def q(t):return (canon(t,(0,6)),canon(partner(t,tuple(range(7))),(0,6)),sigma(t)==t)
    need(x!=y and q(x)==q(y),'seven-vertex reset compression witness')
    a=join(x,sigma(x));b=join(y,sigma(x))
    need(sigma(a)==a and sigma(b)!=b,'blue join distinguishes compressed twins')
    print(json.dumps({'literal_reset_quotients':out,'blue_join_hostile':{'left_rows':[x,y],'right_rows':sigma(x),'blue_outputs':[True,False],'same_input_joint_pairs_and_colour':True}},sort_keys=True,indent=2))
if __name__=='__main__':main()
