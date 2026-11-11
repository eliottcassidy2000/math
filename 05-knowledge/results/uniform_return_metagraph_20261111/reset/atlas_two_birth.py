"""Complete homogeneous-prefix two-birth dictionary, plus actual C4 chart composition."""
from pathlib import Path
import json
from classify_modules import H,strong,need,inherited as a,two_birth as b
from atlas_returns import ATLAS,kernel,cycle
HERE=Path(__file__).resolve().parent


def construct(tail,drop,cut1,cut2,k):
    source=a.ready_source(tail,k);kept=[v for v in range(len(source)) if v!=k+drop]
    parent=a.induced(source,kept)
    def lift(out):return ((out>>1)<<k)|(((1<<k)-1) if out&1 else 0)
    first=a.attach(parent,lift(cut1));last=a.attach(first,lift(cut2))
    return source,kept,parent,first,last


def main():
    report=[];total=candidates=passes=0
    for word,A,B in ATLAS:
        tail=tuple({'00':0,'01':2,'10':1,'11':3}[s] for s in word.split())
        source=a.ready_source(tail,1);need(H(source)==2*A+B,'atlas input')
        row={'tail':word,'A':A,'B':B,'base_scalar_candidates':0,'failure_channels':{},'operators':[]}
        for drop in range(len(tail)):
            kept=[v for v in range(len(source)) if v!=drop+1];P=a.induced(source,kept)
            for cut1 in range(1<<len(P)):
                first=a.attach(P,cut1);vals=b.all_birth_counts(first);total+=len(vals)
                for cut2,h in enumerate(vals):
                    offset=h-3*(2*A+B)
                    if offset not in(2,4):continue
                    candidates+=1;row['base_scalar_candidates']+=1
                    T=a.attach(first,cut2);need(H(T)==h,'independent literal base candidate')
                    co=b.module_runs(T)
                    bad=next((j for j in range(3,len(co)) if co[j]),None)
                    if co[2]!=3*A or bad is not None:
                        failure=('A2' if co[2]!=3*A else 'A'+str(bad))
                        row['failure_channels'][failure]=row['failure_channels'].get(failure,0)+1
                        continue
                    passes+=1;need(strong(T),'all successful targets are strong')
                    root=len(P)-1;K=kernel(T,root)
                    need([v for v in range(len(T)) if T[v]>>root&1]==[root-2,root-1],
                         'exact retained old inports')
                    need(K['E']==3*A and K['q']-K['C']==3*B+offset,'independent all-k kernel law')
                    controls=[]
                    for k in range(1,6):
                        S,ids,W,U,V=construct(tail,drop,cut1,cut2,k)
                        need(H(S)==A*(1<<k)+B and H(V)==3*H(S)+offset,'literal all-k branch')
                        need(a.induced(V,range(len(W)))==a.induced(S,ids),'all retained old arcs')
                        need(V==b.inflate(T,0,k) and strong(V),'fixed prefix template')
                        p=a.root_profiles(V)[len(W)-1]
                        need((p['s'],p['q'],p['gap'],p['delta'])==(K['q'],K['q'],K['C'],K['C']),
                             'full old-root interface')
                        controls.append({'k':k,'kept_old_labels':ids,'trace':[H(S),H(W),H(U),H(V)],'profile':p})
                    row['operators'].append({'drop':drop,'cut1':cut1,'cut2':cut2,'offset':offset,
                        'base_rows':T,'module_runs':co,'kernel':K,'spanning_cycle':cycle(T),'controls':controls})
        report.append(row)
    need((total,candidates,passes)==(270336,2750,30),'complete universe totals')
    composed=[]
    for k in range(1,7):
        S,ids,W,U,V=construct((0,1),0,4,8,k)
        atlas=a.ready_source((1,1,2),k)
        # Output rows have old ports/root before the two newborn tail vertices.
        mapping=list(range(k+1))+[k+3,k+4,k+5,k+1,k+2]
        need(a.induced(atlas,mapping)==V,'actual internal C4+4 graph map')
        # Apply next instruction in its named atlas coordinates.
        S2,ids2,W2,U2,V2=construct((1,1,2),0,4,14,k)
        need(S2==atlas and H(V2)==9*H(S)+14,'common actual intermediate and two-step count')
        composed.append({'k':k,'first_target_mapping_to_atlas':mapping,
                         'counts':[H(S),H(V),H(V2)],'second_retained_atlas_labels':ids2})
    data={'status':'PASS','scope':'All16 neutral4 tails, exactly one tail deletion, two arbitrary homogeneous-prefix births',
          'ordered_cut_pairs':total,'scalar_candidates':candidates,'uniform_instructions':passes,
          'atlas':report,'internal_C4_plus4_then_plus2':composed}
    (HERE/'atlas_two_birth.json').write_text(json.dumps(data,sort_keys=True,indent=2)+'\n')
    print('PASS:270336 cut pairs,2750 scalar candidates,30 all-k instructions; actual C4+4 then+2 composition')


if __name__=='__main__':main()
