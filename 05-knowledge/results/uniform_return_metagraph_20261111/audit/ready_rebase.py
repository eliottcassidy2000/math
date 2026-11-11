"""Fresh named-graph check of the eight-retained 425 return and future order."""
from pathlib import Path
import json
import retained_reset as audit


def frame(names,edges,r):
    rows=audit.masks(names,edges)
    hp,ends=audit.tables(rows)
    rev=tuple(sum(1 << j for j in range(len(rows)) if rows[j]>>i&1) for i in range(len(rows)))
    _,starts=audit.tables(rev)
    full=(1 << len(rows))-1; root=names.index(r); without=full ^ (1 << root)
    h=hp[full]; e=ends[full][root]; s=starts[full][root]
    u=sum(starts[m][root]*hp[full ^ m] for m in range(1,full+1))
    v=sum(ends[m][root]*hp[full ^ m] for m in range(1,full+1))
    d=hp[without]
    k=sum(hp[m]*hp[without ^ m] for m in range(full+1) if m & ~without==0)
    return (h,e,s,u,v,d,k)


def birth(names,edges,name,only_in):
    result=set(edges)
    result.update((old,name) if old==only_in else (name,old) for old in names)
    return names+(name,),result


def main():
    traces=[]
    for k in range(1,7):
        original,original_edges=audit.source(k)
        names=tuple(x for x in original if x!='M')
        edges={(a,b) for a,b in original_edges if a!='M' and b!='M'}
        names,edges=audit.gadget(names,edges,1)
        before=frame(names,edges,'r')
        keep=tuple(x for x in names if x not in {'b','1x'})
        inherited={(a,b) for a,b in edges if a in keep and b in keep}
        parent_rows=audit.masks(keep,inherited)
        hp,_=audit.tables(parent_rows)
        audit.need(hp[-1]==11*(1<<k)+9,'rebase parent formula')
        after_names,after_edges=birth(keep,inherited,'v','1z')
        after=frame(after_names,after_edges,'v')
        h=24*3**k+53*2**k-3
        audit.need(after[0]==h and after[4]-h==1 and after[5]==11*2**k+9,'ready response')
        audit.need(after[2]==after[5]-1,'unique-start defect')
        audit.need(audit.strong(audit.masks(after_names,after_edges)),'strong returned source')
        audit.need({(a,b) for a,b in after_edges if a in keep and b in keep}==inherited,'all eight retained arcs')
        one_names,one_edges=birth(after_names,after_edges,'u','v')
        one=frame(one_names,one_edges,'u')
        two_names,two_edges=birth(one_names,one_edges,'w','u')
        two=frame(two_names,two_edges,'w')
        audit.need(two[0]==3*h+2,'actual next ternary branch')
        traces.append(dict(k=k,parent_H=hp[-1],before_frame=before,ready_frame=after,
                           next_counts=[one[0],two[0]],retained_names=keep))
    old=traces[1]['before_frame']; new=traces[1]['ready_frame']
    audit.need(old==(425,4,41,2542,449,41,2848),'old seven frame')
    audit.need(new==(425,15,52,1491,426,53,1544),'returned seven frame')
    # Imported cone's exact inequalities; its all-context theorem is cited separately.
    audit.need(all(old[j]>=new[j] for j in (0,3,4,6)) and
               all(old[j]<=new[j] for j in (1,2,5)),'two-sided amalgamation order cone')
    report=dict(status='PASS', scope='Six literal source scales; eight-retained425 return; exactseven-coordinate comparison',
                no_construction_lane_imports=True,traces=traces)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('PASS: six ready returns, actual continuations, and the two-sided rooted response cone')


if __name__=='__main__':main()
