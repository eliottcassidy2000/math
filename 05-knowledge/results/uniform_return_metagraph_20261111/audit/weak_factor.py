"""Independent all-labelled small control of rooted SCC defect factorization."""
from pathlib import Path
import json
import retained_reset as audit
from template_controls import graph,residual_sccs


def subsets(mask):
    sub=mask
    while True:
        yield sub
        if sub==0:return
        sub=(sub-1)&mask


def main():
    counts=dict(graphs=0, rooted_factorizations=0, branching_cases=0, leading_prefix_bound_cases=0)
    for n in range(1,6):
        for code in range(1 << (n*(n-1)//2)):
            rows=graph(n,code)
            hp,end=audit.tables(rows)
            full=(1 << n)-1
            components=residual_sccs(rows,full)
            counts['graphs']+=1
            for r in range(n):
                c=next(x for x in components if x>>r&1)
                first=(c & -c).bit_length()-1
                a=sum(x for x in components if x!=c and rows[(x & -x).bit_length()-1]>>first&1)
                b=full ^ a ^ c
                f_c=sum(end[m][r]*hp[c ^ m] for m in subsets(c))
                f=sum(end[m][r]*hp[full ^ m] for m in range(1,full+1))
                q_a=sum(hp[m]*hp[a ^ m] for m in subsets(a))
                defect=f-hp[full]
                d=sum(bool(row>>r&1) for row in rows)
                audit.need(f==q_a*f_c*hp[b], 'rooted ordered-cover product')
                audit.need(hp[full]==hp[a]*hp[c]*hp[b], 'ordinary Hamilton product')
                audit.need(defect==hp[b]*(q_a*f_c-hp[a]*hp[c]), 'defect factor')
                counts['rooted_factorizations']+=1
                if d>=2:
                    audit.need(defect>=2**d-2, 'weak exponential indegree bound')
                    counts['branching_cases']+=1
                    if a:
                        audit.need(hp[full]<=defect, 'leading prefix bounds total H')
                        counts['leading_prefix_bound_cases']+=1
    report=dict(status='PASS', universe='All labelled tournaments of orders1..5 and every root',
                method='Independent literal SCC reachability and endpoint-subset tables', **counts)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(report)


if __name__=='__main__':main()
