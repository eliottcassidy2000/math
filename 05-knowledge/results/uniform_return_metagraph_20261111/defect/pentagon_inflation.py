"""Exact five-tail responses under arbitrary uniformly attached prefix blocks."""
from pathlib import Path
import json
import two_inport as two
import cycle_residual as cyc


def reverse(rows):
    return tuple(sum(1 << j for j in range(len(rows)) if rows[j] >> i & 1)
                 for i in range(len(rows)))


def response(tail):
    rows=two.rows_for(tail)
    # Remove root; a,b are the first two remaining vertices.
    w=tuple(row >> 1 for row in rows[1:])
    hp,ends=two.path_tables(w)
    _,starts=two.path_tables(reverse(w))
    full=(1 << len(w))-1
    e=sum((ends[m][0]+ends[m][1])*hp[full ^ m] for m in range(1,full+1))
    c=sum(ends[m][0]*starts[full ^ m][1]+ends[m][1]*starts[full ^ m][0]
          for m in range(1,full))
    return dict(tail=tail, W_rows=w, E=e, q=hp[full], C=c, B=hp[full]-c)


def main():
    entries=[]
    tests=0
    for tail in ((),(1,),(2,),(1,2),(2,1)):
        co=response(tail)
        two.require(co['C']==2,'ready pentagon response')
        for a in range(4):
            for b in range(4-a):
                graph=cyc.literal((3,)*a,(3,3,3),(3,)*b+tail)
                defect,h=cyc.literal_defect(graph)
                want=12*co['E']*(1 << (a+b))+3*co['B']
                two.require((h,defect)==(want,6),'triangle inflation')
                tests+=1
        entries.append(co)
    report=dict(status='PASS', scope='Five uniform-prefix tail responses and50 literal triangle-gap instances',
                formula='H=E*Q_P(2)+B*H(P); q=s=q_W*H(P); defect=2*H(P)',
                tests=tests, tails=entries)
    (Path(__file__).with_suffix('.json')).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
