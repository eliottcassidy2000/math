"""Check the integer separating hyperplane for seed18, using certified normalized words on stdin."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import sys

# Entry j applies to a cell used by between 4*j and 4*j+3 normalized words.
WEIGHTS = (-17,-17,-17,-15,-17,-17,-3,-3,-2,-2,-3,-1,0,3,2,3,3,7,4,8,
           9,10,2,11,17,6,6,5,16,3,10,1,2,-7,10,17,10,17,4)


def evaluate(words):
    words=sorted(words)
    assert len(words)==len(set(words))==2479
    pairs=list(combinations(range(18),2))
    assert all(len(w)==18 and w[0]==0 and all(0<=v<56 for v in w) for w in words)
    differences=[tuple((w[j]-w[i])%56 for i,j in pairs) for w in words]
    cells=[Counter(col) for col in zip(*differences)]
    assert all(len(c)==39 for c in cells)
    occurrences=[v for c in cells for v in c.values()]
    assert (min(occurrences),max(occurrences))==(14,152)
    tables=[{d:WEIGHTS[n//4] for d,n in c.items()} for c in cells]
    scores=[sum(t[d] for t,d in zip(tables,row)) for row in differences]
    total=sum(sum(t.values()) for t in tables)
    assert (len(pairs),min(scores),max(scores),total)==(153,503,956,19598)
    assert 39*min(scores)-total==19
    print("All 2479 scores are at least",min(scores),"while the balanced mean is",str(total)+"/39.")
    print("The prescribed total is",56*total,"and all legal allocations total at least",56*39*min(scores))
    file=Path(__file__).with_name("seed18_balance_weights.txt")
    data=[list(map(int,s.split())) for s in file.read_text().splitlines() if not s.startswith("#")]
    upper=[]
    for (i,j),freq,row in zip(pairs,cells,data):
        assert len(row)==41 and row[:2]==[i,j] and min(row[2:])>=0
        upper.append(dict(zip(sorted(freq),row[2:])))
    total=sum(sum(t.values()) for t in upper)
    cost=min(sum(t[d] for t,d in zip(upper,word)) for word in differences)
    assert (total,cost)==(2865315,99994)
    assert 56*total//cost==1604
    allocations=dict(data[len(pairs):])
    assert len(allocations)==len(data)-len(pairs)==973
    loads=[Counter() for _ in pairs]
    for k,weight in allocations.items():
        assert weight>0 and 0<=k<len(words)
        for load,d in zip(loads,differences[k]):
            load[d]+=weight
    maxload=max(max(c.values()) for c in loads)
    primal=sum(allocations.values())
    assert (maxload,primal)==(9999991,286532033)
    assert 56*primal>1604*maxload and 56*total<1605*cost
    print("Fractional capacity lies between 56*%d/%d and 56*%d/%d." %
          (primal,maxload,total,cost))


if __name__=="__main__":
    evaluate([tuple(map(int,line.split())) for line in sys.stdin if line.strip()])
