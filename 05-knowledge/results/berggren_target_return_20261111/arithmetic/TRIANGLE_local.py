"""Independent local cover counts, including the repeated-contact hostile."""
from itertools import product,permutations
from pathlib import Path
import json

def count(q,start=None,end=None):
    total=0
    for assignment in product(range(q),repeat=3):
        if start is not None and assignment[start]!=0:continue
        if end is not None and assignment[end]!=0:continue
        weight=1
        for bin_index in range(q):
            vertices=[v for v in range(3) if assignment[v]==bin_index]
            ways=0
            for order in permutations(vertices):
                if bin_index==0 and start is not None and (not order or order[0]!=start):continue
                if bin_index==0 and end is not None and (not order or order[-1]!=end):continue
                if all((a+1)%3==b for a,b in zip(order,order[1:])):ways+=1
            weight*=ways
        total+=weight
    return total

def need(ok,why):
    if not ok:raise RuntimeError(why)

def main():
    rows=[]
    for q in range(1,13):
        free=count(q);start=count(q,start=0);end=count(q,end=0)
        need(free==q**3+2*q,('free',q,free))
        need(start==end==q*q-q+1,('one event',q,start,end))
        need((start%3==0)==(q%3==2),('local divisibility',q))
        if q%3==2:need(start%9!=0,('exact valuation',q))
        rows.append(dict(q=q,unmarked=free,fixed_start=start,fixed_end=end))
    hostile=count(2,start=0,end=1)
    need(hostile==1,'repeated contact hostile')
    source=[58,60,57,16,32,8]
    extensions={}
    for mask in range(64):
        graph=[row|((1<<6) if not(mask>>v&1) else 0) for v,row in enumerate(source)]+[mask]
        h=sum(all(graph[a]>>b&1 for a,b in zip(path,path[1:])) for path in permutations(range(7)))
        extensions[mask]=h
    spectrum=sorted(set(extensions.values()))
    need(spectrum==[9,15,31,41,45,61,75,99],'complete adjacent-triangle one-birth spectrum')
    need(extensions[25]==31 and 29 not in spectrum,'content versus prescribed target')
    report=dict(status='FINITE-EXACT independent assignment/permutation audit of local formulas, not an all-k proof',
                q_values=rows,repeated_contact=dict(q=2,start=0,end=1,count=hostile,
                                                  sole_cover=[[0,1],[2]]),
                adjacent_triangles_one_birth=dict(source=source,all_cut_counts=extensions,spectrum=spectrum,
                                                  plus4_mask=25,plus2_target_29_absent=True))
    Path(__file__).with_name('TRIANGLE_local.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(result='PASS',q_range=[1,12],hostile=hostile)))

if __name__=='__main__':main()
