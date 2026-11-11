"""Literal controls for the all-d minimum branching-defect classification."""
from itertools import combinations, product
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)

def transitive(n):
    return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def build(inports,word):
    d=len(inports); n=d+len(word)+1
    rows=[sum(1<<v for v in range(d+1,n))]
    rows += [(row<<1)|1 for row in inports]
    rows += [0]*len(word)
    for i,letter in enumerate(word):
        v=d+1+i
        rows[v] |= sum(1<<w for w in range(v+1,n))
        for port in range(d):
            if letter>>port&1:
                rows[v] |= 1<<(port+1)
            else:
                rows[port+1] |= 1<<v
    return tuple(rows)

def path_data(rows):
    n=len(rows); full=(1<<n)-1
    end=[[0]*n for _ in range(full+1)]
    for v in range(n):end[1<<v][v]=1
    for mask in range(1,full+1):
        for v,count in enumerate(end[mask]):
            if count:
                rest=rows[v] & (full^mask)
                while rest:
                    bit=rest & -rest;rest-=bit
                    end[mask|bit][bit.bit_length()-1]+=count
    h=[1]+[sum(row) for row in end[1:]]
    return h,end

def profile(rows):
    h,end=path_data(rows)
    full=len(h)-1
    f=sum(end[mask][0]*h[full^mask] for mask in range(1,full+1) if mask&1)
    return h[full],f-h[full]

def strong(rows):
    reach=[row|(1<<v) for v,row in enumerate(rows)]
    for k in range(len(rows)):
        for v in range(len(rows)):
            if reach[v]>>k&1:reach[v]|=reach[k]
    return all(row==(1<<len(rows))-1 for row in reach)

def all_labelled(n):
    pairs=list(combinations(range(n),2))
    for code in range(1<<len(pairs)):
        rows=[0]*n
        for j,(a,b) in enumerate(pairs):
            if code>>j&1:rows[a]|=1<<b
            else:rows[b]|=1<<a
        yield tuple(rows)

def main():
    reports=[]
    for name,I,max_length in (("TT3",transitive(3),4),("C3",(2,4,1),4),
                              ("TT4",transitive(4),3),("TT5",transitive(5),2)):
        d=len(I);budget=(1<<d)-2;queries=strong_count=equalities=0;by_length=[]
        for length in range(max_length+1):
            equal=[]
            for word in product(range(1<<d),repeat=length):
                rows=build(I,word);h,b=profile(rows);sc=strong(rows)
                need(b>=budget,f"minimum bound {name} {word}")
                if sc:
                    strong_count+=1
                    expected=bool(word) and all(letter==(1<<d)-1 for letter in word)
                    need((b==budget)==expected,f"minimum classification {name} {word}")
                    if b==budget:equal.append({"word":word,"H":h});equalities+=1
                queries+=1
            by_length.append({"length":length,"equality_sources":equal})
        reports.append({"inports":name,"d":d,"budget":budget,
                        "queries":queries,"strong_sources":strong_count,
                        "equalities":equalities,"by_length":by_length})
    # All64 labelled I of order4 and all16 one-backbone contact patterns.
    labelled_queries=0
    for I in all_labelled(4):
        Ih,_=path_data(I);is_transitive=Ih[-1]==1
        for letter in range(16):
            rows=build(I,(letter,));h,b=profile(rows)
            need(b>=14,"all-inport lower bound")
            if strong(rows):
                need((b==14)==(is_transitive and letter==15),"all-inport equality")
            labelled_queries+=1
    collar_controls=[]
    for d in range(2,7):
        for I in (transitive(d),)+(((2,4,1),) if d==3 else ()):
            hs,_=path_data(I);HI=hs[-1];full=len(hs)-1
            Q=sum(hs[a]*hs[full^a] for a in range(full+1))
            for k in range(1,5):
                rows=build(I,((1<<d)-1,)*k);h,b=profile(rows)
                need(h==(Q-HI)*(1<<k)+3*HI-Q,"collar H formula")
                need(b==Q-2*HI,"collar defect formula")
                collar_controls.append({"d":d,"I_H":HI,"k":k,"H":h,"defect":b})
    # The symmetric C3 collar gives an actual source of H141.
    source=build((2,4,1),(7,7,7,7));n=len(source)
    need(profile(source)==(141,6),"141 source")
    targets={425:[],427:[]}
    for mask in range(1<<n):
        child=tuple(row|(0 if mask>>v&1 else 1<<n) for v,row in enumerate(source))+(mask,)
        h,_=path_data(child)
        if h[-1] in targets:targets[h[-1]].append(mask)
    need(not any(targets.values()),"unexpected ternary repair")
    report={"status":"PASS","word_experiments":reports,
            "labelled_order4_inport_queries":labelled_queries,
            "collar_formula_controls":collar_controls,
            "H141_source":source,"H141_one_birth_targets":targets,
            "scope":"Finite literal controls; all-order equality follows from positive tagged five-cycle proof in MINIMAL_DEFECT.md."}
    (HERE/"minimal_defect.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","word_queries":sum(x["queries"] for x in reports),
                      "labelled_inport_queries":labelled_queries,
                      "collar_controls":len(collar_controls),"H141_birth_queries":1<<n}))

if __name__=="__main__":main()
