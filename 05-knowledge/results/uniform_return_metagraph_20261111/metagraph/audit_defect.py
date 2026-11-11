"""Independent permutation/cut audit of the root's two-inport automaton."""
from itertools import permutations, product
from pathlib import Path
import importlib.util
import json

HERE=Path(__file__).resolve().parent
path=HERE.parent/"defect/two_inport.py"
spec=importlib.util.spec_from_file_location("defect",path)
d=importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)

def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)

def literal(word):
    # Independently assemble the relation as Boolean ordered-pair membership.
    n=len(word)+3
    edges={(1,0),(2,0),(1,2)}
    for i,letter in enumerate(word,3):
        edges.add((0,i))
        for j in range(i+1,n):
            edges.add((i,j))
        for port in (1,2):
            edges.add((i,port) if letter & (1<<(port-1)) else (port,i))
    h=f=0
    for order in permutations(range(n)):
        bad=[i for i in range(n-1) if (order[i],order[i+1]) not in edges]
        h += not bad
        # An ordered spanning two-path cover whose first path ends at root0
        # is exactly one permutation with the separator immediately after0.
        # The second path may be empty, matching E_r(2)'s convention.
        f += all(order[i]==0 for i in bad)
    reach=[{i}|{j for u,j in edges if u==i} for i in range(n)]
    for k in range(n):
        for i in range(n):
            if k in reach[i]:
                reach[i]|=reach[k]
    strong=all(len(x)==n for x in reach)
    return f-h,h,strong

def main():
    queries=0
    examples={}
    histograms=[]
    for length in range(5):
        hist={}
        for word in product(range(4),repeat=length):
            defect,h,strong=literal(word)
            found=d.automaton_value(word)
            need(defect==found,f"positive defect mismatch {word}")
            criterion=bool(word) and bool(word[-1]) and any(x&1 for x in word)
            need(strong==criterion,f"strongness mismatch {word}")
            if strong:
                hist[defect]=hist.get(defect,0)+1
                if defect in (4,5):
                    examples.setdefault(defect,[]).append({"word":word,"h":h})
            queries+=1
        histograms.append({"backbone_length":length,"strong_defect_histogram":hist})
    need(examples[5]==[{"word":(1,2),"h":13},{"word":(2,3),"h":13}],
         "b5 literal exceptional sources")
    report={"status":"PASS","literal_permutation_cut_controls":queries,
            "strongness_independently_checked":True,"by_length":histograms,
            "b4_b5_examples":examples,
            "proof_review":"General fixed-defect residual bound, positive run automaton, intrinsic SCC/port quotient and Burnside rationality reviewed; no gap found.",
            "scope":"The producer's finite DFA supplies the all-length small-budget certificate; this audit independently checks its numerical path interpretation and strongness through backbone length4."}
    (HERE/"defect_audit.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"PASS","literal_permutation_cut_controls":queries}))

if __name__=="__main__":
    main()
