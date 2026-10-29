"""Tree-code conjugacy and the exact span of unrestricted tree complements."""
from pathlib import Path
from itertools import product,combinations
import json
import probes as p
HERE=Path(__file__).resolve().parent
tree=p.r.old.source('selector_tree_codec','05-knowledge/results/session90_even_graphs_odd_descent_20261029/trijection/check_tree_compression.py')


def insert(basis,x):
    while x:
        i=x.bit_length()-1
        if i not in basis:basis[i]=x;return
        x^=basis[i]


def main():
    out=[]
    for n in range(2,8):
        edges=list(combinations(range(n),2));index={e:i for i,e in enumerate(edges)}
        full=(1<<len(edges))-1;basis={};difference={};first=None;count=0;codec_checks=0
        for word in product(range(n),repeat=n-2):
            F=tree.tree_from_prufer(n,word);f=sum(1<<index[e] for e in F);move=full^f
            insert(basis,move)
            if first is None:first=move
            insert(difference,move^first);count+=1
            if n<=4:
                smaller_full=(1<<((n-1)*(n-2)//2))-1
                for x in range(1<<len(edges)):
                    p.need(tree.compress(n,F,x^move)==(tree.compress(n,F,x)^smaller_full),'tree complement conjugates to converse')
                    codec_checks+=1
        expected=len(edges) if ((n-1)*(n-2)//2)%2 else len(edges)-1
        p.need(len(basis)==expected and len(difference)==len(edges)-1,'full/even translation rank')
        out.append(dict(order=n,trees=count,pair_bits=len(edges),complement_rank=len(basis),
                        difference_rank=len(difference),labelled_components=1<<(len(edges)-len(basis)),codec_checks=codec_checks))
    result=dict(status='FINITE-EXACT',scope='All underlying spanning trees allowed; current tree orientation is induced, no Hamilton legality restriction.',orders=out)
    (HERE/'tree-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result));return result


if __name__=='__main__':main()
