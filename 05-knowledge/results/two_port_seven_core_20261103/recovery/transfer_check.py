"""Two-port transfer on arbitrary digraphs and undirected square-sum graphs."""
from argparse import ArgumentParser
from hashlib import sha256
from itertools import combinations, permutations, product
from math import isqrt
from pathlib import Path
import json

SOURCE=Path(__file__).read_bytes()


def need(condition,message):
    if not condition: raise RuntimeError(message)


def polynomial_paths(rows):
    n=len(rows);size=1<<n
    g=[[[0]*size for _ in rows] for _ in rows]
    for length in range(1,n+1):
        for word in permutations(range(n),length):
            if all(rows[u]>>v&1 for u,v in zip(word,word[1:])):
                g[word[0]][word[-1]][sum(1<<u for u in word)]+=1
    return g


def multiply(a,b):
    out=[0]*len(a)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y and not i&j:out[i|j]+=x*y
    return out


def formula(g,n,a,b,c,d,alpha,beta):
    size=1<<n
    def end(mask):
        return [sum(g[i][j][s] for i in range(n) for j in range(n) if mask>>j&1) for s in range(size)]
    def start(mask):
        return [sum(g[i][j][s] for i in range(n) for j in range(n) if mask>>i&1) for s in range(size)]
    def middle(left,right):
        return [sum(g[i][j][s] for i in range(n) for j in range(n) if left>>i&1 and right>>j&1) for s in range(size)]
    def term(a,b,c,d,arc):
        left,center,right=end(a),middle(b,c),start(d)
        left[0]+=1;center[0]+=arc;right[0]+=1
        return multiply(multiply(left,center),right)[-1]
    return term(a,b,c,d,alpha)+term(c,d,a,b,beta)


def literal(rows,a,b,c,d,alpha,beta):
    n=len(rows);out=rows[:]+[b,d]
    for u in range(n):
        if a>>u&1:out[u]|=1<<n
        if c>>u&1:out[u]|=1<<(n+1)
    if alpha:out[n]|=1<<(n+1)
    if beta:out[n+1]|=1<<n
    return sum(all(out[u]>>v&1 for u,v in zip(word,word[1:]))
               for word in permutations(range(n+2)))


def graph(n,mask,directed):
    rows=[0]*n
    pairs=[(u,v) for u in range(n) for v in range(n) if u!=v] if directed else list(combinations(range(n),2))
    for bit,(u,v) in enumerate(pairs):
        if mask>>bit&1:
            rows[u]|=1<<v
            if not directed:rows[v]|=1<<u
    return rows


def main():
    parser=ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    directed=undirected=0
    for n in range(3):
        for code in range(1<<(n*(n-1))):
            rows=graph(n,code,True);g=polynomial_paths(rows)
            for a,b,c,d in product(range(1<<n),repeat=4):
                for alpha,beta in product(range(2),repeat=2):
                    need(formula(g,n,a,b,c,d,alpha,beta)==literal(rows,a,b,c,d,alpha,beta),'arbitrary-digraph response')
                    directed+=1
    for n in range(4):
        for code in range(1<<(n*(n-1)//2)):
            rows=graph(n,code,False);g=polynomial_paths(rows)
            for a,c in product(range(1<<n),repeat=2):
                for arc in range(2):
                    value=formula(g,n,a,a,c,c,arc,arc)
                    need(value==literal(rows,a,a,c,c,arc,arc) and value%2==0,'undirected response and reversal pairing')
                    undirected+=1
    # The square-sum path1--3--6 closes to the five-cycle after adding10,15.
    labels=[1,3,6]
    square=lambda value:isqrt(value)**2==value
    rows=[sum(1<<v for v in range(3) if u!=v and square(labels[u]+labels[v])) for u in range(3)]
    need(rows==[2,5,2] and sum(1<<u for u in range(3) if square(labels[u]+10))==4
         and sum(1<<u for u in range(3) if square(labels[u]+15))==1 and square(10+15),'square-sum incidences')
    g=polynomial_paths(rows)
    value=formula(g,3,4,4,1,1,1,1)
    need(value==literal(rows,4,4,1,1,1,1)==10,'actual square-sum C5 construction')
    need(Path(__file__).read_bytes()==SOURCE,'source stable')
    result={'status':'FINITE-EXACT PASS','source_sha256':sha256(SOURCE).hexdigest(),
            'all_digraph_controls_orders0_to2':directed,'all_undirected_controls_orders0_to3':undirected,
            'square_sum_example':{'old_labels':[1,3,6],'new_labels':[10,15],
                                  'directed_hamilton_words':value,'unoriented_hamilton_paths':value//2}}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text,encoding='utf8',newline='\n')
    else:print(text,end='')


if __name__=='__main__':main()
