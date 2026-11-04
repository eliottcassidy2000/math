"""Literal squarefree-path controls and hostile tests for parity shortcuts."""
from argparse import ArgumentParser
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path
import gzip
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE=Path(__file__).read_bytes()
CHECKS=0


def require(ok, why):
    global CHECKS
    CHECKS+=1
    if not ok: raise ValueError(why)


def rows(n, code):
    a=[0]*n
    for bit,(u,v) in enumerate(combinations(range(n),2)):
        if code>>bit&1:a[u]|=1<<v
        else:a[v]|=1<<u
    return tuple(a)


def masks(poly):
    while poly:
        bit=poly&-poly
        yield bit.bit_length()-1
        poly-=bit


def multiply(a,b):
    answer=0
    for i in masks(a):
        for j in masks(b):
            if not i&j: answer^=1<<(i|j)
    return answer


def xor(items):
    answer=0
    for item in items:answer^=item
    return answer


def literal(a):
    n=len(a);full=(1<<n)-1
    G=[[0]*n for _ in range(n)]
    S=[0]*n;E=[0]*n;C=[[0]*n for _ in range(n)]
    for k in range(1,n+1):
        for word in permutations(range(n),k):
            if all(a[u]>>v&1 for u,v in zip(word,word[1:])):
                mask=sum(1<<u for u in word)
                G[word[0]][word[-1]]^=1<<mask
                if k==n:S[word[0]]+=1;E[word[-1]]+=1
    for word in permutations(range(n)):
        bad=[i for i,(u,v) in enumerate(zip(word,word[1:])) if not a[u]>>v&1]
        if not bad:
            for u,v in zip(word,word[1:]):C[u][v]+=1
        elif len(bad)==1:
            i=bad[0];C[word[i]][word[i+1]]+=1
    s=[xor(row) for row in G]
    e=[xor(G[i][j] for i in range(n)) for j in range(n)]
    U=(1<<(1<<n))-1
    require(xor(s)==U^1,'Redei on every subset')
    for u in range(n):
        require(e[u]==multiply(U,s[u]),'whole endpoint zeta transform')
        require(s[u]==multiply(U,e[u]),'inverse endpoint zeta transform')
        for v in range(n):
            require(G[u][v]^G[v][u]==multiply(s[u],s[v]),'whole squarefree matrix identity')
            require(((multiply(e[u],s[v])>>full)&1)==C[u][v]%2,'direct two-cover parity')
            require((G[u][v]^G[v][u]).bit_count()%2==C[u][v]%2,'all-simple-path interpretation')
    return {'rows':a,'S':S,'E':E,'H':sum(S),'C':C}


def main():
    parser=ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    counts=[]
    for n in range(1,6):
        count=1<<(n*(n-1)//2)
        for code in range(count):literal(rows(n,code))
        counts.append([n,count])
    bank=ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz'
    raw=bank.read_bytes();members=json.loads(gzip.decompress(raw))['members']
    for item in members['6']:literal(rows(6,item['code']))
    scout=json.loads((HERE/'endpoints.json').read_text())
    hostile_codes=scout['orders'][7]['even_parity_balanced_codes']
    require(len(hostile_codes)==2,'complete order8 parity hostiles')
    hostile=[]
    for code in hostile_codes:
        record=literal(rows(8,code));record['code']=code
        require(record['S']!=record['E'] and all((x-y)%2==0 for x,y in zip(record['S'],record['E'])),
                'parity balance does not imply exact balance')
        hostile.append(record)
    dot=literal(rows(6,689));dot['code']=689
    require(sum(x*y for x,y in zip(dot['S'],dot['E']))%2==1,'even-order dot shortcut is false')
    require(Path(__file__).read_bytes()==SOURCE,'source stable')
    result={'status':'FINITE-EXACT PASS','source_sha256':sha256(SOURCE).hexdigest(),
            'checks':CHECKS,'all_labelled_orders':counts,'all_class_controls_order6':len(members['6']),
            'bank_sha256':sha256(raw).hexdigest(),'even_parity_hostiles':hostile,'dot_parity_hostile':dot,
            'scope':'Every coefficient of every path-matrix identity through labelled order5; all56 classes at6; two exact order8 hostiles. Universal proofs are separate.'}
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text,encoding='utf8',newline='\n')
    else:print(text,end='')


if __name__=='__main__':main()
