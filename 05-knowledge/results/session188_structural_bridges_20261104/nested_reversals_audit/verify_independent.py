"""Independent literal-word averaging and nested-pair attachment checks."""
from itertools import combinations,permutations,product
from functools import lru_cache
from pathlib import Path
import hashlib,json,random


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def rows(n,code):
    T=[0]*n
    for k,(u,v) in enumerate(combinations(range(n),2)):
        if code>>k&1:T[u]|=1<<v
        else:T[v]|=1<<u
    return T


def flip(T,A):
    U=list(T)
    for u,v in combinations(sorted(A),2):U[u]^=1<<v;U[v]^=1<<u
    return U


def extend(T,a,c):
    n=len(T);U=list(T)+[1<<(n+1),0]
    for u in range(n):
        if a>>u&1:U[u]|=1<<n
        else:U[n]|=1<<u
        if c>>u&1:U[u]|=1<<(n+1)
        else:U[n+1]|=1<<u
    return U


def hp(T,start=None):
    n=len(T);full=(1<<n)-1
    @lru_cache(None)
    def finish(used,last):
        if used==full:return 1
        rem=T[last]&~used;total=0
        while rem:
            b=rem&-rem;rem-=b;v=b.bit_length()-1
            total+=finish(used|b,v)
        return total
    return sum(finish(1<<v,v) for v in (range(n) if start is None else [start]))


def op(T):return [((1<<len(T))-1)^r^(1<<i) for i,r in enumerate(T)]


def is_module(T,A):
    mask=sum(1<<v for v in A)
    return all(T[v]&mask in (0,mask) for v in range(len(T)) if v not in A)


def reverse_runs(word,A):
    out=[];i=0
    while i<len(word):
        if word[i] not in A:out.append(word[i]);i+=1;continue
        j=i
        while j<len(word) and word[j] in A:j+=1
        out.extend(reversed(word[i:j]));i=j
    return tuple(out)


def partial_valid(word,T,A,x):
    for u,v in zip(word,word[1:]):
        if (u==x and v in A) or (v==x and u in A):continue
        if not T[u]>>v&1:return False
    return True


def weight(word,A,x):
    i=word.index(x);left=i>0 and word[i-1] in A;right=i+1<len(word) and word[i+1] in A
    if left:return 1
    return len(A)-1 if right else len(A)


def averaging_words():
    cases=valid_words=0
    for n in (1,2,3,4):
        # Complete through n=3; a declared 10-code slice at n=4.
        codes=range(1<<(n*(n-1)//2)) if n<=3 else range(0,64,7)
        words=list(permutations(range(n+2)))
        for code in codes:
            T=rows(n,code)
            for mask in range(1,1<<n):
                A={i for i in range(n) if mask>>i&1}
                if not is_module(T,A):continue
                V=[i for i in range(n) if i not in A];U=flip(T,A)
                paired=[(w,reverse_runs(w,A)) for w in words]
                for ax,cy,allc in product(range(1<<len(V)),range(1<<len(V)),(False,True)):
                    a=sum(1<<v for j,v in enumerate(V) if ax>>j&1)
                    c=sum(1<<v for j,v in enumerate(V) if cy>>j&1)+(mask if allc else 0)
                    TX,UX=extend(T,a,c),extend(U,a,c);left=right=0
                    for w,z in paired:
                        need(reverse_runs(z,A)==w,'run map not involutive')
                        pv=partial_valid(w,TX,A,n);qv=partial_valid(z,UX,A,n)
                        need(pv==qv,'partially tested word not preserved')
                        need(weight(w,A,n)==weight(z,A,n),'contact weight changed')
                        if pv:left+=weight(w,A,n);valid_words+=1
                        if qv:right+=weight(z,A,n)
                    need(left==right,'weighted query sum differs');cases+=1
    return {'cases':cases,'valid_partial_words':valid_words,
            'universe':'All tournaments through order 3; order-4 codes 0,7,...,63; every nonempty module and all allowed exterior/second-port masks.'}


def substitute(Q,A):
    m=len(A);T=list(A)+[0]*(len(Q)-1)
    for u in range(m):
        for j in range(1,len(Q)):
            if Q[0]>>j&1:T[u]|=1<<(m+j-1)
            else:T[m+j-1]|=1<<u
    for i,j in combinations(range(1,len(Q)),2):
        a,b=m+i-1,m+j-1
        if Q[i]>>j&1:T[a]|=1<<b
        else:T[b]|=1<<a
    return T


def nested_batteries():
    rng=random.Random(188812);records=[]
    for n in range(3,9):
        for trial in range(4):
            m=2+trial%(n-2);Q=rows(n-m+1,rng.randrange(1<<((n-m+1)*(n-m)//2)))
            A=rows(m,rng.randrange(1<<(m*(m-1)//2)));B=substitute(Q,A)
            v=m+trial%(n-m);b=(hp(B,v)+hp(op(B),v))%2;d=B[v].bit_count()
            # Keep the old B names first; add up to one literal outside vertex.
            T=list(B)
            if trial%2:
                toward=bool(rng.randrange(2));T.append((1<<n)-1 if toward else 0)
                if not toward:
                    for u in range(n):T[u]|=1<<n
            outer=set(range(n));inner=set(range(m));only=flip(T,outer);both=flip(only,inner)
            ext=sum(1<<u for u in range(n,len(T)) if T[u]&1)
            queries=[]
            for u in range(n):
                if u!=v:queries.append((ext|(1<<u)|((1<<v) if B[v]>>u&1 else 0),ext|(1<<v)))
            if d%2:queries.append((ext|(1<<v),ext|(1<<v)))
            if (n-1+d)%2:queries.append((ext,ext|(1<<v)))
            need(len(queries)==len(set(queries)),'duplicate battery')
            counts=[(a,c,hp(extend(T,a,c)),hp(extend(only,a,c)),hp(extend(both,a,c))) for a,c in queries]
            need(sum(x[3] for x in counts)==sum(x[4] for x in counts),'exact nested averaging failed')
            need(sum(x[2]-x[4] for x in counts)%4==2*(1+b)%4,'odd/even battery residue')
            endpoint=(ext|(1<<v),ext)
            ends=[hp(extend(S,*endpoint)) for S in (T,only,both)]
            need(ends[1]==ends[2],'endpoint query not invariant')
            if b:need((ends[0]-ends[2])%4==2,'endpoint separator failed')
            else:need(any((x[2]-x[4])%4==2 for x in counts),'battery has no literal separator')
            records.append({'core_order':n,'inner_size':m,'T':T,'single_outer':only,'simultaneous':both,
                            'anchor':v,'endpoint_parity':b,'battery':counts,'endpoint_query':endpoint,'endpoint_counts':ends})
    return records


def main():
    source=Path(__file__).resolve().parent.parent/'simultaneous_response'/'ANCHORED.md'
    print(json.dumps({'status':'PASS','averaging_word_bijection':averaging_words(),
                     'nested_literal_batteries':nested_batteries(),
                     'audited_proof_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                     'method':'No producer imports. Literal partially tested Hamilton words and actual full attachment counts retain all degrees.'},indent=2,sort_keys=True))


if __name__=='__main__':main()
