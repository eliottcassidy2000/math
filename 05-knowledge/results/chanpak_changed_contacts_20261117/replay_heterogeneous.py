from itertools import permutations, combinations, combinations

def arc(rows,i,j):return rows[i]>>j&1

def skeleton(rows,r):
    parts=[i for i in range(len(rows)) if i!=PORT]
    answer=0
    for seq in permutations(parts):
        for inds in combinations(range(len(seq)+1),r):
            if r>2:break
            mask=sum(1<<j for j in inds)
            w=[]
            for j in range(len(seq)+1):
                if mask>>j&1:w.append(PORT)
                if j<len(seq):w.append(seq[j])
            if all(arc(rows,i,j) for i,j in zip(w,w[1:])):answer+=1
    return answer

def inflate(rows,k):
    parent=[x for x in range(len(rows)) if x!=PORT]+[PORT]*k
    return [sum(1<<j for j,b in enumerate(parent)
      if (a!=b and arc(rows,a,b)) or (a==b==PORT and i<j)) for i,a in enumerate(parent)]

def check_tournament(rows):
    n=len(rows)
    for i in range(n):
        if arc(rows,i,i): raise AssertionError(("loop", i))
        for j in range(i):
            if arc(rows,i,j)+arc(rows,j,i)!=1: raise AssertionError(("pair", i, j))
    for i in range(n):
        seen=1<<i
        for _ in range(n):
            seen|=sum(1<<j for j in range(n) if any(seen>>h&1 and arc(rows,h,j) for h in range(n)))
        if seen!=(1<<n)-1: raise AssertionError(("strong", i, seen))

def hamilton(rows):
    n=len(rows);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for m in range(1,1<<n):
        for i in range(n):
            if not d[m][i]:continue
            outs=rows[i]&~m
            while outs:
                b=outs&-outs;outs-=b
                d[m|b][b.bit_length()-1]+=d[m][i]
    return sum(d[-1])
for CHILD,PORT,EXPECTED,COEFF in (
    ((252,445,504,32,136,208,154,8,249),3,((9,29,15),(27,147,105)),(15,-1)),
    ((350,300,376,112,98,449,2,351,88),6,((37,91,50),(207,581,350)),(50,-9)),
):
    SOURCE=tuple(row&127 for row in CHILD[:7])
    for rows,expected in zip((SOURCE,CHILD),EXPECTED):
        check_tournament(rows)
        if rows[PORT].bit_count()!=1:
            raise AssertionError("not a single exit")
        actual=tuple(skeleton(rows,r) for r in range(3))
        if actual!=expected:
            raise AssertionError((actual,expected))
        print(len(rows),actual)
    for k in range(1,6):
        a,b=inflate(SOURCE,k),inflate(CHILD,k)
        check_tournament(a);check_tournament(b)
        h,g=hamilton(a),hamilton(b)
        if h!=COEFF[0]*2**k+COEFF[1] or g!=7*(h-8):
            raise AssertionError((k,h,g))
        print(k,h,g,g//7)
