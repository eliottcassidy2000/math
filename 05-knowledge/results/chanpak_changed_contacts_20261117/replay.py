from itertools import permutations, combinations
CHILD=(250,228,249,242,194,16,32,96)
SOURCE=tuple(row&63 for row in CHILD[:6])
PORT=5

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
for rows in (SOURCE,CHILD):
    check_tournament(rows)
    vals=tuple(skeleton(rows,r) for r in range((len(rows)+2)//2))
    expected=(9,15,5) if rows==SOURCE else (9,49,35)
    if vals[:3]!=expected or any(vals[3:]):raise AssertionError(vals)
    print(len(rows), vals)
for k in range(1,8):
    a,b=inflate(SOURCE,k),inflate(CHILD,k)
    check_tournament(a);check_tournament(b)
    h,g=hamilton(a),hamilton(b)
    if (h,g)!=(5*2**k+5,35*2**k-21):raise AssertionError((k,h,g))
    print(k,h,g,g//7)

def cut(w):
    if not w: return 0
    r,b = divmod(w,2)
    return cut(r) + ((r+(r&1))//2 if not b else 0)

# Different construction: build D outright, with q at index m.
for w in range(106):
    m = max(1,w.bit_length())
    out = [0]*(m+1)
    for i in range(m):
        for j in range(i+1,m):
            out[i] |= 1<<j
        if w >> (m-1-i) & 1:
            out[m] |= 1<<i
        else:
            out[i] |= 1<<m
    if hamilton(out) != 1+2*cut(w):
        raise AssertionError(("binary cut",w))
    for t in (2,4):
        b = w&1
        if -6*(2*cut(w)+b)+(1-b)+(7*t-13)+(t-1)*w == -56:
            raise AssertionError(("transitive cover inverse",t,w))
print("TT3 transitive-deletion exclusion: all bounded exact addresses checked")

for branch in (1,2):
    vals=[]
    for t in (2,4,8,16):
        for j in range(1,56,2):
            for z in range(0 if branch==1 else 1,56):
                cost=((4*t+1)*j+5*(t-1)*z if branch==1 else 5*j+3*(t-1)*z)
                if cost<=56:
                    vals.append((t+1)*j+t*z)
    peak=max(vals)
    if peak!=(21 if branch==1 else 71):
        raise AssertionError(("source bound", branch, peak))
    print("loose source bound",branch,peak)
