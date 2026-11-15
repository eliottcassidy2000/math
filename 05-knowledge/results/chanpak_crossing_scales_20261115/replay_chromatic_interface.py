from collections import defaultdict

def hvalues(rows):
    n=len(rows)
    full=(1<<n)-1
    totals=[1]+[0]*full
    paths=[[0]*n for _ in range(full+1)]
    for s in range(1,full+1):
        for i in range(n):
            if not (s>>i)&1: continue
            t=s^(1<<i)
            nb=t&rows[i]
            if not t: paths[s][i]=1
            while nb:
                bj=nb&-nb
                paths[s][i]+=paths[t][bj.bit_length()-1]
                nb^=bj
        totals[s]=sum(paths[s])
    return totals

def profile(rows):
    h=hvalues(rows)
    return h[-1],sum(a*h[(len(h)-1)^s] for s,a in enumerate(h))

def trans(n):
    return [(1<<n)-(1<<(i+1)) for i in range(n)]

def generator(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for m in range(1<<len(pairs)):
        rows=[0]*n
        for r,(i,j) in enumerate(pairs):
            if m>>r&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
        yield rows

def a_count(v):
    b,h,prev=0,1,0
    for s in v:
        c=int(s)
        if not c: h+=b+prev
        b=2*b+c
        prev=c
    return h

def signature(up,wp,v):
    u,qu=up
    h,q=wp
    x=qu-u
    t=u*h
    k=len(v)
    r=int(v[::-1],2)
    d=x*q
    l=d*(1<<k)+t*((1<<k)-r)
    c=x*(q-h)+2*t+d*int(v,2)+t*(a_count(v)-1)
    return l,c

def glue(word,W,U):
    k=len(word)
    rs=trans(k)+[0]
    for i,bit in enumerate(word):
        if bit=='1':rs[k]|=1<<i
        else:rs[i]|=1<<k
    l,r,s=len(rs),len(W),len(U)
    wm=((1<<r)-1)<<l
    um=((1<<s)-1)<<(l+r)
    return [a|wm|(um if i!=k else 0) for i,a in enumerate(rs)]+[(a<<l)|um for a in W]+[(a<<(l+r))|(1<<k) for a in U]

def main():
    C3=[4,1,2]
    S4=[8,1,3,6]
    packs=[(S4,trans(1),'0',118,37),
           (trans(3),trans(1),'10',59,39),
           (trans(3),trans(1),'110',117,97),
           (C3,C3,'0',234,99)]
    for U,W,v,l,c in packs:
        if signature(profile(U),profile(W),v)!=(l,c):raise AssertionError((v,l,c))
        for g in (1,2):
            actual=hvalues(glue('1'+'0'*g+v,W,U))[-1]
            if actual!=(l<<g)+c:raise AssertionError((v,g,actual))
    print('8 literal whole-glue checks on the two persistent +2 pairs')

    limit=99
    profs={n:{profile(rows) for rows in generator(n)} for n in range(1,7)}
    palette=defaultdict(set)
    for a in profs:
        for b in profs:
            if a+b>7:continue
            for u,qu in profs[a]:
                for h,q in profs[b]:
                    x=qu-u
                    G=x*(q-h)+2*u*h
                    d=x*q
                    for value in range((limit-G)//d+1):
                        v=bin(value)[2:]
                        l,c=signature((u,qu),(h,q),v)
                        if c<=limit:palette[c].add(l//(l&-l))
    links={c:palette[c]&palette[c+2] for c in range(1,limit-1,2) if palette.get(c,set())&palette.get(c+2,set())}
    metrics=len(palette),sum(map(len,palette.values())),links
    if metrics!=(34,95,{37:{59},97:{117}}):raise AssertionError(metrics)
    print('Core order 1–6 profile counts:',[len(profs[n]) for n in profs])
    print('Checksum lists and incidences through 99:',metrics[:2])
    print('Exact adjacent intersections:',{c:sorted(links[c]) for c in links})
if __name__=='__main__':main()
