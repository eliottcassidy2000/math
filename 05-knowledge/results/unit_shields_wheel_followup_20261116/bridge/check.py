from itertools import combinations, product

def require(ok, why):
    if not ok:
        raise AssertionError(why)

def prof(w):
    h,A,p=1,0,0
    for e in map(int,w):
        h,A,p=(h,2*A+1,1) if e else (h+A+p,2*A,0)
    return h,A+1

def rc(w):
    return ''.join('1' if b=='0' else '0' for b in reversed(w))

def spine_word(u,v):
    # Bits of the two apices' out-neighbor indicators in forward spine order:
    # y loses all U-internal spine arcs; x beats all V-internal spine arcs.
    return [(int(a),0) for a in u]+[(1,int(b)) for b in v]

def actual_h(bits):
    # Literal Hamilton path enumeration as a 3-run automaton. Every
    # full path has a choice of apex order; spine vertices are assigned
    # to before/between/after. The individual runs increase in spine
    # order, so a start condition and final last-vertex condition are
    # sufficient. No use of odd cycles or bridge formulas.
    answer=0
    for alpha,beta in [(0,1),(1,0)]:
        table={(0,0,0): 1}
        # Per run: 0 empty; 1 occupied with invalid current endpoint;
        # 2 occupied with valid current endpoint.
        for flags in bits:
            nt={}
            for key,n in table.items():
                for r in range(3):
                    if key[r]==0 and ((r==1 and not flags[alpha]) or
                                      (r==2 and not flags[beta])):
                        continue
                    closed=(r==0 and not flags[alpha] or
                            r==1 and not flags[beta] or r==2)
                    new=key[:r]+(1+closed,)+key[r+1:]
                    nt[new]=nt.get(new,0)+n
            table=nt
        for key,count in table.items():
            if 1 in key:
                continue
            if key[1]==0 and (alpha,beta)!=(1,0):
                continue
            answer+=count
    return answer

def direct_out(bits):
    m=len(bits)
    rows=[set(range(i+1,m)) for i in range(m)]+[set(),{m}]
    for i,pair in enumerate(bits):
        for k,b in enumerate(pair):
            if b:
                rows[m+k].add(i)
            else:
                rows[i].add(m+k)
    return rows

def check_strong(bits):
    rows=direct_out(bits)
    for start in range(len(rows)):
        q=[start];seen={start}
        for a in q:
            for b in rows[a]-seen:
                q.append(b);seen.add(b)
        require(len(seen)==len(rows),'not strong')

def check():
    u0,u1='1101011000','1111111110'
    require(prof(u0)==(769,857) and prof(u1)==(513,1023),'left profiles')
    require([i for i,(a,b) in enumerate(zip(u0,u1)) if a!=b]==[2,4,7,8],'edits')
    for r in range(1,4):
        for j in combinations(range(10),r):
            for signs in product((-1,1),repeat=r):
                require(sum(s*2**x for s,x in zip(signs,j))!=166,'smaller edit')
    checked=0
    for t in list(range(1,25))+[40,60,80,120]:
        D=2**t
        for prefix,m,left in [('11',3,520569),('101',5,867615)]:
            w=prefix+'0'*t+'1011010'
            j,Q=prof(w)
            require((j,Q)==(m*83*D+59,m*128*D+91),('right profile',w))
            require(-128*j+83*Q==1,('unit',w))
            v=rc(w)
            want=left*D+123357
            for u,h in [(u0,want),(u1,want+2)]:
                bits=spine_word(u,v)
                require(actual_h(bits)==h,('literal count',u,v))
                if t<=3:
                    check_strong(bits)
                checked+=1
    for n in range(2,15):
        for tail in product('01',repeat=n-2):
            w='1'+''.join(tail)+'0'
            h,Q=prof(w)
            wanted=any(w==prefix+'0'*(n-7-len(prefix))+'1011010'
                       for prefix in ('11','101') if n>7+len(prefix))
            require((-128*h+83*Q==1)==wanted,('inverse',w))
    print(f'{checked} direct physical three-run counts; unit inverse and edit size replayed.')

if __name__=='__main__':
    check()
