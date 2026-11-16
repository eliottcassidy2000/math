# Exact model of TT_m -> B, one head x with independent contacts to both.
CONFIGS = (
    (73,(0,5,17,7,43,15,63,127)),
    (89,(0,1,11,35,15,87,31,127)),
    (77,(2,8,3,21,7,31,63,127)),
    (105,(0,1,11,19,39,79,31,127)),
)

def literal_count(m, word, x, b):
    """Independent exact path enumeration by split at the head.
    Prefix/suffix inside the upper tape are increasing in position.
    DP implements every possible path inside each subset of B."""
    p=len(b)
    full=(1<<p)-1
    ends=[[0]*p for _ in range(full+1)]
    starts=[[0]*p for _ in range(full+1)]
    for i in range(p):
        ends[1<<i][i]=1
        starts[1<<i][i]=1
    for S in range(1,full+1):
        for i in range(p):
            if not (S>>i&1):
                continue
            for j in range(p):
                if S>>j&1:
                    continue
                V=S|(1<<j)
                if b[i]>>j&1:
                    ends[V][j]+=ends[S][i]
                else:
                    starts[V][j]+=starts[S][i]
    cases=[0]*8
    tape_mask=(1<<m)-1
    word=int(format(word,f'0{m}b')[::-1],2)
    for S in range(tape_mask+1):
        Z=tape_mask^S
        incoming_ok=not S or not(word>>(S.bit_length()-1)&1)
        after_nonempty=bool(Z)
        outgoing_ok=not Z or bool(word& (Z&-Z))
        key=incoming_ok+2*after_nonempty+4*outgoing_ok
        cases[key]+=1
    answer=0
    for S in range(full+1):
        Z=full^S
        e=sum(ends[S][i] for i in range(p) if not x>>i&1)
        s=sum(starts[Z][i] for i in range(p) if x>>i&1)
        h=sum(ends[Z])
        if Z==0:h=s=1
        for k, t in enumerate(cases):
            before=e if S else (k&1)
            after=(h if k&4 else 0) if k&2 else s
            answer+=t*before*after
    return answer

def upper_half(A):
    N=p=b=0
    for z in format(A,'b'):
        if z=='0':N+=(b+p)//2
        p=int(z);b=2*b+p
    return N

def offset(A):
    return 9*upper_half(A)+245*(A//2)+127*(A%2)

def run(m,words,base,configs=CONFIGS):
    for j,A in enumerate(words):
        vals=[literal_count(m,A,x,b) for x,b in configs]
        assert vals==[base+2*i+6*j for i in range(4)]
        print(m,A,*vals)

if __name__=='__main__':
    run(7,(108,111),27323)
    run(15,(28085,27854,28109,28157),6962193)

    second=((15,(0,5,17,7,43,79,31,127)),
            (31,(2,8,3,21,39,15,63,127)),
            (47,(0,33,11,19,7,29,63,127)),
            (55,(2,8,3,21,39,15,63,127)))
    run(8,(129,135),16611,second)
    run(19,(258408,284439,265290,262716),34801935,second)
