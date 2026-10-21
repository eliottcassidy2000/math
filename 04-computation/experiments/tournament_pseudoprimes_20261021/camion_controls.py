# camion_controls.py (corrected after the independent audit, MISTAKE-651)
# (i)  strong semicomplete digraphs (3^C(m,2), linear condensation): quasi-period p-1 with (3/p), m >= 1, p >= 5.
# (ii) connected labelled graphs (2^C(m,2), unordered components): quasi-period p-1 with (2/p) for m >= 2
#      (audit Lemma C; Mani-Stones 2016); the law FAILS only at m = 1.  The pre-audit version tested from m = 1
#      and wrongly concluded "no period".
# (iii) genuine negative control: the rooted count m*s(m) has quasi-period p(p-1), not p-1.
def strong_mod(base,M,p):
    T=[pow(base,j*(j-1)//2,p) for j in range(M+1)]; s=[0]*(M+1); row=[1]
    for k in range(1,M+1):
        row=[1]+[(row[i-1]+row[i])%p for i in range(1,k)]+[1]
        s[k]=(T[k]-sum(row[i]*s[i]*T[k-i] for i in range(1,k)))%p
    return s
def connected_mod(M,p):
    g=[pow(2,j*(j-1)//2,p) for j in range(M+1)]; c=[0]*(M+1); row=[1]
    for m in range(1,M+1):
        if m>1: row=[1]+[(row[i-1]+row[i])%p for i in range(1,m-1)]+[1]
        c[m]=(g[m]-sum(row[k-1]*c[k]*g[m-k] for k in range(1,m)))%p
    return c
def leg(a,p): return 1 if pow(a,(p-1)//2,p)==1 else -1
def law(seq,L,e,p,m0,M): return all((seq[m+L]-e*seq[m])%p==0 for m in range(m0,M-L+1))
for p in [5,7,11,13,17,19,23,29,31]:
    M=3*p*(p-1)+10 if p<=13 else 6*p
    s3=strong_mod(3,M,p); e3=leg(3,p)
    cg=connected_mod(M,p); e2=leg(2,p)
    s2=strong_mod(2,M,p); rooted=[(m*s2[m])%p for m in range(M+1)]
    out=[f'p={p}:',
         f'semicomplete (3/p)={e3:+d} from m>=1: {law(s3,p-1,e3,p,1,M)}',
         f'connected (2/p) from m>=2: {law(cg,p-1,e2,p,2,M)}, from m>=1: {law(cg,p-1,e2,p,1,M)}',
         f'rooted m*s(m) period p-1: {law(rooted,p-1,e2,p,1,M)}']
    if p<=13: out.append(f'rooted period p(p-1): {law(rooted,p*(p-1),e2,p,1,M)}')
    print('  '.join(out))
