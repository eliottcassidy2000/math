exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'w34_twist_bias.py')).read().split("rng=random.Random(5)")[0])
MF={k:v.astype(float) for k,v in MI.items()}
def has_rank1_factor(word, maxlen=10):
    n=len(word)
    for i in range(1,n):            # skip the first letter (start effects)
        M=np.identity(dim)
        for j in range(i,min(n-1,i+maxlen)):   # skip the last letter
            M=MF[L.index(word[j])]@M
            sv=np.linalg.svd(M,compute_uv=False)
            if sv[0]>0 and sv[1]<1e-9*sv[0]: return True
    return False
rng=random.Random(9); n=int(sys.argv[1]); trials=int(sys.argv[2])
st={'strong, has rank-one factor':[0,0,0],'strong, bottleneck-free (<=10)':[0,0,0]}
for t in range(trials):
    w=[L[rng.randrange(4)] for _ in range(n)]
    if not is_strong(tournament_from_word(2,w)): continue
    h=Hmod(w); k='strong, has rank-one factor' if has_rank1_factor(w) else 'strong, bottleneck-free (<=10)'
    st[k][0]+=1; st[k][1]+=(h%3==0); st[k][2]+=(h%5==0)
for k,(c,z3,z5) in st.items():
    print(f'{k:32s} n={n}: samples {c:6d}  P(3|H) = {z3/max(c,1):.4f}   P(5|H) = {z5/max(c,1):.4f}')
