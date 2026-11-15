exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'w34_twist_bias.py')).read().split("rng=random.Random(5)")[0])
M3={k:(v%3) for k,v in MI.items()}
def rank_le1_mod3(M):
    rows=[r for r in M if r.any()]
    if not rows: return True
    r0=rows[0]
    for r in rows[1:]:
        if not ((r==r0).all() or (r==(2*r0)%3).all()): return False
    return True
def has_mod3_bottleneck(word,maxlen=10):
    n=len(word)
    for i in range(1,n):
        M=np.identity(dim,dtype=np.int64)
        for j in range(i,min(n-1,i+maxlen)):
            M=(M3[L.index(word[j])]@M)%3
            if rank_le1_mod3(M): return True
    return False
rng=random.Random(9); n=int(sys.argv[1]); trials=int(sys.argv[2])
st={'strong, mod-3 bottleneck':[0,0],'strong, no mod-3 bottleneck':[0,0]}
for t in range(trials):
    w=[L[rng.randrange(4)] for _ in range(n)]
    if not is_strong(tournament_from_word(2,w)): continue
    h=Hmod(w); k='strong, mod-3 bottleneck' if has_mod3_bottleneck(w) else 'strong, no mod-3 bottleneck'
    st[k][0]+=1; st[k][1]+=(h%3==0)
for k,(c,z3) in st.items(): print(f'{k:30s} n={n}: samples {c:6d}  P(3|H) = {z3/max(c,1):.4f} (uniform 0.3333)')
