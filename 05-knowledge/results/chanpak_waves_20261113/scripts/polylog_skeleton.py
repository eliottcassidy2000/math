"""Wave 25 (THM-5172 (A)): along a fixed block skeleton the counts grow exponentially in the J-runs, so the words
J^a0 JOZ J^a1 (from the root state) with count <= X number O((log X)^2)."""
def step(s,m):
    H,E,L,h=s
    if m=='O': return (H,0,2*L,h)
    if m=='Z': return (H+L-E,h,2*L-h,h)
    return (H+L-E,H-E,H+L,H)
def run(s,w):
    for ch in w: s=step(s,ch)
    return s
for X in (10**3,10**6,10**9,10**12,10**15):
    cnt=0
    for a0 in range(200):
        s0=run((1,1,1,1),'J'*a0)
        if s0[0]>X: break
        s1=run(s0,'JOZ')
        for a1 in range(200):
            s2=run(s1,'J'*a1)
            if s2[0]>X: break
            cnt+=1
    print(f"X = 10^{len(str(X))-1}: words J^a0 JOZ J^a1 with count <= X: {cnt}")
