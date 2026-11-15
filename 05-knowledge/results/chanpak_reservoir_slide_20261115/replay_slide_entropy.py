"""Exact finite branching bounds in two active return components."""
from replay_slide_cocycle import live,require

for h,Q,O in [(13,104,(18,18,0,1)),(45,1440,(517,818,1,1))]:
    F,C,g=live(h,Q,O)
    degrees={u:sum(v in C for v,_ in F[u]) for u in C}
    require(max(degrees.values())<=2 and min(degrees.values())<2 and g>0)
    f={u:1 for u in C}
    for j in range(1,13):
        f={u:sum(f[v] for v,_ in F[u] if v in C) for u in C}
    M=max(f.values())
    require(0<M<2**12)
    print('h',h,'internal branch histogram',[list(degrees.values()).count(i) for i in range(3)],
          'max length-12 words',M,'out of',2**12, 'power exponent bound',f"{(M-1).bit_length()}/12")
