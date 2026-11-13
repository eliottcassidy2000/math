"""Digits of every vertex: sigma_a = B_a - h (sole-exit / out-birth convention, THM-4865) and
   tau_a = E_a - h (in-birth convention: z beaten only by a).  C2 at a gives 3h + 2*digit."""
import sys; sys.path.insert(0,'scripts')
def subset_tables(adj, n):
    N = 1 << n
    # end[S][v]: paths covering S ending at v ; start[S][v]: starting at v
    end = [[0]*n for _ in range(N)]; start = [[0]*n for _ in range(N)]
    for v in range(n): end[1<<v][v] = 1; start[1<<v][v] = 1
    for S in range(1, N):
        for v in range(n):
            c = end[S][v]
            if c:
                out = adj[v] & ~S
                while out:
                    u = (out & -out).bit_length()-1; out &= out-1
                    end[S|(1<<u)][u] += c
    radj = [0]*n
    for v in range(n):
        for u in range(n):
            if adj[v]>>u&1: radj[u] |= 1<<v
    for S in range(1, N):
        for v in range(n):
            c = start[S][v]
            if c:
                inn = radj[v] & ~S
                while inn:
                    u = (inn & -inn).bit_length()-1; inn &= inn-1
                    start[S|(1<<u)][u] += c
    h = [sum(end[S]) for S in range(N)]; h[0] = 1
    return end, start, h
def digits(adj, n):
    end, start, h = subset_tables(adj, n)
    N = 1 << n; full = N-1; H = h[full]
    sig = []; tau = []
    for a in range(n):
        B = 0; E = 0
        rest = full & ~(1<<a)
        A = rest
        while True:
            comp = full & ~A   # contains a
            B += h[A] * start[comp][a]
            E += end[comp][a] * h[A]
            if A == 0: break
            A = (A-1) & rest
        sig.append(B - H); tau.append(E - H)
    return H, sig, tau
if __name__ == '__main__':
    from tern_check import parse, blowup
    for s,n,u in [('1110011111111111001111111111',8,3),('111001111111111111001111111111111110',9,3),('111010111111111111100111111111111111',9,3)]:
        adj = parse(s,n)
        print('skeleton', s, 'u', u)
        for m in range(1,5):
            A,n2,block,newid = blowup(adj,n,u,m)
            H,sig,tau = digits(A,n2)
            inv = {v:k for k,v in newid.items()}
            lab = [str(inv.get(x, 'u%d'%(x-block[0]+1))) for x in range(n2)]
            print('  j=%d H=%d'%(m,H), 'sigma:', dict(zip(lab,sig)), '\n        tau:', dict(zip(lab,tau)))
