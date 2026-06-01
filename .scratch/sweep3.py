from lrc_coloring_core import *
import collections

# WHY chi=omega: arcs each length 1/n, n arcs, total length =1. A clique = points all
# pairwise within 1/n. Claim: omega(G(t)) = max over arcs of "how many points lie in a
# half-open window"... Let's TEST a cleaner characterization:
#   omega = max_{window of length 1/n, half-open} #points in it ?  (sliding window)
# and chi = same? If arcs total length = circumference exactly, the circular-arc graph is
# actually an interval graph after cutting at any gap, OR chi=omega by a covering argument.

def sliding_window_max(speeds, n, t):
    pts = sorted(positions(speeds, n, t) if False else positions(speeds, t))
    thr = F(1,n)
    N=len(pts)
    # points are NOT necessarily distinct; treat as multiset on circle
    arr = sorted(pts)
    best=0
    for i in range(N):
        # window [arr[i], arr[i]+thr) , count points within (strictly within 1/n means dist<thr
        # but two points can be adjacent if circ_dist<thr. A clique needs PAIRWISE <thr.
        # pairwise<thr for a set in a window: max gap between extreme members <thr in circular sense
        c=0
        for j in range(N):
            d=circ_dist(arr[i],arr[j])
            # only count those on the "forward" side within thr to define window
        # simpler: count max set pairwise within thr -> that's the true clique. skip window heuristic
        pass
    return None

# Direct check: does omega(G) equal max number of points in any OPEN arc of length 1/n?
def max_points_in_open_arc(speeds, n, t):
    pts = sorted(positions(speeds, t))
    thr=F(1,n)
    N=len(pts)
    # doubled for wrap
    ext = pts + [p+1 for p in pts]
    best=0
    for i in range(N):
        # arc (pts[i]-eps, pts[i]+thr-eps): count points x with pts[i] <= x < pts[i]+thr (mod)
        c=0
        base=pts[i]
        for x in ext:
            if base <= x < base+thr:
                c+=1
        best=max(best,c)
    return best

bad=0; tested=0
for speeds in [(1,2,3,4),(1,3,4,7),(2,3,5,7),(1,2,4,8),(1,2,3,4,5),(1,3,4,5,9),(1,3,5,7,9,11)]:
    n=len(speeds)+1
    crit=critical_times(speeds,n); mids=cell_midpoints(crit)
    for t in mids:
        tested+=1
        w=clique_number_arc(speeds,n,t)
        a=max_points_in_open_arc(speeds,n,t)
        if w!=a:
            bad+=1
            if bad<=5: print("MISMATCH",speeds,t,"omega",w,"arc",a)
print("clique == max-points-in-arc-of-length-1/n :", "ALL MATCH" if bad==0 else f"{bad} mismatch", "of",tested)
print()

# Time-coloring obstruction (part 2): N(t)=#runners within 1/n of observer.
# observer-arc = open arc of length 2*(1/n)=2/n centered at 0? No: runner within 1/n of obs
# means runner position in (-1/n,1/n) => arc length 2/n. So N(t)=# runners in that window.
# LRC: some t has N(t)=0. The OBSERVER is in EVERY clique it's part of. observer isolated
# (degree 0) <=> N(t)=0 <=> observer's own color is "free".
# Coloring obstruction question: relate N(t)=0 to chi. observer isolated => removing it,
# chi(G) = chi(G - obs). When N(t)=0 the danger graph restricted... let's tabulate.
print("relation N(t)=0  vs  observer isolated vs chi:")
for speeds in [(1,2,3,5),(2,3,5,7),(1,2,3,4,6)]:
    n=len(speeds)+1
    crit=critical_times(speeds,n); mids=cell_midpoints(crit)
    for t in mids:
        Nt=observer_danger_count(speeds,n,t)
        if Nt==0:
            Ncount,edges,pts=danger_graph(speeds,n,t)
            obs_deg=sum(1 for e in edges if 0 in e)
            k,_=chromatic_number(speeds,n,t)
            print(f"  {speeds} t={t} N=0 obs_deg={obs_deg} chi={k}")
            break
