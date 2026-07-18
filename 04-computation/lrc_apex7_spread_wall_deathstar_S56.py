# Why the covering argument STALLS at the classification: the apex-7 wall.
# To bound the spread I must rule out a "small cluster S (k speeds) + large outliers L (13-k speeds)".
# L must cover the good set of S (measure mu_S). L's danger density = (13-k)/7.
# Contradiction needs mu_S > (13-k)/7. But union bound gives mu_S >= 1-k/7, which is <=0 for k>=7.
for k in range(1,14):
    L=13-k
    dens=L/7                      # large-cluster danger density (measure of L-danger)
    ub=max(0.0,1-k/7)            # union-bound lower bound on S-good measure
    note=""
    if L<=6 and ub<=dens: note=" <- union bound gives mu_S>=0 only; covering NOT excluded (apex-7 wall)"
    if L>=7: note=" (>=7 large speeds: density>=1, they can tile - no obstruction)"
    print(f" k={k:2d} small, {L:2d} large: L-danger density={dens:.3f}, union-bd mu_S>= {ub:.3f}{note}")
print("\nThe wall: at k=7 (7 small vs 6 large), 7*(1/7)=1 so the union bound on mu_S collapses to 0,")
print("and 6 large speeds have density 6/7<1 (can't tile) yet CAN cover a small-measure good set.")
print("Deciding whether they do = the apex-7 / Fraenkel tiling core = the LRC(14) hard problem itself.")
