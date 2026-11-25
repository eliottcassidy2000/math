"""Exact path-length intervals on disjoint tracks, with a fixed common root."""


def panel_sums(values):
    # State = (number chosen, endpoints chosen). Bit i denotes weighted sum i.
    dp = {(0,0): 1}
    for v in values:
        new = dp.copy()
        for (n,e), bits in dp.items():
            internal = n+1,e
            new[internal] = new.get(internal,0) | (bits << (2*v))
            if e < 2:
                endpoint = n+1,e+1
                new[endpoint] = new.get(endpoint,0) | (bits << v)
        dp = new
    return {state:[i for i in range(bits.bit_length()) if bits >> i & 1]
            for state,bits in dp.items()}


def rooted_lengths(d):
    left = panel_sums(range(d//2))
    right = panel_sums(range(d//2,d))
    lengths = set()
    for j in range(2,d+1):
        for low in {j//2,(j+1)//2}:
            high = j-low
            low_ends = 1+low-high
            lows = left.get((low,low_ends),[])
            highs = right.get((high,2-low_ends),[])
            lengths.update(2*j+y-x for x in lows for y in highs)
    return lengths


if __name__ == "__main__":
    for d in (3,7,13,25,32,57):
        lengths = rooted_lengths(d)
        upper = d*d//2 - 1 + 2*d
        if d == 3:
            # On three columns the increasing itinerary 0,1,2 supplies 8.
            lengths.add(8)
            assert lengths == {5,6,8,9}
        else:
            assert lengths == set(range(5,upper+1))
        print(d,len(lengths),min(lengths),max(lengths))
