# Barrier walk W_k = sum_{i<=k} (2 v_i - 3) along the Antihydra odd orbit of 27; b at odd terms = 5 + W_k.  Also the
# golden-ratio hitting probability under the Bernoulli valuation model and Collatz comparison numbers.
import sys, json, math
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
a, b = 8, 0
W = 0; minW = 0; argminW = 0; k = 0; step = 0
hist = {}
min_b_after_first_odd = None
while step < N:
    if a & 1:
        if b == 0: print('HALT', step); sys.exit()
        v = ((3*a-1) & -(3*a-1)).bit_length() - 1
        a = (3*a-1) >> 1; b -= 1
        if k > 0 or step > 3:
            pass
        k += 1
        if k >= 2:   # W uses the valuation of the previous odd term
            W += 2*v_prev - 3
            if W < minW: minW, argminW = W, k
        v_prev = v
        hist[v] = hist.get(v,0)+1
        if min_b_after_first_odd is None or b < min_b_after_first_odd: min_b_after_first_odd = b
    else:
        a = (3*a) >> 1; b += 2
    step += 1
# model: r = P(walk with steps -1 (1/2), +1 (1/4), +3 (1/8), ... ever decreases by 1) solves r^3 - 2r + 1 = 0 -> r = (sqrt5-1)/2
phi = (1+5**0.5)/2
r = 1/phi
print(json.dumps({'steps': N, 'odd_terms': k, 'W_final': W, 'W_min': minW, 'at_odd_index': argminW, 'min_b_after_first_odd': min_b_after_first_odd,
                  'model_r': r, 'check_cubic': r**3-2*r+1, 'P_halt_model_from_B0_5': r**6, 'phi^-6': phi**-6,
                  'expected_step': sum((2*v-3)*2.0**-v for v in range(1,60)),
                  'collatz_drift_log2_per_odd_step': math.log2(3) - 2, 'antihydra_slope_threshold': 1.5, 'collatz_slope_threshold': math.log2(3),
                  'empirical_v_freq': {v: round(c/k,5) for v,c in sorted(hist.items())[:6]}}, indent=0))
