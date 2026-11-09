"""S19 audit: miss density per window for a cover2 output (used for the order <= 6 bank to 10^7).
Usage: python density_windows.py COVER2_OUTPUT"""
import sys
s = open(sys.argv[1]).read()
nums = list(map(int, s.split('misses:')[1].split('count')[0].split()))
print("misses listed:", len(nums), "(permanent gaps included)")
for a, b in [(1, 10**5), (10**5, 10**6), (10**6, 2 * 10**6), (2 * 10**6, 5 * 10**6), (5 * 10**6, 10**7)]:
    c = sum(1 for v in nums if a <= v < b)
    print(f"  [{a}, {b}): {c} misses, {c / ((b - a) / 2):.2e} per odd integer")
