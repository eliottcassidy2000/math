"""Discover a triangulation; output is untrusted until exact verification."""
import json
from collections import Counter
from scipy.spatial import Delaunay
from explore import N, GOOD
import numpy as np

vertices = [0]+GOOD
points = np.array([[s >> i & 1 for i in range(N)] for s in vertices])
tri = Delaunay(points)
simplices = sorted(set(tuple(sorted(map(int, row))) for row in tri.simplices))
hist = Counter()
kept = []
for row in simplices:
    a = points[list(row)]
    det = round(np.linalg.det(a[1:]-a[0]))
    hist[abs(det)] += 1
    if det:
        kept.append(row)
print(json.dumps({'simplices': len(simplices), 'det_histogram_discovery_only': dict(hist)}))
from pathlib import Path
Path(__file__).with_name('triangulation.json').write_text(json.dumps({'vertex_masks': vertices, 'simplices': kept}, separators=(',', ':'))+'\n', encoding='utf-8')
