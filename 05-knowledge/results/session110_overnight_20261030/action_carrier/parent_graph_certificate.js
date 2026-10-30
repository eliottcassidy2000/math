function requireTrue(ok, message) {
  if (!ok) throw new Error(message);
}
function census(nonabelian) {
  const n = 21;
  const multiply = (x, y) => nonabelian
    ? 3 * ((Math.floor(x / 3) + [1, 2, 4][x % 3] * Math.floor(y / 3)) % 7) + (x % 3 + y % 3) % 3
    : (x + y) % n;
  const inverses = Array.from({ length: n }, (_, x) => {
    for (let y = 0; y < n; y++) if (multiply(x, y) === 0 && multiply(y, x) === 0) return y;
    throw new Error("Missing inverse");
  });
  for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) for (let c = 0; c < n; c++)
    requireTrue(multiply(multiply(a, b), c) === multiply(a, multiply(b, c)), "Associativity");
  const pairs = [], seen = new Set([0]);
  for (let x = 1; x < n; x++) if (!seen.has(x)) {
    const y = inverses[x];
    requireTrue(y !== x, "Nonidentity involution");
    pairs.push([x, y]); seen.add(x); seen.add(y);
  }
  const all = (1 << n) - 1;
  let prime = 0;
  const moduleSizes = {};
  for (let mask = 0; mask < (1 << pairs.length); mask++) {
    const S = new Set(pairs.map((pair, i) => pair[(mask >> i) & 1]));
    const rows = Array.from({ length: n }, (_, x) => {
      let row = 0;
      for (let y = 0; y < n; y++) if (S.has(multiply(inverses[x], y))) row |= 1 << y;
      return row;
    });
    for (let x = 0; x < n; x++) for (let y = 0; y < n; y++)
      requireTrue(x === y ? !(rows[x] & (1 << y)) : !!(rows[x] & (1 << y)) !== !!(rows[y] & (1 << x)), "Tournament");
    const proper = [];
    for (let v = 1; v < n; v++) {
      let M = 1 | (1 << v), changed = true;
      while (changed && M !== all) {
        changed = false;
        for (let z = 0; z < n; z++) if (!(M & (1 << z))) {
          const forward = rows[z] & M;
          if (forward !== 0 && forward !== M) { M |= 1 << z; changed = true; }
        }
      }
      if (M !== all) proper.push(M);
    }
    if (proper.length === 0) prime++;
    else {
      const sizes = [...new Set(proper.map(M => M.toString(2).replace(/0/g, "").length))].sort((a, b) => a - b).join(",");
      moduleSizes[sizes] = (moduleSizes[sizes] || 0) + 1;
    }
  }
  requireTrue(prime === (nonabelian ? 980 : 992), "Prime census");
  return { group: nonabelian ? "C7 semidirect C3" : "C21", labelledStates: 1024,
    prime, composite: 1024 - prime, properModuleSizes: moduleSizes };
}
const result = { status: "PASS", method: "Direct arcs and smallest-module closure; no subgroup recurrence",
  cases: [census(false), census(true)] };

console.log(JSON.stringify(result, null, 2));
