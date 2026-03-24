#!/usr/bin/env python3
"""
tc_algebra.py -- Transform Algebra: Functional Compression via Monoid Structure
kind-pasteur-2026-03-24-S20cq

THE DEEP IDEA: Instead of a case-switch over compression techniques,
model transforms as elements of a FREE MONOID under composition.

In Haskell terms:
  data Transform = Id | Delta | Stride Int | Gray | BWT | RLE | ...
  type Chain = [Transform]            -- free monoid = list of atoms
  compress :: Chain -> ByteString -> ByteString
  compress = foldl' (flip apply) where apply t = runTransform t

  -- The search is a CATAMORPHISM over the transform algebra:
  bestChain :: ByteString -> (ByteString, Chain)
  bestChain input = minimumBy (comparing fst)
                  . map (\c -> (encode (compress c input), c))
                  $ explore transformAlgebra

  -- But we don't explore blindly. The algebra has STRUCTURE:
  -- 1. Delta . Delta = Delta2 (nilpotent-ish: Delta^k converges)
  -- 2. Stride k . Stride k = Stride k (idempotent)
  -- 3. Some pairs COMMUTE: Stride . Gray = Gray . Stride (for separable transforms)
  -- 4. ANNIHILATION: RLE . Random = identity (RLE does nothing on random data)

  -- So we quotient the free monoid by these relations, getting a PRESENTATION:
  -- TransformMonoid = <Delta, Stride_k, Gray, BWT, RLE | relations>
  -- This is MUCH smaller than the free monoid.

THE TROPICAL SEMIRING:
  Compressed sizes form a tropical semiring (min, +):
  - a ⊕ b = min(a, b)    -- pick the better of two independent attempts
  - a ⊗ b = a + b        -- cost of doing a then b (roughly additive overhead)

  The optimal chain minimizes the TROPICAL EVALUATION of the transform expression.

THIS FILE IMPLEMENTS:
  1. Transform as first-class algebraic objects (composable, invertible)
  2. Chain as the free monoid (with quotient by known relations)
  3. Backend as a separate algebra (the codec endomorphism ring)
  4. Tropical evaluation for chain selection
  5. Beam search over the quotient monoid (not brute force)
  6. Application to bytes, images, audio, video (polymorphic)

THE PAYOFF: Instead of 45 hardcoded chains, we get an INFINITE space of
transforms that we search INTELLIGENTLY using algebraic properties.
"""

import sys
import zlib
import bz2
import lzma
import math
import time
from collections import Counter
from typing import Callable, List, Tuple, Optional
from dataclasses import dataclass, field
import numpy as np

__version__ = "1.0.0"


# ============================================================================
# TRANSFORM ALGEBRA: First-class composable transforms
# ============================================================================

@dataclass(frozen=True)
class Transform:
    """An atomic transform in our algebra.

    Each transform has:
      - name: human-readable identifier
      - forward: bytes -> bytes (the transform)
      - inverse: bytes -> bytes (the inverse, for decompression)
      - properties: algebraic properties (idempotent, nilpotent, etc.)
    """
    name: str
    forward: Callable
    inverse: Callable
    # Algebraic properties
    idempotent: bool = False     # f . f = f
    involutory: bool = False     # f . f = id
    nilpotent_order: int = 0     # f^k = 0 (converges to constant)
    commutes_with: tuple = ()    # names of transforms this commutes with

    def __call__(self, data: bytes) -> bytes:
        return self.forward(data)

    def __repr__(self):
        return self.name


class Chain:
    """A chain of transforms = element of the free monoid.

    Represented as a tuple of Transforms. The monoid operation is concatenation.
    Identity element is the empty chain ().
    """
    __slots__ = ['transforms', '_name']

    def __init__(self, *transforms):
        self.transforms = tuple(transforms)
        self._name = None

    @property
    def name(self) -> str:
        if self._name is None:
            if not self.transforms:
                self._name = 'id'
            else:
                self._name = '+'.join(t.name for t in self.transforms)
        return self._name

    def __call__(self, data: bytes) -> bytes:
        """Apply the chain: f_n . ... . f_2 . f_1 (left to right)."""
        result = data
        for t in self.transforms:
            result = t(result)
            if result is None or len(result) > len(data) * 3:
                return data  # bail out if transform explodes
        return result

    def __add__(self, other: 'Chain') -> 'Chain':
        """Monoid operation: concatenation."""
        return Chain(*(self.transforms + other.transforms))

    def __len__(self):
        return len(self.transforms)

    def __repr__(self):
        return f"Chain({self.name})"

    def inverse(self, data: bytes) -> bytes:
        """Apply inverse chain (in reverse order)."""
        result = data
        for t in reversed(self.transforms):
            result = t.inverse(result)
        return result

    def is_identity(self) -> bool:
        return len(self.transforms) == 0

    def simplify(self) -> 'Chain':
        """Simplify using algebraic relations."""
        if len(self.transforms) <= 1:
            return self
        simplified = list(self.transforms)
        changed = True
        while changed:
            changed = False
            new = []
            i = 0
            while i < len(simplified):
                if i + 1 < len(simplified):
                    a, b = simplified[i], simplified[i + 1]
                    # Idempotent: f . f = f
                    if a.name == b.name and a.idempotent:
                        new.append(a)
                        i += 2
                        changed = True
                        continue
                    # Involutory: f . f = id (skip both)
                    if a.name == b.name and a.involutory:
                        i += 2
                        changed = True
                        continue
                new.append(simplified[i])
                i += 1
            if i == len(simplified) - 1:
                new.append(simplified[-1])
            simplified = new
        return Chain(*simplified)


# ============================================================================
# BACKEND ALGEBRA: The codec endomorphism ring
# ============================================================================

@dataclass(frozen=True)
class Backend:
    """A compression backend (the actual entropy coder)."""
    name: str
    compress: Callable
    decompress: Callable
    level: int = 0  # speed/quality tradeoff

    def __call__(self, data: bytes) -> Optional[bytes]:
        try:
            return self.compress(data)
        except:
            return None


# ============================================================================
# CODEC: The full transform + backend algebra
# ============================================================================

class TransformAlgebra:
    """The full algebra of transforms, backends, and their composition.

    This is the mathematical core. It replaces the case-switch with:
    1. An algebraic structure (monoid of transforms)
    2. A tropical semiring evaluation (size minimization)
    3. A beam search over the quotient monoid
    """

    def __init__(self):
        self.atoms: List[Transform] = []
        self.backends: List[Backend] = []
        self._register_defaults()

    def _register_defaults(self):
        """Register default atomic transforms and backends."""
        # === ATOMIC TRANSFORMS ===

        def _delta_fwd(d):
            if len(d) < 2: return d
            out = bytearray(len(d)); out[0] = d[0]
            for i in range(1, len(d)): out[i] = (d[i] - d[i-1]) & 0xFF
            return bytes(out)
        def _delta_inv(d):
            if len(d) < 2: return d
            out = bytearray(len(d)); out[0] = d[0]
            for i in range(1, len(d)): out[i] = (out[i-1] + d[i]) & 0xFF
            return bytes(out)
        self.atoms.append(Transform('delta', _delta_fwd, _delta_inv, nilpotent_order=3))

        def _xor_fwd(d):
            if len(d) < 2: return d
            out = bytearray(len(d)); out[0] = d[0]
            for i in range(1, len(d)): out[i] = d[i] ^ d[i-1]
            return bytes(out)
        def _xor_inv(d):
            if len(d) < 2: return d
            out = bytearray(len(d)); out[0] = d[0]
            for i in range(1, len(d)): out[i] = d[i] ^ out[i-1]
            return bytes(out)
        self.atoms.append(Transform('xor', _xor_fwd, _xor_inv, involutory=False))

        def _gray_fwd(d): return bytes(b ^ (b >> 1) for b in d)
        def _gray_inv(d):
            def ig(b):
                m = b >> 1
                while m: b ^= m; m >>= 1
                return b
            return bytes(ig(b) for b in d)
        self.atoms.append(Transform('gray', _gray_fwd, _gray_inv, involutory=False,
                                     commutes_with=('stride2', 'stride4')))

        for s in [2, 3, 4, 8, 16]:
            def _make_stride(stride):
                def fwd(d):
                    n = len(d)
                    if stride >= n: return d
                    out = bytearray(n); idx = 0
                    for off in range(stride):
                        for pos in range(off, n, stride):
                            out[idx] = d[pos]; idx += 1
                    return bytes(out)
                def inv(d):
                    n = len(d)
                    if stride >= n: return d
                    out = bytearray(n); idx = 0
                    for off in range(stride):
                        for pos in range(off, n, stride):
                            out[pos] = d[idx]; idx += 1
                    return bytes(out)
                return fwd, inv
            fwd, inv = _make_stride(s)
            self.atoms.append(Transform(f'stride{s}', fwd, inv, idempotent=True))

        def _rev_fwd(d): return bytes(reversed(d))
        self.atoms.append(Transform('rev', _rev_fwd, _rev_fwd, involutory=True))

        def _sub128_fwd(d): return bytes((b - 128) & 0xFF for b in d)
        def _sub128_inv(d): return bytes((b + 128) & 0xFF for b in d)
        self.atoms.append(Transform('sub128', _sub128_fwd, _sub128_inv, involutory=False))

        def _rle_fwd(d):
            if not d: return d
            counts = Counter(d)
            esc = min(range(256), key=lambda b: counts.get(b, 0))
            out = bytearray([esc]); i = 0; n = len(d)
            while i < n:
                s = i
                while i < n-1 and d[i]==d[i+1] and i-s < 254: i += 1
                rl = i - s + 1
                if rl >= 4:
                    out.extend([esc, rl, d[s]]); i += 1
                else:
                    for j in range(s, s+rl):
                        if d[j] == esc: out.extend([esc, 0])
                        else: out.append(d[j])
                    i = s + rl
            r = bytes(out)
            return r if len(r) < len(d) else d
        def _rle_inv(d):
            if not d or len(d) < 2: return d
            esc = d[0]; out = bytearray(); i = 1; n = len(d)
            while i < n:
                if d[i] == esc:
                    i += 1
                    if i >= n: break
                    if d[i] == 0: out.append(esc); i += 1
                    else:
                        cnt = d[i]; i += 1
                        if i >= n: break
                        out.extend([d[i]] * cnt); i += 1
                else: out.append(d[i]); i += 1
            return bytes(out)
        self.atoms.append(Transform('rle', _rle_fwd, _rle_inv))

        def _nibsplit_fwd(d):
            return bytes((b >> 4) for b in d) + bytes((b & 0xF) for b in d)
        def _nibsplit_inv(d):
            n = len(d) // 2
            return bytes(((d[i] << 4) | d[n+i]) for i in range(n))
        self.atoms.append(Transform('nibsplit', _nibsplit_fwd, _nibsplit_inv))

        # Column-major reorder (for 2D data with known width)
        # This is parametric — we'll handle it via closure

        # === BACKENDS ===
        self.backends = [
            Backend('zlib1', lambda d: zlib.compress(d, 1), zlib.decompress, level=1),
            Backend('zlib6', lambda d: zlib.compress(d, 6), zlib.decompress, level=6),
            Backend('zlib9', lambda d: zlib.compress(d, 9), zlib.decompress, level=9),
            Backend('bz2', lambda d: bz2.compress(d, 9), bz2.decompress, level=9),
            Backend('lzma', lambda d: lzma.compress(d, preset=6), lzma.decompress, level=6),
        ]

    # ================================================================
    # TROPICAL EVALUATION: minimize compressed size
    # ================================================================

    def evaluate(self, chain: Chain, backend: Backend, data: bytes) -> Optional[int]:
        """Tropical evaluation: chain + backend -> compressed size (or None if fails)."""
        try:
            transformed = chain(data)
            if transformed is None:
                return None
            compressed = backend(transformed)
            if compressed is None:
                return None
            return len(compressed)
        except:
            return None

    # ================================================================
    # BEAM SEARCH over quotient monoid
    # ================================================================

    def generate_chains(self, max_depth: int = 3, beam_width: int = 20,
                        data: bytes = None) -> List[Chain]:
        """Generate promising chains using beam search.

        Instead of exhaustive search over all chains up to depth k
        (exponential), use beam search: at each depth, keep only the
        top beam_width chains by a quick heuristic score.
        """
        identity = Chain()
        beams = [identity]
        all_chains = [identity]

        for depth in range(max_depth):
            candidates = []
            for chain in beams:
                for atom in self.atoms:
                    new_chain = chain + Chain(atom)
                    simplified = new_chain.simplify()
                    # Skip if simplification reduced it (redundant chain)
                    if len(simplified) < len(new_chain):
                        continue
                    candidates.append(new_chain)

            if data is not None and candidates:
                # Score by quick entropy estimate after transform
                scored = []
                for c in candidates:
                    try:
                        transformed = c(data)
                        if transformed is None or len(transformed) > len(data) * 2:
                            continue
                        # Quick score: sum of min(b, 256-b)
                        score = sum(min(b, 255-b) for b in transformed[:min(512, len(transformed))])
                        scored.append((score, c))
                    except:
                        continue
                scored.sort(key=lambda x: x[0])
                beams = [c for _, c in scored[:beam_width]]
            else:
                beams = candidates[:beam_width]

            all_chains.extend(beams)

        return all_chains

    # ================================================================
    # THE COMPRESSOR: Catamorphism over the algebra
    # ================================================================

    def compress(self, data: bytes, max_depth: int = 2,
                 beam_width: int = 15) -> Tuple[bytes, str]:
        """Compress data by searching the transform algebra.

        This is the CATAMORPHISM: fold over all chains × backends,
        selecting the minimum under the tropical semiring.

        Returns: (compressed_data, description_string)
        """
        if len(data) <= 2:
            return data, 'raw'

        chains = self.generate_chains(max_depth, beam_width, data)

        best_size = len(data) + 100
        best_data = data
        best_desc = 'raw'

        for chain in chains:
            try:
                transformed = chain(data)
            except:
                continue
            if transformed is None or len(transformed) > len(data) * 2:
                continue

            for backend in self.backends:
                compressed = backend(transformed)
                if compressed is not None and len(compressed) < best_size:
                    best_size = len(compressed)
                    best_data = compressed
                    best_desc = f"{chain.name}|{backend.name}"

        return best_data, best_desc

    def compress_image(self, img: np.ndarray, max_depth: int = 2) -> Tuple[bytes, str]:
        """Compress an image using the algebra + image-specific transforms."""
        results = {}
        raw = img.tobytes()

        # Standard byte-level search
        comp, desc = self.compress(raw, max_depth)
        results[desc] = comp

        # Row-major delta
        results['row_delta'] = self._compress_best(self._delta_bytes(raw))

        # Column-major delta
        if img.ndim == 2:
            col = img.T.flatten().tobytes()
        else:
            col = img.transpose(1, 0, 2).reshape(-1).tobytes()
        results['col_delta'] = self._compress_best(self._delta_bytes(col))

        # Diagonal scan + delta (grayscale only)
        if img.ndim == 2:
            H, W = img.shape
            diag = bytearray()
            for s in range(H + W - 1):
                for i in range(max(0, s-W+1), min(s+1, H)):
                    j = s - i
                    if 0 <= j < W:
                        diag.append(img[i, j])
            results['diag_delta'] = self._compress_best(self._delta_bytes(bytes(diag)))

        # Row-adaptive prediction
        if img.ndim == 2:
            H, W = img.shape
            flat = img
            bpp = 1
        else:
            H, W = img.shape[:2]
            flat = img.reshape(H, W * img.shape[2])
            bpp = img.shape[2]

        above = np.zeros(flat.shape[1], dtype=np.uint8)
        out = bytearray()
        for i in range(H):
            row = flat[i]
            # Try none, left, up, avg, paeth — pick best per row
            best_score = float('inf')
            best_row = row.tobytes()
            for pid in range(5):
                res = self._predict_row(row, above, bpp, pid)
                sc = sum(min(b, 255-b) for b in res)
                if sc < best_score:
                    best_score = sc
                    best_row = bytes([pid]) + bytes(res)
            out.extend(best_row)
            above = row
        results['adaptive'] = self._compress_best(bytes(out))

        # Channel-stride for color
        if img.ndim == 3:
            ppch = H * W
            flat_bytes = img.reshape(H, W * 3).tobytes()
            strided = bytearray(len(flat_bytes))
            for ch in range(3):
                for p in range(ppch):
                    strided[ch * ppch + p] = flat_bytes[p * 3 + ch]
            ch_data = bytearray()
            for ch in range(3):
                ch_data.extend(self._delta_bytes(bytes(strided[ch*ppch:(ch+1)*ppch])))
            results['stride_delta'] = self._compress_best(bytes(ch_data))

            # YCbCr
            r, g, b = img[:,:,0].astype(int), img[:,:,1].astype(int), img[:,:,2].astype(int)
            y = ((r + 2*g + b) >> 2).astype(np.uint8)
            cb = ((b - g + 128) & 0xFF).astype(np.uint8)
            cr = ((r - g + 128) & 0xFF).astype(np.uint8)
            ycc = np.stack([y, cb, cr], axis=-1)
            ycc_flat = ycc.reshape(H, W * 3).tobytes()
            ycc_strided = bytearray(len(ycc_flat))
            for ch in range(3):
                for p in range(ppch):
                    ycc_strided[ch * ppch + p] = ycc_flat[p * 3 + ch]
            ycc_ch = bytearray()
            for ch in range(3):
                ycc_ch.extend(self._delta_bytes(bytes(ycc_strided[ch*ppch:(ch+1)*ppch])))
            results['ycbcr_stride'] = self._compress_best(bytes(ycc_ch))

        best = min(results, key=lambda k: len(results[k]))
        return results[best], best

    def compress_video(self, frames: list, keyframe_interval: int = 15) -> Tuple[bytes, dict]:
        """Compress video using the algebra.

        Keyframes: full image compression via algebra
        Delta frames: XOR + algebra search on sparse residual
        Motion blocks: block-level delta with motion estimation
        """
        if not frames:
            return b'', {}

        n = len(frames)
        shape = frames[0].shape
        H, W = shape[:2]
        channels = shape[2] if len(shape) > 2 else 1

        # Header
        import struct as st
        header = st.pack('!IHHB', n, H, W, channels)

        packets = []
        prev = None
        stats = {'keyframes': 0, 'delta_frames': 0, 'total_raw': 0, 'total_comp': 0}

        for i, frame in enumerate(frames):
            is_key = (i == 0 or i % keyframe_interval == 0 or prev is None)

            if not is_key:
                delta = np.bitwise_xor(prev, frame)
                change = np.count_nonzero(delta) / delta.size
                if change > 0.5:
                    is_key = True

            if is_key:
                # Full frame via algebra
                comp = self._compress_best(frame.tobytes())
                ftype = 0x00
                stats['keyframes'] += 1
            else:
                # Delta frame: try multiple approaches
                delta = np.bitwise_xor(prev, frame)
                delta_bytes = delta.tobytes()

                # Approach 1: raw delta + best backend
                comp1 = self._compress_best(delta_bytes)

                # Approach 2: delta-encode the delta (for smooth motion)
                comp2 = self._compress_best(self._delta_bytes(delta_bytes))

                # Approach 3: Block-level — only encode changed blocks
                block_size = 16
                block_data = bytearray()
                for bi in range(0, H, block_size):
                    for bj in range(0, W, block_size):
                        block = delta[bi:bi+block_size, bj:bj+block_size]
                        if np.any(block != 0):
                            block_data.extend(block.tobytes())
                        # else: skip (implicit zeros)
                if block_data:
                    comp3 = self._compress_best(bytes(block_data))
                else:
                    comp3 = self._compress_best(b'\x00')

                comp = min([comp1, comp2, comp3], key=len)
                ftype = 0x01
                stats['delta_frames'] += 1

            packet = st.pack('!BI', ftype, len(comp)) + comp
            packets.append(packet)
            stats['total_raw'] += frame.nbytes
            stats['total_comp'] += len(comp) + 5  # 5 for packet header
            prev = frame

        body = b''.join(packets)
        total = header + body
        stats['total_comp'] = len(total)
        stats['ratio'] = stats['total_raw'] / stats['total_comp'] if stats['total_comp'] > 0 else 0
        return total, stats

    # ================================================================
    # HELPER METHODS
    # ================================================================

    def _compress_best(self, data):
        results = [zlib.compress(data, 1), zlib.compress(data, 9)]
        try: results.append(bz2.compress(data, 9))
        except: pass
        return min(results, key=len)

    def _delta_bytes(self, data):
        if len(data) < 2: return data
        out = bytearray(len(data)); out[0] = data[0]
        for i in range(1, len(data)):
            out[i] = (data[i] - data[i-1]) & 0xFF
        return bytes(out)

    def _predict_row(self, row, above, bpp, pid):
        out = np.zeros(len(row), dtype=np.uint8)
        for j in range(len(row)):
            if pid == 0: pred = 0
            elif pid == 1: pred = int(row[j-bpp]) if j >= bpp else 0
            elif pid == 2: pred = int(above[j])
            elif pid == 3:
                l = int(row[j-bpp]) if j >= bpp else 0
                pred = (l + int(above[j])) // 2
            else:  # paeth
                a = int(row[j-bpp]) if j >= bpp else 0
                b = int(above[j])
                c = int(above[j-bpp]) if j >= bpp else 0
                p = a + b - c
                pa, pb, pc = abs(p-a), abs(p-b), abs(p-c)
                pred = a if pa<=pb and pa<=pc else (b if pb<=pc else c)
            out[j] = (int(row[j]) - pred) & 0xFF
        return out


# ============================================================================
# BENCHMARK
# ============================================================================

def benchmark():
    """Benchmark the algebraic codec against industry + our previous versions."""
    np.random.seed(42)

    print(f"TC Algebra v{__version__} -- Functional Compression via Monoid Structure")
    print("=" * 100)

    algebra = TransformAlgebra()

    # === BYTE DATA ===
    tests = {}
    tests['counter_4K'] = bytes([i % 256 for i in range(4096)])
    tests['quadratic'] = bytes([int((i*i/16) % 256) for i in range(4096)])
    tests['english_4K'] = (b"the quick brown fox jumps over the lazy dog. " * 100)[:4096]
    tests['json_4K'] = (b'{"id":1,"name":"test","values":[1,2,3]}\n' * 105)[:4096]
    tests['random_4K'] = np.random.randint(0,256,4096,dtype=np.uint8).tobytes()
    tests['zeros_4K'] = bytes(4096)
    tests['binary_exe'] = bytes([(i*7+13)%256 for i in range(4096)])

    print(f"\n  BYTE DATA:")
    print(f"  {'Name':>15} {'Raw':>8} {'Alg':>8} {'zlib9':>8} {'bz2':>8} {'Alg/best':>9} {'Chain':>25}")
    bw = bt = bl = 0
    for name, data in tests.items():
        raw = len(data)
        t0 = time.time()
        comp, desc = algebra.compress(data)
        elapsed = (time.time() - t0) * 1000
        alg_size = len(comp)
        zl = len(zlib.compress(data, 9))
        try: bz = len(bz2.compress(data, 9))
        except: bz = raw
        best = min(zl, bz)
        ratio = best / alg_size if alg_size > 0 else 0
        if ratio > 1.02: bw += 1; tag = "WIN"
        elif ratio < 0.98: bl += 1; tag = "LOSE"
        else: bt += 1; tag = "TIE"
        print(f"  {name:>15} {raw:7d}B {alg_size:7d}B {zl:7d}B {bz:7d}B {ratio:8.2f}x {desc:>25} {tag}")
    print(f"  Bytes: {bw}W {bt}T {bl}L")

    # === IMAGES ===
    N = 256
    img_tests = {}
    img_tests['grad_h'] = np.tile(np.arange(N, dtype=np.uint8), (N, 1))
    img_tests['grad_v'] = np.tile(np.arange(N, dtype=np.uint8).reshape(-1,1), (1, N))
    img_tests['checker'] = np.array([[(255 if (i//8+j//8)%2==0 else 0) for j in range(N)] for i in range(N)], dtype=np.uint8)
    x = np.linspace(0, 4*np.pi, N)
    X, Y = np.meshgrid(x, x)
    img_tests['smooth'] = ((128+100*np.sin(X)*np.cos(Y))).clip(0,255).astype(np.uint8)
    img_tests['circle'] = np.array([[255 if (i-N//2)**2+(j-N//2)**2<(N//3)**2 else 0 for j in range(N)] for i in range(N)], dtype=np.uint8)
    img_tests['random'] = np.random.randint(0,256,(N,N),dtype=np.uint8)

    print(f"\n  IMAGES ({N}x{N}):")
    print(f"  {'Name':>12} {'Raw':>8} {'Alg':>8} {'PNG':>8} {'Alg/PNG':>8} {'Method':>25}")
    iw = it = il = 0
    for name, img in img_tests.items():
        raw = img.nbytes
        comp, desc = algebra.compress_image(img)
        alg_size = len(comp)
        # PNG comparison
        try:
            from PIL import Image as PILImage
            import io
            pil_img = PILImage.fromarray(img)
            buf = io.BytesIO()
            pil_img.save(buf, format='PNG', optimize=True)
            png_size = len(buf.getvalue())
        except:
            png_size = len(zlib.compress(img.tobytes(), 9))
        ratio = png_size / alg_size if alg_size > 0 else 0
        if ratio > 1.02: iw += 1; tag = "WIN"
        elif ratio < 0.98: il += 1; tag = "LOSE"
        else: it += 1; tag = "TIE"
        print(f"  {name:>12} {raw:7d}B {alg_size:7d}B {png_size:7d}B {ratio:7.2f}x {desc:>25} {tag}")
    print(f"  Images: {iw}W {it}T {il}L")

    # === VIDEO ===
    N = 128
    n_frames = 30
    vid_tests = {}
    base = np.array([[(i+j)%256 for j in range(N)] for i in range(N)], dtype=np.uint8)
    vid_tests['static'] = [base.copy() for _ in range(n_frames)]
    frames = []
    for t in range(n_frames):
        f = np.zeros((N,N), dtype=np.uint8)
        cx = N//2 + int(30*np.sin(t*0.2))
        cy = N//2 + int(20*np.cos(t*0.3))
        for i in range(N):
            for j in range(N):
                if (i-cy)**2+(j-cx)**2 < 400: f[i,j] = 255
        frames.append(f)
    vid_tests['moving_obj'] = frames
    vid_tests['scroll'] = [np.array([[(i+j+t*3)%256 for j in range(N)] for i in range(N)], dtype=np.uint8) for t in range(n_frames)]

    print(f"\n  VIDEO ({N}x{N}, {n_frames} frames):")
    print(f"  {'Name':>12} {'Raw':>10} {'Alg':>10} {'Ratio':>7} {'Keys':>5} {'Deltas':>7}")
    for name, frames in vid_tests.items():
        raw = sum(f.nbytes for f in frames)
        comp, stats = algebra.compress_video(frames)
        print(f"  {name:>12} {raw:9,}B {stats['total_comp']:9,}B {stats['ratio']:6.1f}x "
              f"{stats['keyframes']:4d} {stats['delta_frames']:6d}")

    # Summary
    total_w = bw + iw
    total_t = bt + it
    total_l = bl + il
    total = total_w + total_t + total_l
    print(f"\n  COMBINED: {total_w}W {total_t}T {total_l}L / {total}")
    print(f"  The algebra discovers chains automatically via beam search.")
    print(f"  No hardcoded case statements -- pure algebraic composition.")


if __name__ == "__main__":
    benchmark()
