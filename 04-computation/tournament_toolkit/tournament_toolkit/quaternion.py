"""
quaternion.py — Quaternion linear algebra for attention heads.

Provides QuaternionLinear: a drop-in replacement for standard
linear layers that achieves 75% parameter savings by coupling
four components via the Hamilton product.

The Hamilton product rules (ij=k, jk=i, ki=j) encode the
tournament 3-cycle structure at the algebraic level.

Usage:
    from tournament_toolkit.quaternion import QuaternionLinear, QuaternionAttention

    # Drop-in linear replacement (75% fewer params)
    layer = QuaternionLinear(d_in=256, d_out=64)
    output = layer(input_array)  # input shape: (batch, 256)

    # Full attention head
    attn = QuaternionAttention(d_model=256, d_head=64)
    output = attn(sequence)  # shape: (seq_len, 256)

Based on:
    Tay et al., ACL 2019 (Quaternion Transformers)
    Grassucci et al., ICLR 2021 (Quaternion Neural Networks)
"""

import numpy as np


class QuaternionLinear:
    """Linear layer using Hamilton product coupling.

    Instead of a single d_in × d_out weight matrix (d_in * d_out params),
    uses four (d_in/4) × (d_out/4) component matrices coupled via
    the Hamilton product (d_in * d_out / 4 params = 75% savings).

    Args:
        d_in: Input dimension (must be divisible by 4)
        d_out: Output dimension (must be divisible by 4)
        scale: Initialization scale (default: 0.02)
    """

    def __init__(self, d_in, d_out, scale=0.02):
        if d_in % 4 != 0:
            raise ValueError(f"d_in must be divisible by 4, got {d_in}")
        if d_out % 4 != 0:
            raise ValueError(f"d_out must be divisible by 4, got {d_out}")

        self.d_in = d_in
        self.d_out = d_out
        self.d_q_in = d_in // 4
        self.d_q_out = d_out // 4

        # Four quaternion component weight matrices
        self.W_r = np.random.randn(self.d_q_in, self.d_q_out) * scale
        self.W_i = np.random.randn(self.d_q_in, self.d_q_out) * scale
        self.W_j = np.random.randn(self.d_q_in, self.d_q_out) * scale
        self.W_k = np.random.randn(self.d_q_in, self.d_q_out) * scale

    @property
    def num_params(self):
        return 4 * self.d_q_in * self.d_q_out

    @property
    def standard_params(self):
        return self.d_in * self.d_out

    @property
    def savings_fraction(self):
        return 1.0 - self.num_params / self.standard_params

    def forward(self, x):
        """Apply quaternion linear transformation.

        Args:
            x: Input array of shape (..., d_in)

        Returns:
            Output array of shape (..., d_out)
        """
        d = self.d_q_in
        r = x[..., :d]
        i = x[..., d:2*d]
        j = x[..., 2*d:3*d]
        k = x[..., 3*d:4*d]

        # Hamilton product: (W_r + W_i·i + W_j·j + W_k·k)(r + i·i + j·j + k·k)
        o_r = r @ self.W_r - i @ self.W_i - j @ self.W_j - k @ self.W_k
        o_i = r @ self.W_i + i @ self.W_r + j @ self.W_k - k @ self.W_j
        o_j = r @ self.W_j - i @ self.W_k + j @ self.W_r + k @ self.W_i
        o_k = r @ self.W_k + i @ self.W_j - j @ self.W_i + k @ self.W_r

        return np.concatenate([o_r, o_i, o_j, o_k], axis=-1)

    def __call__(self, x):
        return self.forward(x)

    def __repr__(self):
        return (f"QuaternionLinear({self.d_in}, {self.d_out}, "
                f"params={self.num_params} vs standard {self.standard_params}, "
                f"savings={self.savings_fraction:.0%})")


class QuaternionAttention:
    """Single attention head using quaternion weight coupling.

    Uses QuaternionLinear for Q, K, V projections and output projection.
    The Hamilton product couples all four components, saving 75% of
    per-head parameters compared to standard attention.

    Args:
        d_model: Model/embedding dimension (must be divisible by 4)
        d_head: Head dimension (must be divisible by 4)
    """

    def __init__(self, d_model, d_head, scale=0.02):
        self.d_model = d_model
        self.d_head = d_head

        self.W_Q = QuaternionLinear(d_model, d_head, scale)
        self.W_K = QuaternionLinear(d_model, d_head, scale)
        self.W_V = QuaternionLinear(d_model, d_head, scale)
        self.W_O = QuaternionLinear(d_head, d_model, scale)

    @property
    def num_params(self):
        return sum(w.num_params for w in [self.W_Q, self.W_K, self.W_V, self.W_O])

    @property
    def standard_params(self):
        return sum(w.standard_params for w in [self.W_Q, self.W_K, self.W_V, self.W_O])

    def forward(self, x):
        """Compute quaternion attention.

        Args:
            x: Input array of shape (seq_len, d_model)

        Returns:
            output: Array of shape (seq_len, d_model)
            attn_weights: Array of shape (seq_len, seq_len)
        """
        Q = self.W_Q(x)
        K = self.W_K(x)
        V = self.W_V(x)

        # Scaled dot-product attention
        scores = Q @ K.T / np.sqrt(self.d_head)

        # Softmax
        scores_exp = np.exp(scores - scores.max(axis=-1, keepdims=True))
        attn = scores_exp / scores_exp.sum(axis=-1, keepdims=True)

        # Apply attention and output projection
        context = attn @ V
        output = self.W_O(context)

        return output, attn

    def __call__(self, x):
        return self.forward(x)

    def __repr__(self):
        return (f"QuaternionAttention(d_model={self.d_model}, d_head={self.d_head}, "
                f"params={self.num_params} vs standard {self.standard_params}, "
                f"savings={1 - self.num_params/self.standard_params:.0%})")


class StandardAttention:
    """Standard attention head for comparison.

    Uses four independent weight matrices (no Hamilton coupling).
    """

    def __init__(self, d_model, d_head, scale=0.02):
        self.d_model = d_model
        self.d_head = d_head

        self.W_Q = np.random.randn(d_model, d_head) * scale
        self.W_K = np.random.randn(d_model, d_head) * scale
        self.W_V = np.random.randn(d_model, d_head) * scale
        self.W_O = np.random.randn(d_head, d_model) * scale

    @property
    def num_params(self):
        return sum(w.size for w in [self.W_Q, self.W_K, self.W_V, self.W_O])

    def forward(self, x):
        Q = x @ self.W_Q
        K = x @ self.W_K
        V = x @ self.W_V

        scores = Q @ K.T / np.sqrt(self.d_head)
        scores_exp = np.exp(scores - scores.max(axis=-1, keepdims=True))
        attn = scores_exp / scores_exp.sum(axis=-1, keepdims=True)

        context = attn @ V
        output = context @ self.W_O

        return output, attn

    def __call__(self, x):
        return self.forward(x)


def cartan_decompose(A):
    """Decompose matrix A into antisymmetric + symmetric traceless + scalar.

    Returns:
        A_anti: Antisymmetric part (tournament/directional)
        A_sym: Symmetric traceless part (cooperation/undirected)
        scalar: Scalar (uniform attention)
    """
    n = A.shape[0]
    A_anti = (A - A.T) / 2
    A_sym_full = (A + A.T) / 2
    scalar = np.trace(A_sym_full) / n
    A_sym = A_sym_full - scalar * np.eye(n)
    return A_anti, A_sym, scalar


def attention_to_tournament(attn):
    """Convert attention matrix to tournament adjacency.

    Arc i→j exists iff attention(i,j) > attention(j,i).

    Args:
        attn: (n, n) attention weight matrix

    Returns:
        adj: (n, n) binary tournament adjacency matrix
    """
    n = attn.shape[0]
    adj = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if attn[i, j] > attn[j, i]:
                adj[i, j] = 1
            else:
                adj[j, i] = 1
    return adj


def count_3cycles(adj):
    """Count directed 3-cycles in tournament adjacency matrix."""
    n = adj.shape[0]
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i,j] and adj[j,k] and adj[k,i]) or \
                   (adj[i,k] and adj[k,j] and adj[j,i]):
                    count += 1
    return count
