# True:  target = v^{2^s} for random bits v, hence witness exists and decryption can be checked
# False: target is a random non-square, hence no witness (this is the challenge setup)
VALID_WITNESS_MODE = False

# rng seed
SEED = 42

# BN254 scalar field
# ~ DEFAULT_PRIME = 21888242871839275222246405745257275088548364400416034343698204186575808495617
DEFAULT_PRIME=65537

# Number of bitcheck gates
W = 3

def to_hex(x):
    return f"0x{int(x):x}"

def make_witness(F, m, bits):
    # Witness for the given bits, computed through the squaring chain
    s = m - W
    witness = [F(1)] + [F(b) for b in bits] + [F(0)] * (s - 1)
    cur = sum(F(b) * (1 << i) for i, b in enumerate(bits))
    for j in range(1, s):
        cur = cur ** 2
        witness[W + j] = cur
    return witness

def build_circuit(F, m):
    # Constraint matrices A, B, C, D such that gate i encodes A_i(x)*B_i(x) = C_i(x)*D_i(x)
    # Columns: 0 = x_0 (=1), 1..W = bits, W+j = squaring chain values
    # Gates: 0..W-1 = bitchecks, W..m-1 = squarings (last one is equality check)
    assert m >= 4
    s = m - W
    n = m - 1

    if VALID_WITNESS_MODE:
        bits = [randint(0, 1) for _ in range(W)]
        v = sum(bits[i] * (1 << i) for i in range(W))
        target = F(v) ** (2 ** s)
        witness = make_witness(F, m, bits)
    else:
        target = F.random_element()
        while target == 0 or target.is_square():
            target = F.random_element()
        witness = None

    A, B, C, D = (matrix(F, m, n + 1) for _ in range(4))

    # Bitchecks: a = x_0, b = c = d = b_i so that b_i^2 = x_0 * b_i
    for i in range(W):
        A[i, 0] = 1
        B[i, 1 + i] = 1
        C[i, 1 + i] = 1
        D[i, 1 + i] = 1

    # Squaring chain: gate W+j-1 enforces y_{j-1}^2 = y_j (or = target if j == s)
    # First squaring (j=1) takes the linear form sum 2^k * b_k as input
    for j in range(1, s + 1):
        i = W + j - 1
        D[i, 0] = 1
        if j == 1:
            for k in range(W):
                A[i, 1 + k] = 1 << k # coefficients 2^k
                B[i, 1 + k] = 1 << k
        else:
            A[i, W + j - 1] = 1
            B[i, W + j - 1] = 1
        if j == s:
            C[i, 0] = target # x_0 * target = target at last gate
        else:
            C[i, W + j] = 1

    return A, B, C, D, n, target, witness


def U_matrix(A, B, C, D, xi, i, j, F):
    # 4x4 U_i^{(j)} matrix
    return matrix(F, 4, 4, [
        [A[i, j], C[i, j], -xi[i][j], 0       ],
        [D[i, j], B[i, j],  0,        xi[i][j]],
        [0,       0,        B[i, j],  C[i, j] ],
        [0,       0,        D[i, j],  A[i, j] ],
    ])

def compute_M(L, R, xi, A, B, C, D, m, n, k, F):
    # Public matrices M^{(j)} = sum_i L_i U_i^{(j)} R_i for j = 0, ..., n
    M = [matrix(F, k, k) for _ in range(n + 1)]
    for j in range(n + 1):
        for i in range(m):
            M[j] += L[i] * U_matrix(A, B, C, D, xi, i, j, F) * R[i]
    return M

def setup(F, m, normalized=True, seed=None, r0=0, c0=0, r1=0, r1p=1, r_eq=2):
    # Build a fresh AADP instance
    # normalized=True puts the secret (L, R, xi) into canonical form via normalize()
    # per-gate randomness is reject-sampled so normalize() never bails on a precondition
    # pass seed for reproducible instance
    if seed is not None:
        set_random_seed(seed)
    A, B, C, D, n, target, witness = build_circuit(F, m)
    k = 2 * m + 1
    EI = E_indicators(F)

    def sample_gate(i):
        # reject anything that would later make normalize() blow up
        # u = L[r_eq, 0] / L[r_eq, 1] is invariant under the 6 standard steps
        # (cols 0 and 1 of L only get rescaled by alpha at step 6)
        # so checking 1 - target * u^2 != 0 here is the same check normalize() runs
        while True:
            xi_i = [F.random_element() for _ in range(n + 1)]
            L_i = random_matrix(F, k, 4)
            R_i = random_matrix(F, 4, k)
            if (L_i * EI["xi"] * R_i)[r0, c0] == 0:
                continue
            if L_i[r1, 0] == 0:
                continue
            if i < W and L_i[r1p, 3] == 0:
                continue
            if i == m - 1:
                if L_i[r_eq, 1] == 0:
                    continue
                u = L_i[r_eq, 0] / L_i[r_eq, 1]
                if 1 - target * u * u == 0:
                    continue
            return xi_i, L_i, R_i

    xi, L, R = [], [], []
    for i in range(m):
        xi_i, L_i, R_i = sample_gate(i)
        xi.append(xi_i)
        L.append(L_i)
        R.append(R_i)
    M = compute_M(L, R, xi, A, B, C, D, m, n, k, F)
    params = {
        "F": F, "m": m, "n": n, "k": k,
        "A": A, "B": B, "C": C, "D": D,
        "target": target, "witness": witness, "M": M,
        "xi": xi, "L": L, "R": R,
    }
    if normalized:
        params = normalize(params, r0=r0, c0=c0, r1=r1, r1p=r1p, r_eq=r_eq)
    return params


def encrypt(params, msg):
    # Add msg to the bottom-right entry of M^{(0)}
    F, k = params["F"], params["k"]
    M_hat = [matrix(F, Mj) for Mj in params["M"]]
    M_hat[0][k - 1, k - 1] += F(msg)
    return M_hat

def decrypt(M_hat, witness):
    # Delta is the rank-1 matrix with a 1 at the bottom-right
    # det(E - t * Delta) is linear in t
    F = M_hat[0].base_ring()
    k = M_hat[0].nrows()
    n = len(M_hat) - 1
    E = sum(witness[j] * M_hat[j] for j in range(n + 1))
    Rt = PolynomialRing(F, "t")
    t = Rt.gen()
    Delta = matrix(Rt, k, k)
    Delta[k - 1, k - 1] = Rt.one()
    det_poly = (E.change_ring(Rt) - t * Delta).det()
    if det_poly == 0:
        raise ValueError("det(E - t * Delta) is identically zero")
    roots = det_poly.roots()
    if len(roots) != 1:
        raise ValueError(f"expected one root, got {roots}")
    return roots[0][0]

# Normalization helpers
# Total freedoms: 6m + 4 (6 standard per gate + 1 extra at each of the 3 bitchecks and the equality)
def E_indicators(F):
    # Indicator matrices E^(X) such that U_i^{(j)} = sum_X X_{i,j} * E^(X) + xi_{i,j} * E^(xi)
    def E(r, c):
        M = matrix(F, 4, 4)
        M[r, c] = 1
        return M
    return {
        "A":  E(0, 0) + E(3, 3),
        "B":  E(1, 1) + E(2, 2),
        "C":  E(0, 1) + E(2, 3),
        "D":  E(1, 0) + E(3, 2),
        "xi": -E(0, 2) + E(1, 3),
    }


def P_bitcheck(F):
    # Matrix P for which P^2 = P for the bitcheck symmetry
    return matrix(F, 4, 4, [
        [0, 1, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, -1, 1],
    ])


def apply_standard_freedom(sym_id, i, param, L, R, xi, A, B, C, D, m, n, F):
    # Apply one of the 6 standard freedoms at gate i with the given parameter
    L = [matrix(F, Lj) for Lj in L]
    R = [matrix(F, Rj) for Rj in R]
    xi = [list(x) for x in xi]
    I4 = identity_matrix(F, 4)
    def E(r, c):
        M = matrix(F, 4, 4); M[r, c] = 1; return M

    if sym_id == 1:
        # Rescale: T^(xi) -> T^(xi) / c
        c = param
        L[i] = L[i] * matrix(F, 4, 4, [[1,0,0,0],[0,1,0,0],[0,0,c,0],[0,0,0,c]])
        R[i] = matrix(F, 4, 4, [[1,0,0,0],[0,1,0,0],[0,0,1/c,0],[0,0,0,1/c]]) * R[i]
        xi[i] = [c * xi[i][j] for j in range(n + 1)]
    elif sym_id == 2:
        # Shift A: T^(A) -> T^(A) + lam * T^(xi)
        lam = param
        L[i] = L[i] * (I4 + lam * E(1, 3))
        R[i] = (I4 - lam * E(0, 2)) * R[i]
        xi[i] = [xi[i][j] - lam * A[i, j] for j in range(n + 1)]
    elif sym_id == 3:
        # Shift B
        lam = param
        L[i] = L[i] * (I4 - lam * E(0, 2))
        R[i] = (I4 + lam * E(1, 3)) * R[i]
        xi[i] = [xi[i][j] - lam * B[i, j] for j in range(n + 1)]
    elif sym_id == 4:
        # Shift C
        lam = param
        L[i] = L[i] * (I4 + lam * E(1, 2))
        R[i] = (I4 - lam * E(1, 2)) * R[i]
        xi[i] = [xi[i][j] - lam * C[i, j] for j in range(n + 1)]
    elif sym_id == 5:
        # Shift D
        lam = param
        L[i] = L[i] * (I4 + lam * E(0, 3))
        R[i] = (I4 - lam * E(0, 3)) * R[i]
        xi[i] = [xi[i][j] + lam * D[i, j] for j in range(n + 1)]
    elif sym_id == 6:
        # alpha I: rescale L by alpha, R by 1/alpha
        alpha = param
        L[i] = alpha * L[i]
        R[i] = (1 / alpha) * R[i]
    else:
        raise ValueError(f"unknown sym_id {sym_id}")
    return L, R, xi


def apply_bitcheck_fix(i, t, L, R, F):
    # Apply fix for bitcheck symmetry at gate i with parameter t
    P = P_bitcheck(F)
    I4 = identity_matrix(F, 4)
    ML = I4 + (t - 1) * P
    MR = I4 + (1 / t - 1) * P.transpose()
    L = list(L); R = list(R)
    L[i] = L[i] * ML
    R[i] = MR * R[i]
    return L, R


def apply_equality_fix(i, u, target, L, R, F):
    # Apply fix for equality symmetry at gate i with parameter u
    # M_R(u) = M_L(u)^{-1}; singular when 1 - target * u^2 = 0
    D = 1 - target * u * u
    if D == 0:
        raise ValueError(f"singular parameter u^2 = 1/target")
    T = target
    ML = matrix(F, 4, 4, [
        [1,  -T*u, 0,    0  ],
        [-u, 1,    0,    0  ],
        [0,  0,    1,    T*u],
        [0,  0,    u,    1  ],
    ])
    MR = matrix(F, 4, 4, [
        [1/D,  T*u/D,  0,      0     ],
        [u/D,  1/D,    0,      0     ],
        [0,    0,      1/D,    -T*u/D],
        [0,    0,      -u/D,   1/D   ],
    ])
    L = list(L); R = list(R)
    L[i] = L[i] * ML
    R[i] = MR * R[i]
    return L, R


def apply_standard_norm_step(i, L, R, xi, A, B, C, D, m, n, k, F, r0, c0, r1, EI):
    # Apply steps 1-6 sequentially at gate i to enforce the 6 standard conditions
    T_xi = L[i] * EI["xi"] * R[i]
    if T_xi[r0, c0] == 0:
        raise ValueError(f"T_{i}^(xi)[{r0}, {c0}] = 0; pick different (r0, c0)")
    L, R, xi = apply_standard_freedom(1, i, T_xi[r0, c0], L, R, xi, A, B, C, D, m, n, F)
    T_A = L[i] * EI["A"] * R[i]
    L, R, xi = apply_standard_freedom(2, i, -T_A[r0, c0], L, R, xi, A, B, C, D, m, n, F)
    T_B = L[i] * EI["B"] * R[i]
    L, R, xi = apply_standard_freedom(3, i, -T_B[r0, c0], L, R, xi, A, B, C, D, m, n, F)
    T_C = L[i] * EI["C"] * R[i]
    L, R, xi = apply_standard_freedom(4, i, -T_C[r0, c0], L, R, xi, A, B, C, D, m, n, F)
    T_D = L[i] * EI["D"] * R[i]
    L, R, xi = apply_standard_freedom(5, i, T_D[r0, c0], L, R, xi, A, B, C, D, m, n, F)
    if L[i][r1, 0] == 0:
        raise ValueError(f"L_{i}[{r1}, 0] = 0; pick different r1")
    L, R, xi = apply_standard_freedom(6, i, 1 / L[i][r1, 0], L, R, xi, A, B, C, D, m, n, F)
    return L, R, xi

def normalize(params, r0=0, c0=0, r1=0, r1p=1, r_eq=2):
    # apply 6m standard + 3 bitcheck-extra + 1 equality-extra conditions:
    # - standard 6m, per gate: T^(xi)[r0,c0]=1, T^(A,B,C,D)[r0,c0]=0, L[r1,0]=1
    # - bitcheck-extra, 3 gates: L[r1p, 3] = 1
    # - equality-extra, 1 gate: L[r_eq, 0] = 0
    # all hold at once when we return
    F = params["F"]
    m, n, k = params["m"], params["n"], params["k"]
    A, B, C, D = params["A"], params["B"], params["C"], params["D"]
    target = params["target"]
    L = [matrix(F, Li) for Li in params["L"]]
    R = [matrix(F, Ri) for Ri in params["R"]]
    xi = [list(x) for x in params["xi"]]
    EI = E_indicators(F)

    M_orig = compute_M(L, R, xi, A, B, C, D, m, n, k, F)

    # Standard 6m conditions, gate by gate
    for i in range(m):
        L, R, xi = apply_standard_norm_step(i, L, R, xi, A, B, C, D, m, n, k, F, r0, c0, r1, EI)

    # Bitcheck-extra at the 3 bitcheck gates: scale so L[r1p, 3] = 1
    for i in range(W):
        if L[i][r1p, 3] == 0:
            raise ValueError(f"L_{i}[{r1p}, 3] = 0; pick different r1p")
        t = 1 / L[i][r1p, 3]
        L, R = apply_bitcheck_fix(i, t, L, R, F)

    # equality-extra at gate m-1: pick u so L[r_eq, 0] = 0 at the end
    # equality-extra changes L[:, 0] linearly in u
    eq = m - 1
    if L[eq][r_eq, 1] == 0:
        raise ValueError(f"L_{eq}[{r_eq}, 1] = 0; pick different r_eq")
    u = L[eq][r_eq, 0] / L[eq][r_eq, 1]
    L, R = apply_equality_fix(eq, u, target, L, R, F)

    # Re-run standard-6 at gate m-1 since equality symmetry fix destroyed them
    L, R, xi = apply_standard_norm_step(eq, L, R, xi, A, B, C, D, m, n, k, F, r0, c0, r1, EI)

    # Sanity check: M must match exactly
    M_final = compute_M(L, R, xi, A, B, C, D, m, n, k, F)
    if not all(M_final[j] == M_orig[j] for j in range(n + 1)):
        raise RuntimeError("normalization changed M (bug)")

    out = {**params}
    out["L"] = L
    out["R"] = R
    out["xi"] = xi
    out["M"] = M_final
    return out

def main():
    if len(sys.argv) < 2:
        print(f"Usage: sage aadp_construction.sage <m>")
        print(f"  m: number of gates (>= 4); 3 bitchecks + (m - 3) squaring gates")
        return
    m = int(sys.argv[1])
    F = GF(DEFAULT_PRIME)
    print(f"--- AADP construction (m={m}, valid_witness={VALID_WITNESS_MODE}, seed={SEED}) ---")
    params = setup(F, m, seed=SEED)
    msg = F.random_element()
    M_hat = encrypt(params, msg)
    print(f"  target    : {to_hex(params['target'])}")
    print(f"  msg       : {to_hex(msg)}")
    if VALID_WITNESS_MODE:
        recovered = decrypt(M_hat, params["witness"])
        print(f"  recovered : {to_hex(recovered)}")
        print(f"  decrypt OK: {msg == recovered}")
    else:
        print(f"  (no witness, decryption skipped)")


if __name__ == "__main__":
    main()
