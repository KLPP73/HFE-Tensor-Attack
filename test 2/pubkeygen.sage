import time
import random
import sys

q = 3

n = 12
m = 2*n-1 #must be <= 2n-1
argc = len(sys.argv)
if argc > 2:
    n = int(sys.argv[1])
    m = int(sys.argv[2])
    if m > 2*n-1:
        raise ValueError("m cannot be greater than 2n-1")

print(f"[*] Parameters :\n\t- q = {q}\n\t- n = {n}\n\t- m = {m}")


C2_list = Combinations(n, 2).list()
C3_list = Combinations(n, 3).list()
idx_C3 = {tuple(c): i for i, c in enumerate(C3_list)}

K = GF(q,repr="int")
F = K^n
f = F.basis() #f[0] to f[n-1] are the canonical basis vectors of F

E = ExteriorAlgebra(F)
e = E.basis(1).values() #e[0] to e[n-1] are the canonical basis vectors of F within E

#note: elements of F are coerced to E whenever added/multiplied with elements of E
#e.g. f[0]*f[1] = 0 is a dot product, but f[0]*e[1] = e[0]*e[1] is a wedge product

#option to fix random seed for reproducibility
#import datetime as dt
#s = int(dt.datetime.now().timestamp())
s = 11
set_random_seed(s)

# --- Matrix helper functions

def randomMatrix(a,b):
    #returns random matrix of maximal rank in K^(a*b)
    while True:
        M = matrix.random(K,a,b)
        if M.rank() == min(a,b):
            break
    return M
    
def randomSymmatrix():
    #returns random symmetric matrix of maximal rank in K^(n*n)
    while True:
        M = matrix.random(K,n)
        M = M + M.transpose()
        for i in range(n):
            M[i,i] = K.random_element() #in case we are in characteristic 2
        if M.is_invertible():
            break
    return M
    
def getHankelBasis():
    #returns a list of 2*n-1 matrices forming a basis of n-dimensional Hankel matrices
    return [matrix(K,n,lambda i,j:i+j==k) for k in range(2*n-1)]
    
def randomHankel():
    #returns random Hankel matrix in K^(n*n)
    return sum([H*K.random_element() for H in getHankelBasis()])
    
def isHankel(M):
    #returns True iff M is Hankel
    for i in range(n):
        for j in range(i+1):
            if (M[i-j,j] != M[i,0]) or (M[n-1-(i-j),n-1-j] != M[n-1-i,n-1]):
                return False
    return True
    
def flatten(M):
    #flattens a list of lists [[a,b],[c],..] into a single list [a,b,c,...]
    #if the input is a matrix, this turns it into the concatenation of its rows
    return [entry for row in M for entry in row]

def unflatten(v,a,b):
    #maps a vector in K^(a*b) to a matrix in K^(a x b), inverting flatten
    return matrix([v[b*i:b*(i+1)] for i in range(a)])
    
# --- Exterior algebra helper functions

def vectorToEk(k,w):
    #given a vector w of dimension binom(n,k), returns the corresponding degree-k element in E
    C = Combinations(n,k)
    L = C.list()
    x = E.zero()
    for i in range(C.cardinality()):
        x += w[i]*prod(e[j] for j in L[i])
    return x

def EkToVector(k,x):
    #inverse map of the previous one: given a degree-k element of E, returns its vector
    C = Combinations(n,k)
    w = [0]*C.cardinality()
    for (b,val) in x.list():
        i = C.rank(list(b))
        w[i] = val
    return vector(w)
    
def hodge(x):
    #given an element x of E returns its Hodge dual
    y = E(0)
    for (b,val) in x.list():
        y += val*prod(e[i] for i in b.complement())*prod(e[i] for i in range(b.capacity(),n))
    return y

def contract(x,y):
    #returns x contracted by y
    return hodge(hodge(x)*y)

I = matrix.identity(K,n)
J = matrix(K,n,lambda i,j:i+j==n-1)
    
# --- X1 key generation

def convolution(x,y):
    #returns convolution of vectors x and y of dimension n
    return vector([sum([x[j]*y[i-j] for j in range(max(0,i+1-n),min(n,i+1))]) for i in range(2*n-1)])
    

#Public key
def computePK(S,T):
    #returns list of public key matrices defining the m symmetric bilinear maps
    PK = [matrix(K,n) for i in range(m)]
    for i in range(n):
        for j in range(i+1):
            z = T*convolution(S*f[i],S*f[j])
            for k in range(m):
                PK[k][i,j] = PK[k][j,i] = z[k]
    return PK


# --- attack implementation

def phi_bil(A,x):
    #takes as input A in K^(n*n) and (x_1,...,x_n) in (Wedge^2(F))^n
    #outputs sum_i x_i*a_i where a_i is the i-th row of A
    return sum([x[i]*A[i] for i in range(n)])

def phi(PK,x):
    return [phi_bil(PK[k],x) for k in range(m)]

def phi_matrix(PK):
    m = len(PK)
    d2 = len(C2_list)
    d3 = len(C3_list)
    entries = {}
    for i in range(n):
        for idx_c2, (c1, c2) in enumerate(C2_list):
            row_idx = i * d2 + idx_c2
            for k in range(m):
                A = PK[k]
                col_offset = k * d3
                for j in range(n):
                    if j == c1 or j == c2:
                        continue         
                    val = A[i, j]
                    if val == 0:
                        continue
                    if j < c1:
                        c_sorted = (j, c1, c2)
                        sign = 1
                    elif j < c2:
                        c_sorted = (c1, j, c2)
                        sign = -1
                    else:
                        c_sorted = (c1, c2, j)
                        sign = 1              
                    col_idx = col_offset + idx_C3[c_sorted]
                    entries[(row_idx, col_idx)] = sign * val
    return matrix(K, n * d2, m * d3, entries, sparse=False)

def BtoE(b):
    #takes as input a vector in K^(binomial(n,2)*n) (e.g. an element of B)
    #outputs the corresponding list of n elements of Wedge^2(F)
    d = binomial(n,2)
    return [vectorToEk(2,b[i*d:(i+1)*d]) for i in range(n)]

def EtoB(L):
    #inverse map of BtoE
    return vector(flatten([EkToVector(2,L[i]) for i in range(n)]))
    
def phip(B,x):
    #second linear system from the overleaf writeup
    return [phi_bil(x,BtoE(b)) for b in B]

def phip_matrix(B):
    d2 = len(C2_list)
    d3 = len(C3_list)
    entries = {}
    for i in range(n):
        for j in range(n):
            row_idx = i * n + j
            for idx_B, b in enumerate(B):
                col_offset = idx_B * d3
                # Le bloc i du vecteur b contient les coefficients de e_{c1} ^ e_{c2}
                for idx_c2, (c1, c2) in enumerate(C2_list):
                    if j == c1 or j == c2:
                        continue    
                    val = b[i * d2 + idx_c2]
                    if val == 0:
                        continue
                    if j < c1:
                        c_sorted = (j, c1, c2)
                        sign = 1
                    elif j < c2:
                        c_sorted = (c1, j, c2)
                        sign = -1
                    else:
                        c_sorted = (c1, c2, j)
                        sign = 1
                    col_idx = col_offset + idx_C3[c_sorted]
                    entries[(row_idx, col_idx)] = sign * val
    return matrix(K, n * n, len(B) * d3, entries, sparse=False)

# -- Function to expand the public key using the distinguisher
def expand_public_key(PK):
    print(f"[*] Computing kernel of phi")
    Kerphi = phi_matrix(PK).kernel() # -- Computation of the orthogonal
    B = Kerphi.basis()
    print(f"[*] Computing kernel of phi_matrix")
    Kerphip = phip_matrix(B).kernel() # -- Computation of the biorthogonal
    print(f"[*] Kernel found of dimension {Kerphip.dimension()}")
    #Bp = Kerphip.basis()
    Bp = []
    # Getting a random basis 
    while len(Bp) < 2*n-1:
        cand = Kerphip.random_element()
        if cand not in span(Bp, K):
            Bp.append(cand)


    PKp = [unflatten(P, n, n) for P in Bp] 
    return PKp


def sidelnikov_shestakov(G, max_retries=50):
    """
    Recovers an equivalent private key (evaluation points and partial multipliers)
    from a GRS generator matrix G using the Sidelnikov-Shestakov cross-ratio method.
    """
    F = G.base_ring()
    k = G.nrows()
    n = G.ncols()
    r = n - k  
    
    for attempt in range(max_retries):
        # 1. Random permutation of the columns
        P = list(range(n))
        random.shuffle(P)
        G_perm = G.matrix_from_columns(P)
        
        # 2. Bring G_perm into systematic form [I_k | Y]
        try:
            G_sys = G_perm.echelon_form()
        except ValueError:
            continue
            
        Y = G_sys[:, k:]
        
        if 0 in Y.list():
            continue
            
        # 3. Fix 3 parameters to anchor the affine transformation.
        gamma = F(0)
        for elem in F:
            if elem not in [F(0), F(1)]:
                gamma = elem
                break
                
        x = [F(0)] * k
        y = [F(0)] * r
        x[0] = F(0)
        x[1] = F(1)
        y[0] = gamma
        
        C = (F(1) - gamma) / gamma
        
        degenerate_case = False
        
        for j in range(1, r):
            Z_1j = (Y[1, j] * Y[0, 0]) / (Y[1, 0] * Y[0, j])
            if Z_1j + C == F(0):
                degenerate_case = True
                break
            y[j] = Z_1j / (Z_1j + C)
            
        if degenerate_case: continue
            
        D = y[1] / gamma
        for i in range(2, k):
            Z_i1 = (Y[i, 1] * Y[0, 0]) / (Y[i, 0] * Y[0, 1])
            if Z_i1 - D == F(0):  # Protection ajoutée ici aussi !
                degenerate_case = True
                break
            x[i] = (Z_i1 * y[1] - D * gamma) / (Z_i1 - D)
            
        if degenerate_case: continue
            
        c = [F(0)] * k
        d = [F(0)] * r
        
        for i in range(k):
            c[i] = Y[i, 0] * (x[i] - y[0])
            
        for j in range(r):
            d[j] = (Y[0, j] * (x[0] - y[j])) / c[0]
            
        evaluation_points_perm = x + y
        multipliers_perm = c + d
        
        evaluation_points = [None] * n
        multipliers = [None] * n
        
        for new_idx, orig_idx in enumerate(P):
            evaluation_points[orig_idx] = evaluation_points_perm[new_idx]
            multipliers[orig_idx] = multipliers_perm[new_idx]
            
        c_orig = multipliers[:k]
        d_orig = multipliers[k:]
        
        return evaluation_points, c_orig, d_orig

    raise ValueError("Failure: impossible to find a valid configuration after {} trials.".format(max_retries))

# -- The attack function on the extended public key.
# -- Necessitates first to expand the public key to size 2n-1
def attack_public_key(PK):
    assert len(PK) == 2*n-1
    
    load("fast_rank1.spyx")  
    print(f"[*] Looking for rank 1 matrices")

    target_count = min(n+2,q+1)

    rank1_matrices = find_rank1_matrices(PK, K, n, 2*n-1, target_count)

    print(f"[*] Rank 1 matrices found !")

    if q > n:
        pts = []
        for mat in rank1_matrices:
            pivot_idx = -1
            for j in range(n):
                if mat[j,j] != 0:
                    pivot_idx = j
                    break
            if pivot_idx == -1:
                print(f"[!] Error : Matrix {idx} as a completely zero diagonal")
                continue
            w_j = mat[pivot_idx, pivot_idx].sqrt()
            col_j = mat.column(pivot_idx)

            w = col_j / w_j
            pts.append(w)


        print(f"[*] Succesfully extracted {len(pts)} points. Constructing generator matrix")

        G = column_matrix(pts) 
        
        ev_points, c, d = sidelnikov_shestakov(G)

        recovered_pairs = [(pts[i], ev_points[i]) for i in range(len(pts))]
        
        basis_pairs = recovered_pairs[:n]
        target_pair = recovered_pairs[n]

        w_target = target_pair[0]
        r_target = target_pair[1]

    # Build the Vandermonde vector for the target root
        v_target = vector(K, [r_target**i for i in range(n)])

    # Build the W and V matrices from the basis pairs
        W_basis = column_matrix([p[0] for p in basis_pairs])
        V_basis = column_matrix([vector(K, [p[1]**i for i in range(n)]) for p in basis_pairs])

        print("[*] Calculating linear dependencies for the (n+1)-th vector...")

    # 2. Find the linear combination coefficients in both spaces
    # W_basis * c = w_target  =>  c = W_basis^-1 * w_target
        c_coeffs = W_basis.inverse() * w_target

    # V_basis * d = v_target  =>  d = V_basis^-1 * v_target
        d_coeffs = V_basis.inverse() * v_target

    # 3. Resolve the scaling factors and construct the true linear W matrix
        W_scaled_cols = []
        for i in range(n):
            # The projective scaling factor gamma_i = d_i / c_i
            gamma_i = d_coeffs[i] / c_coeffs[i]
            
            # Scale the extracted vector to exactly match the linear transformation
            w_true = basis_pairs[i][0] / gamma_i
            W_scaled_cols.append(w_true)

        W_scaled = column_matrix(W_scaled_cols)
        print("[*] Projective scaling resolved. Matrices strictly aligned.")

        # 4. Recover S'
        # W_scaled = S'^T * V_basis  =>  S'^T = W_scaled * V_basis^-1
        S_prime_T = W_scaled * V_basis.inverse()
        S_prime = S_prime_T.transpose()

        print("\n[+] Equivalent Secret Matrix S' successfully computed!")

        return S_prime
    
    if q <= n:
        # On teste toutes les matrices pour voir laquelle donne un vecteur proportionnel à S^T e1
        for rank1_mat in rank1_matrices:
            # On décompose rank1_mat sous la forme c*v*v^t
            c = 0
            v = vector(rank1_mat.base_ring(), rank1_mat.nrows()) # Vecteur nul par défaut
            
            for i in range(rank1_mat.nrows()):
                if rank1_mat[i, i] != 0:
                    v = rank1_mat.column(i)
                    c = 1 / rank1_mat[i, i]
                    break
            
            if c != 0:
                k = len(PK)
                F = PK[0].base_ring()
                
                contraintes_M = {(0,1): 1}

                colonnes_systeme = []
                
                for r in range(n):
                    for s in range(n):
                        M_base = matrix(F, n, n)
                        M_base[r, s] = 1
                        colonne = []
                        for A_i in PK:
                            res = A_i * M_base.transpose() - M_base * A_i
                            colonne.extend(res.list()) 
                        colonnes_systeme.append(colonne)
                        
                # (Impact des y_i)
                for idx in range(k):
                    for j in range(n):
                        y_base = vector(F, n)
                        y_base[j] = 1
                        colonne = []
                        for i in range(k):
                            if i == idx:
                                terme = - (v.column() * y_base.row() - y_base.column() * v.row())
                                colonne.extend(terme.list())
                            else:
                                colonne.extend([0] * (n**2))
                        colonnes_systeme.append(colonne)
                        
                Matrice_Systeme = matrix(F, colonnes_systeme).transpose()
                
                lignes_systeme = Matrice_Systeme.rows()
                B_list = [0] * len(lignes_systeme) 
                
                #On ajoute les contraintes
                for (r, s), valeur in contraintes_M.items():
                    nouvelle_ligne = [0] * Matrice_Systeme.ncols()
                    index_inconnue = r * n + s 
                    nouvelle_ligne[index_inconnue] = 1 
                    
                    lignes_systeme.append(vector(F, nouvelle_ligne))
                    B_list.append(valeur)
                    
                Matrice_Systeme_Etendue = matrix(F, lignes_systeme)
                B = vector(F, B_list)
                Matrice_Systeme_Etendue = matrix(F, lignes_systeme)
        
                dim_solutions = Matrice_Systeme_Etendue.right_kernel().dimension()
                
                # 4. Résolution du système AX = B
                try:
                    vecteur_solution = Matrice_Systeme_Etendue.solve_right(B)
                    
                    print(f"[*] System solved with constraints!")
                    print(f"[*] Solution space dimension (remaining degrees of freedom)  : {dim_solutions}")
                    
                    if dim_solutions == 0:
                        print("[*] The solution found is unique.")
                    else:
                        print("[*] There are multiple solutions.")
                    
                    # Reconstruction de M
                    M_sol = matrix(F, n, n, vecteur_solution[:n**2].list())
                    
                    if M_sol.rank() >= n-1:
                        # Reconstruction des y_i
                        y_sols = []
                        offset = n**2
                        for i in range(k):
                            y_i_sol = vector(F, vecteur_solution[offset : offset + n].list())
                            y_sols.append(y_i_sol)
                            offset += n

                        S_T_cols = [v]

                        print("M sol rank: ", M_sol.rank())
                        char_poly = M_sol.characteristic_polynomial()
                        coeffs = char_poly.list()
                        v_noise = [-coeffs[n-1-k] for k in range(n)]

                        for k in range(n-1):
                            next_col = M_sol * S_T_cols[-1] - v_noise[k]*v
                            S_T_cols.append(next_col)
                        
                        S_prime_T = column_matrix(K, S_T_cols)
                        S_prime = S_prime_T.transpose()
                        
                        print("S prime rank: ", S_prime.rank())

                        return S_prime


                except ValueError:
                    print(f"[!] The system is incompatible with the constraints (impossible to find a solution)")
                    print(f"[I] The dimension of the kernel matrix with constraints was : {dim_solutions}")

     
    
def full_attack(PK):
    """
    The complete attack
    """
    if len(PK) > 4:
        PK_short = PK[:4]
        m = 4
    PKe = expand_public_key(PK_short)
    return attack_public_key(PKe)


S = randomMatrix(n,n) #this is S_1 in the article
T = randomMatrix(m,2*n-1)



PK = computePK(S,T)

# ~ RK = [randomSymmatrix() for i in range(m)]

# ~ print(PK)

save(PK,"public_key")





