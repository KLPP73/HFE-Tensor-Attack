load("aadp_construction.sage")


# Bytes per field element; BN254 fits in 32 bytes
ENTRY_BYTES = 32


def write_challenge(M_hat, target, output_path):
    # File layout:
    # - 32-byte target
    # - M_hat^(0), ..., M_hat^(n) row-major, each entry 32 bytes big-endian
    # - M_hat^(0) holds the encrypted message in its bottom-right entry
    # - M_hat^(j) for j >= 1 are unchanged from M^(j)
    with open(output_path, "wb") as f:
        f.write(int(target).to_bytes(ENTRY_BYTES, byteorder="big"))
        for Mj in M_hat:
            for r in range(Mj.nrows()):
                for c in range(Mj.ncols()):
                    f.write(int(Mj[r, c]).to_bytes(ENTRY_BYTES, byteorder="big"))


def read_challenge(input_path, m, F):
    # Returns (target, [M_hat^(0), ..., M_hat^(n)])
    n = m - 1
    k = 2 * m + 1
    with open(input_path, "rb") as f:
        data = f.read()
    expected = ENTRY_BYTES + (n + 1) * k * k * ENTRY_BYTES
    if len(data) != expected:
        raise ValueError(f"file size mismatch: got {len(data)}, expected {expected}")
    target = F(int.from_bytes(data[:ENTRY_BYTES], "big"))
    pos = ENTRY_BYTES
    M_hat = []
    for j in range(n + 1):
        Mj = matrix(F, k, k)
        for r in range(k):
            for c in range(k):
                Mj[r, c] = F(int.from_bytes(data[pos:pos + ENTRY_BYTES], "big"))
                pos += ENTRY_BYTES
        M_hat.append(Mj)
    return target, M_hat


def write_secrets(msg, params, output_path):
    # Layout:
    # - 32-byte msg
    # - for each gate i in 0..m-1: L_i row-major (k*4), R_i row-major (4*k), xi_i (n+1)
    # - each field element is 32 bytes big-endian
    L, R, xi = params["L"], params["R"], params["xi"]
    m = params["m"]
    with open(output_path, "wb") as f:
        f.write(int(msg).to_bytes(ENTRY_BYTES, byteorder="big"))
        for i in range(m):
            for r in range(L[i].nrows()):
                for c in range(L[i].ncols()):
                    f.write(int(L[i][r, c]).to_bytes(ENTRY_BYTES, byteorder="big"))
            for r in range(R[i].nrows()):
                for c in range(R[i].ncols()):
                    f.write(int(R[i][r, c]).to_bytes(ENTRY_BYTES, byteorder="big"))
            for j in range(len(xi[i])):
                f.write(int(xi[i][j]).to_bytes(ENTRY_BYTES, byteorder="big"))


def read_secrets(input_path, m, F):
    # Returns (msg, L, R, xi)
    n = m - 1
    k = 2 * m + 1
    expected = ENTRY_BYTES + m * (8 * k + n + 1) * ENTRY_BYTES
    with open(input_path, "rb") as f:
        data = f.read()
    if len(data) != expected:
        raise ValueError(f"secrets file size mismatch: got {len(data)}, expected {expected}")
    pos = 0
    msg = F(int.from_bytes(data[pos:pos + ENTRY_BYTES], "big"))
    pos += ENTRY_BYTES
    L_list, R_list, xi_list = [], [], []
    for i in range(m):
        Li = matrix(F, k, 4)
        for r in range(k):
            for c in range(4):
                Li[r, c] = F(int.from_bytes(data[pos:pos + ENTRY_BYTES], "big"))
                pos += ENTRY_BYTES
        L_list.append(Li)
        Ri = matrix(F, 4, k)
        for r in range(4):
            for c in range(k):
                Ri[r, c] = F(int.from_bytes(data[pos:pos + ENTRY_BYTES], "big"))
                pos += ENTRY_BYTES
        R_list.append(Ri)
        xi_i = []
        for j in range(n + 1):
            xi_i.append(F(int.from_bytes(data[pos:pos + ENTRY_BYTES], "big")))
            pos += ENTRY_BYTES
        xi_list.append(xi_i)
    return msg, L_list, R_list, xi_list


def roundtrip_test(params, msg, M_hat, public_path, secret_path):
    # read both files back and check:
    # - public: target and matrices match byte-for-byte
    # - secret: msg, L, R, xi match byte-for-byte
    # - reconstr: M_hat rebuilt from read-back (msg, L, R, xi) matches the original
    # - decrypt: with a valid witness, decrypt(read-back ciphertext) gives msg back
    F, m, n, k = params["F"], params["m"], params["n"], params["k"]
    A, B, C, D = params["A"], params["B"], params["C"], params["D"]

    target_read, M_hat_read = read_challenge(public_path, m, F)
    target_ok = (target_read == params["target"])
    matrices_ok = all(M_hat_read[j] == M_hat[j] for j in range(n + 1))
    print(f"  public  : target match {target_ok}, matrices match {matrices_ok}")

    msg_read, L_read, R_read, xi_read = read_secrets(secret_path, m, F)
    msg_ok = (msg_read == msg)
    L_ok = all(L_read[i] == params["L"][i] for i in range(m))
    R_ok = all(R_read[i] == params["R"][i] for i in range(m))
    xi_ok = all(xi_read[i] == params["xi"][i] for i in range(m))
    print(f"  secret  : msg match {msg_ok}, L match {L_ok}, R match {R_ok}, xi match {xi_ok}")

    M_rebuilt = compute_M(L_read, R_read, xi_read, A, B, C, D, m, n, k, F)
    M_hat_rebuilt = [matrix(F, Mj) for Mj in M_rebuilt]
    M_hat_rebuilt[0][k - 1, k - 1] += msg_read
    reconstr_ok = all(M_hat_rebuilt[j] == M_hat[j] for j in range(n + 1))
    print(f"  reconstr: M_hat rebuilt from read-back secrets matches original {reconstr_ok}")

    all_ok = target_ok and matrices_ok and msg_ok and L_ok and R_ok and xi_ok and reconstr_ok
    if not all_ok:
        return False
    if not VALID_WITNESS_MODE:
        print(f"  decrypt : skipped (no witness in non-residue mode)")
        return True
    recovered = decrypt(M_hat_read, params["witness"])
    ok = (recovered == msg)
    print(f"  decrypt : msg = {to_hex(msg)}, recovered = {to_hex(recovered)}, match {ok}")
    return ok


def main():
    if len(sys.argv) < 3:
        print("Usage: sage aadp_challenge_writer.sage <m> <seed> [output_path]")
        print("  Writes two files:")
        print("  <output_path>            public ciphertext (default aadp_challenge_m<m>.bin)")
        print("  <output_path>.secret.bin private sidecar with msg + (L, R, xi)")
        return
    m = int(sys.argv[1])
    seed = int(sys.argv[2])
    output_path = sys.argv[3] if len(sys.argv) >= 4 else f"aadp_challenge_m{m}.bin"
    if output_path.endswith(".bin"):
        secret_path = output_path[:-4] + ".secret.bin"
    else:
        secret_path = output_path + ".secret.bin"

    F = GF(DEFAULT_PRIME)
    params = setup(F, m, normalized=True, seed=seed)
    msg = F.random_element()
    M_hat = encrypt(params, msg)
    write_challenge(M_hat, params["target"], output_path)
    write_secrets(msg, params, secret_path)

    n, k = params["n"], params["k"]
    public_bytes = ENTRY_BYTES + (n + 1) * k * k * ENTRY_BYTES
    secret_bytes = ENTRY_BYTES + m * (8 * k + n + 1) * ENTRY_BYTES
    print(f"--- AADP challenge (m={m}, valid_witness={VALID_WITNESS_MODE}, seed={seed}) ---")
    print(f"  target      : {to_hex(params['target'])}")
    print(f"  msg         : {to_hex(msg)}    (decryption answer; also in secret file)")
    print(f"  public file : {output_path}  ({public_bytes} bytes)")
    print(f"  secret file : {secret_path}  ({secret_bytes} bytes)")
    roundtrip_test(params, msg, M_hat, output_path, secret_path)
    if VALID_WITNESS_MODE:
        print(f"NOTE: VALID_WITNESS_MODE=True, so target is a square")
        print(f"flip it in aadp_construction.sage for the actual challenge (non-square target)")


main()

# ~ if __name__ == "__main__":
    # ~ main()
