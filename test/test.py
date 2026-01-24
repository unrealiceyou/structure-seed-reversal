from math import gcd

class LCG:
    def __init__(self, seed: int):
        self.M = 25214903917
        self.A = 11
        self.modulus = 1 << 48
        self.state = (seed ^ self.M) % self.modulus
    
    def nextSeed(self) -> int:
        self.state = (self.state * self.M + self.A) % self.modulus
        return self.state
    
    def next(self, numBits: int) -> int:
        if not 1 <= numBits <= 48:
            raise ValueError("numBits must be between 1 and 48")
        self.nextSeed()
        return self.state >> (48 - numBits)
    
    def nextInt(self, bound: int) -> int:
        if bound <= 0:
            raise ValueError("bound must be positive")
        
        while True:
            bits = self.next(31)  
            value = bits % bound   
            
            if bits - value + (bound - 1) >= 0:
                return value
            print("Warning: nextInt regenerated a value!")
    
    def getCurrentSeed(self) -> int:
        return self.state


def region_lcg(w, x_r, z_r, salt):
    P = 341873128712
    Q = 132897987541

    return LCG(w + P*x_r + Q*z_r + salt)


if __name__ == "__main__":
    M = 25214903917
    A = 11
    P = 341873128712
    Q = 132897987541

    w = 9189798541153775729 % 2**48
    spacing = 34
    separation = 12
    bound = spacing - separation
    salt = 94251327

    # Generating the info tuple:
    # (x_r, z_r, x_c, z_c)
    lcgA = region_lcg(w, 1, 1, salt)
    infoA = (
        1,
        1,
        1*spacing + lcgA.nextInt(bound),
        1*spacing + lcgA.nextInt(bound),
    )

    lcgB = region_lcg(w, 2, 2, salt)
    infoB = (
        2,
        2,
        2*spacing + lcgB.nextInt(bound),
        2*spacing + lcgB.nextInt(bound),
    )


    lcgC = region_lcg(w, 3, 3, salt)
    infoC = (
        3,
        3,
        3*spacing + lcgC.nextInt(bound),
        3*spacing + lcgC.nextInt(bound),
    )

    lcgD = region_lcg(w, 64, 34, salt)
    infoD = (
        64,
        34,
        64*spacing + lcgD.nextInt(bound),
        34*spacing + lcgD.nextInt(bound),
    )

    print(f"Structure seed: {w}")
    info_list = [infoA, infoB, infoC, infoD]
    print(f"Constraints: {info_list}")

    # Testing w_h
    print("="*10 + f" w_h Check. w_h = {w >> 35} " + "="*10)

    x_r, z_r, x_c, z_c = infoA
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0 = 1
    else:
        c_0 = 0

    w_h = ((s_0 >> 35) - (D_h + c_0)) % 2**48
    print(f"infoA w_h = {w_h}, c_0 = {c_0}")

    x_r, z_r, x_c, z_c = infoB
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0 = 1
    else:
        c_0 = 0

    w_h = ((s_0 >> 35) - (D_h + c_0)) % 2**48
    print(f"infoB w_h = {w_h}, c_0 = {c_0}")

    x_r, z_r, x_c, z_c = infoC
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0 = 1
    else:
        c_0 = 0

    w_h = ((s_0 >> 35) - (D_h + c_0)) % 2**48
    print(f"infoC w_h = {w_h}, c_0 = {c_0}")


    x_r, z_r, x_c, z_c = infoD
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0 = 1
    else:
        c_0 = 0

    w_h = ((s_0 >> 35) - (D_h + c_0)) % 2**48
    print(f"infoD w_h = {w_h}, c_0 = {c_0}")


    # Testing the first e_0 constraint formula
    print("=============== Testing the first constraint formula (all of these should be 0)")
    x_r, z_r, x_c, z_c = infoA
    D = (P*x_r + Q*z_r + salt) % 2**48
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    X = ( dx - ((e_0*M + A) >> 17) ) % gcd(2**18, bound)
    print(f"infoA: {X}")

    x_r, z_r, x_c, z_c = infoB
    D = (P*x_r + Q*z_r + salt) % 2**48
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    X = ( dx - ((e_0*M + A) >> 17) ) % gcd(2**18, bound)
    print(f"infoB: {X}")

    x_r, z_r, x_c, z_c = infoC
    D = (P*x_r + Q*z_r + salt) % 2**48
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    X = ( dx - ((e_0*M + A) >> 17) ) % gcd(2**18, bound)
    print(f"infoC: {X}")

    x_r, z_r, x_c, z_c = infoD
    D = (P*x_r + Q*z_r + salt) % 2**48
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    X = ( dx - ((e_0*M + A) >> 17) ) % gcd(2**18, bound)
    print(f"infoD: {X}")

    # Testing the second e_0 constraint formula for all c_0 and c_1
    print("=============== Testing the second constraint formula (X should be equal for all constraints)")
    # Testing infoA
    x_r, z_r, x_c, z_c = infoA
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0_true = 1
    else:
        c_0_true = 0

    w_h = ((s_0 >> 35) - (D_h + c_0_true)) % 2**48

    A_true = (2**35 * M * w_h) % 2**48
    B_true = ((2**35 * (D_h + c_0_true) + e_0) * M + A) % 2**48

    if A_true + B_true >= 2**48:
        c_1_true = 1
    else:
        c_1_true = 0

    X_true = 2**18 * ((M * w_h) % 2**13) % bound

    print(f"infoA: c_0_true = {c_0_true}, c_1_true = {c_1_true}")
    print(f"X_true = {X_true}")
    for c_0 in range(2):
        for c_1 in range(2):
            B = ((2**35 * (D_h + c_0) + e_0) * M + A) % 2**48

            X = (dx - (B >> 17) + c_1*2**31) % bound
            print(f"c_0={c_0},\tc_1={c_1},\tX={X}")
            

    # Testing infoB
    x_r, z_r, x_c, z_c = infoB
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0_true = 1
    else:
        c_0_true = 0

    w_h = ((s_0 >> 35) - (D_h + c_0_true)) % 2**48

    A_true = (2**35 * M * w_h) % 2**48
    B_true = ((2**35 * (D_h + c_0_true) + e_0) * M + A) % 2**48

    if A_true + B_true >= 2**48:
        c_1_true = 1
    else:
        c_1_true = 0

    X_true = 2**18 * ((M * w_h) % 2**13) % bound

    print(f"infoB: c_0_true = {c_0_true}, c_1_true = {c_1_true}")
    print(f"X_true = {X_true}")
    for c_0 in range(2):
        for c_1 in range(2):
            B = ((2**35 * (D_h + c_0) + e_0) * M + A) % 2**48

            X = (dx - (B >> 17) + c_1*2**31) % bound
            print(f"c_0={c_0},\tc_1={c_1},\tX={X}")

    
    # Testing infoC
    x_r, z_r, x_c, z_c = infoC
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0_true = 1
    else:
        c_0_true = 0

    w_h = ((s_0 >> 35) - (D_h + c_0_true)) % 2**48

    A_true = (2**35 * M * w_h) % 2**48
    B_true = ((2**35 * (D_h + c_0_true) + e_0) * M + A) % 2**48

    if A_true + B_true >= 2**48:
        c_1_true = 1
    else:
        c_1_true = 0

    X_true = 2**18 * ((M * w_h) % 2**13) % bound

    print(f"infoC: c_0_true = {c_0_true}, c_1_true = {c_1_true}")
    print(f"X_true = {X_true}")
    for c_0 in range(2):
        for c_1 in range(2):
            B = ((2**35 * (D_h + c_0) + e_0) * M + A) % 2**48

            X = (dx - (B >> 17) + c_1*2**31) % bound
            print(f"c_0={c_0},\tc_1={c_1},\tX={X}")

    
    # Testing infoD
    x_r, z_r, x_c, z_c = infoD
    D = (P*x_r + Q*z_r + salt) % 2**48
    D_h = D >> 35
    
    s_0 = ((w + D) ^ M) % 2**48
    e_0 = s_0 % 2**35
    dx = x_c - x_r*spacing
    dz = z_c - z_r*spacing

    if w % 2**35 + D % 2**35 >= 2**35:
        c_0_true = 1
    else:
        c_0_true = 0

    w_h = ((s_0 >> 35) - (D_h + c_0_true)) % 2**48

    A_true = (2**35 * M * w_h) % 2**48
    B_true = ((2**35 * (D_h + c_0_true) + e_0) * M + A) % 2**48

    if A_true + B_true >= 2**48:
        c_1_true = 1
    else:
        c_1_true = 0

    X_true = 2**18 * ((M * w_h) % 2**13) % bound

    print(f"infoC: c_0_true = {c_0_true}, c_1_true = {c_1_true}")
    print(f"X_true = {X_true}")
    for c_0 in range(2):
        for c_1 in range(2):
            B = ((2**35 * (D_h + c_0) + e_0) * M + A) % 2**48

            X = (dx - (B >> 17) + c_1*2**31) % bound
            print(f"c_0={c_0},\tc_1={c_1},\tX={X}")


