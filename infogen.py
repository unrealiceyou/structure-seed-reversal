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

def info_gen(w, x_r, z_r, salt):
    lcg = region_lcg(w, x_r, z_r, salt)
    info = (
        x_r,
        z_r,
        x_r*spacing + lcg.nextInt(bound),
        z_r*spacing + lcg.nextInt(bound),
    )
    return (info)


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
    infoList = [
        info_gen(w, 1, 1, salt),
        info_gen(w, 2, 1, salt),
        info_gen(w, 3, 1, salt),
        info_gen(w, 4, 1, salt),
        info_gen(w, 5, 1, salt),
        info_gen(w, 6, 1, salt),
        info_gen(w, 7, 1, salt),
        info_gen(w, 8, 1, salt),
        info_gen(w, 9, 1, salt),
        info_gen(w, 10, 1, salt),
        info_gen(w, 11, 1, salt),
    ]
    
    for info in infoList:
        print(info)
