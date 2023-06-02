#from https://github.com/uclcrypto/many01proofs/blob/main/montgomery.py
from gmpy2 import mpz, powmod, xmpz
from math import log, ceil

class Montgomery:
    def __init__(self, p, r_bits):
        self.p = p
        self.r = mpz(2**r_bits)
        self.r_bits = r_bits 
        self.r_mod_mask = (1 << r_bits)-1
        self.r_low_mask = (1 << (r_bits//2))-1

        assert(self.r > self.p)

        self.r_inv = powmod(self.r, p-2, p)
        p_inv = 1
        for i in range(int(log(r_bits)/log(2)+1)):
            p_inv = (p_inv * (2 - self.p * p_inv)) & self.r_mod_mask
        self.p_inv = p_inv
        
        self.ih = self.p_inv >> self.r_bits
        self.il = self.p_inv & ((1 << (self.r_bits//2)) - 1)
        self.ph = self.p >> self.r_bits
        self.pl = self.p & ((1 << (self.r_bits//2)) - 1)


    def repr(self, x):
        return x * self.r % self.p

    def reduce(self, x):
      
        q = ((x & self.r_mod_mask) * self.p_inv) & self.r_mod_mask
        a = (x  - q * self.p) >> self.r_bits
       
        if a < 0:
            a += self.p
        return a

    def mult(self, x, y):
        return self.reduce(x*y)
        

#p = mpz("0x FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF 93C467E3 7DB0C7A4 D1BE3F81 0152CB56 A1CECC3A F65CC019 0C03DF34 709AFFBD 8E4B59FA 03A9F0EE D0649CCB 621057D1 1056AE91 32135A08 E43B4673 D74BAFEA 58DEB878 CC86D733 DBE7BF38 154B36CF 8A96D156 7899AAAE 0C09D4C8 B6B7B86F D2A1EA1D E62FF864 3EC7C271 82797722 5E6AC2F0 BD61C746 961542A3 CE3BEA5D B54FE70E 63E6D09F 8FC28658 E80567A4 7CFDE60E E741E5D8 5A7BD469 31CED822 03655949 64B83989 6FCAABCC C9B31959 C083F22A D3EE591C 32FAB2C7 448F2A05 7DB2DB49 EE52E018 2741E538 65F004CC 8E704B7C 5C40BF30 4C4D8C4F 13EDF604 7C555302 D2238D8C E11DF242 4F1B66C2 C5D238D0 744DB679 AF289048 7031F9C0 AEA1C4BB 6FE9554E E528FDF1 B05E5B25 6223B2F0 9215F371 9F9C7CCC 69DDF172 D0D62342 17FCC003 7F18B93E F5389130 B7A661E5 C26E5421 4068BBCA FEA32A67 818BD307 5AD1F5C7 E9CC3D17 37FB2817 1BAF84DB B6612B78 81C1A48E 439CD03A 92BF5222 5A2B38E6 542E9F72 2BCE15A3 81B5753E A8427633 81CCAE83 512B3051 1B32E5E8 D8036214 9AD030AA BA5F3A57 98BB22AA 7EC1B6D0 F17903F4 E234EA60 34AA8597 3F79A93F FB82A75C 47C03D43 D2F9CA02 D03199BA CEDDD453 34DBC6B5 FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF")