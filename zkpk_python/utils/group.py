from utils.constants import p,q
from utils.montgomery import Montgomery
from gmpy2 import mpz_random,random_state,mpz, powmod
from sys import getsizeof
from time import time

def compute_size_array(a):
    size=0
    for i in range(len(a)):          
        size+= getsizeof(a[i])
    return size  

def compute_size(A,S,T1,T2,tau_x,mu,phi,t_bar,bp1,bp2,e1,l,m):
    a_bp1,b_bp1,Ls_bp1,Rs_bp1=bp1
    a_bp2,b_bp2,Ls_bp2,Rs_bp2=bp2
    a_e1,Ls_e1,Rs_e1=e1

    return getsizeof(A)+getsizeof(S)+getsizeof(T1)+getsizeof(T2)\
        +getsizeof(tau_x)+getsizeof(mu)+getsizeof(phi)+getsizeof(t_bar)\
        +getsizeof(a_bp1)+getsizeof(b_bp1)+compute_size_array(Ls_bp1)+compute_size_array(Rs_bp1)\
        +getsizeof(a_bp2)+getsizeof(b_bp2)+compute_size_array(Ls_bp2)+compute_size_array(Rs_bp2)\
        +getsizeof(a_e1)+compute_size_array(Ls_e1)+compute_size_array(Rs_e1)

def compute_size2(A,S,T0,T1,T2,tau_x,mu,t_hat,bp1,S_prime,T1_prime,mu_prime,tau_x_prime,t_prime,bp2,e1,l,m,n):
    a_bp1,b_bp1,Ls_bp1,Rs_bp1=bp1
    a_bp2,b_bp2,Ls_bp2,Rs_bp2=bp2
    a_e1,Ls_e1,Rs_e1=e1
    return getsizeof(A)+getsizeof(S)+getsizeof(T0)+getsizeof(T1)+getsizeof(T2)\
        +getsizeof(tau_x)+getsizeof(mu)+getsizeof(t_hat)\
        +getsizeof(S_prime)+getsizeof(T1_prime)\
        +getsizeof(mu_prime)+getsizeof(tau_x_prime)+getsizeof(t_prime)\
        +getsizeof(a_bp1)+getsizeof(b_bp1)+compute_size_array(Ls_bp1)+compute_size_array(Rs_bp1)\
        +getsizeof(a_bp2)+getsizeof(b_bp2)+compute_size_array(Ls_bp2)+compute_size_array(Rs_bp2)\
        +getsizeof(a_e1)+compute_size_array(Ls_e1)+compute_size_array(Rs_e1)
          

def compute_size3(S,mu,e,l,m):
    a_e1,Ls_e1,Rs_e1=e
    if l*m != 1 and l!=1:
        return getsizeof(S)+getsizeof(mu)\
          +getsizeof(a_e1)+getsizeof(Ls_e1[0])*len(Ls_e1)+getsizeof(Rs_e1[0])*len(Rs_e1)
    else:
         return getsizeof(S)+getsizeof(mu)\
          +getsizeof(a_e1)+getsizeof(Ls_e1)+getsizeof(Rs_e1)
        
def random_exp():
        """Take uniformly at random an integer in {0,...,q-1} and returns it.

        Therefore, pow(self.g, self.random_exp(), self.p) is a random element
        of the group.
        """
        rs = random_state(int(time() * 100000))
        return mpz_random(rs,q)

#from https://github.com/uclcrypto/many01proofs/blob/main/group.py

# Class supporting exponentiation with precomputation (radix method)
class PowRadix:
    def __init__(self, base, k=1, montgomery = False):

        self.montgomery = montgomery
        self.base=base

        if self.montgomery:
            self.M = Montgomery(p, 4096)
            base = self.M.repr(base)
            self.r_one = self.M.repr(mpz(1))

        # Precomputation in base
        self.table_length = -(-q.bit_length() // k)  # Double negative to take the ceiling
        self.k = k
        table = []
        row_base = base
        running_base = row_base
        for _ in range(self.table_length):
            if self.montgomery:
                row = [self.r_one]
            else:
                row = [1]
            for j in range(1, 2**k):
                row.append(running_base)
                if self.montgomery:
                    running_base = self.M.mult(running_base, row_base)
                else:
                    running_base = running_base * row_base % p
            table.append(row)
            row_base = running_base
        self.table = table
        self.max_exponent = 2**256

    def pow(self, e):
        # Computing pow(base, e, p)
        if not 0 <= e < self.max_exponent:
            return powmod(self.base,e,p)
        #assert 0 <= e < self.max_exponent, 'exponent out of range'
        y = mpz(1)
        if self.montgomery:
            y = self.r_one

        for i in range(self.table_length):
            e_slice = e[i * self.k: (i+1) * self.k]
            if self.montgomery:
                y = self.M.mult(y, self.table[i][e_slice]) 
            else:
                y = y * self.table[i][e_slice] % p
        if self.montgomery:
            y = self.M.reduce(y)
        assert y==powmod(self.base,e,p), "exponentiation problem"
        return y
