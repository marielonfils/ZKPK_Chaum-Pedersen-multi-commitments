from hashlib import sha256
from random import SystemRandom
from constants import p,r
import numpy as np
from gmpy2 import mpz


def generate_generators(n):
    """p=23
    r=2
    q=11"""
    l=np.zeros(n,dtype=object)
    n_generated=0
    while n_generated != n:
        s=SystemRandom().randint(1,p)
        #m = sha256(s).hexdigest()
        t = pow(s,r,p)
        if t!=1:
            l[n_generated]=mpz(t)
            #print("mpz(",t,") ,")
            n_generated+=1
    """for i in range(q):
        print(pow(list[0],i,p))"""
    return l

#128
gs=generate_generators(128)
print("gs")
hs=generate_generators(128)
print("hs")
g_bolds=generate_generators(128*128)
print("g_bolds")
h_bolds=generate_generators(128*128)
print("h_bolds")
np.savez("constants2.npz", gs,hs,g_bolds,h_bolds)
