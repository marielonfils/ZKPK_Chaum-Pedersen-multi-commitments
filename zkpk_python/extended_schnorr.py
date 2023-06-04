import numpy as np
from utils.hash2 import hash_elems
from gmpy2 import powmod as pow,mpz

### Implementation of the Extended Schnorr argument###

#Prover's function
def extended_schnorr_proof(g,P,a,p,q):
    """
    Provides a zero knowledge proof of the relation : g[0]^a[0]...g[n-1]^a[n-1]=P
    inputs
        - g : array of mpz of size n containing elements of a group modulo p
        - a : array of mpz of size n, element of a group modulo q
        - P : mpz, element of a group modulo p
        - p : mpz, prime value, all elements are modulo p
        - q : mpz, exponents are modulo q
        Note that n should be a power of 2.
    outputs:
        - a_s : array of mpz of size log_2(n) containing elements of a group modulo q
        - bs  : array of mpz of size log_2(n) containing elements of a group modulo q
        - Ls  : array of mpz of size log_2(n) containing elements of a group modulo p
        - Rs  : array of mpz of size log_2(n) containing elements of a group modulo p
    """
    Ls=[]
    Rs=[]
    def proof(g,P,a,n):

        if n==1:
            return (a,Ls,Rs)   

        #step 1#
        n_prime=int(n/2)
        a_l=a[:n_prime]
        a_r=a[n_prime:]
        g_l=g[:n_prime]
        g_r=g[n_prime:]

        L=1
        R=1
        for i in range(n_prime):
            L=L*pow(g_r[i],a_l[i],p)%p
            R=R*pow(g_l[i],a_r[i],p)%p
        Ls.append(L)
        Rs.append(R)

        #step 2-3-4 #
        x=mpz(hash_elems(L,R,g,P, q=q))
        rand=0
        while x==0:
            x=mpz(hash_elems(L,R,g,P, q=q))
            rand+=1
        x_1= pow(x, q-2,q)

        #step 5#
        g_prime=[]
        a_prime=[]
        x2=pow(x,2,q)
        x_2=pow(x,q-3,q)
        for i in range(n_prime):
            g_prime.append((pow(g_r[i],x,p)*pow(g_l[i],x_1,p))% p)
            a_prime.append(((a_r[i]*x_1)%q+ (a_l[i]*x)%q)%q)   
        
        P_prime=(pow(L,x2,p)*P)%p * pow(R,x_2,p) %p
        
        return proof(g_prime,P_prime,a_prime,n_prime)
    n=len(g)
    if not((n != 0) and (n & (n-1) == 0)):
        raise ValueError("basis are not a power of 2")
    if len(a)!=n :
        raise ValueError("vectors have not the same length")
    return proof(g,P,a,n)

#Verifier's function
def extended_schnorr_verification(g,P,a,Ls,Rs,p,q):
    """
    Verifies if g[0]^c[0]...g[n-1]^c[n-1]=P
    inputs
        - g : array of mpz of size n containing elements of a group modulo p
        - P : mpz element of a group modulo p
        - a : mpz element of a group modulo q
        - Ls: array of mpz of size log_2(n) containing elements of a group modulo p
        - Rs: array of mpz of size log_2(n) containing elements of a group modulo p
        - p : mpz, prime value, all elements are modulo p
        - q : mpz, exponents are modulo q
        Note that n should be a power of 2.
    outputs
        - True : success, proof verified
        - False: the proof is erroneous
        - ValueError : if len(Rs) != len(Ls) or if n is not a power of 2
    """
    n_prime=len(g)
    n_prime2=len(Ls)
    if not((n_prime != 0) and (n_prime & (n_prime-1) == 0)):
        raise ValueError("basis are not a power of 2")
    if 2**n_prime2!=n_prime:
        raise ValueError("basis are not a power of 2")
    if len(Rs)!=n_prime2:
        raise ValueError("vectors have not the same length")
    
    a_prime=a
    g_prime=g
    P_prime=P
    i=0
    while(n_prime >1):
        #step 1#
        n_prime=int(n_prime/2)

        #step 5#
        g_l=g_prime[:n_prime]
        g_r=g_prime[n_prime:]

        g_prime=[]
        rand=0
        x=mpz(hash_elems(Ls[i],Rs[i],g,P, q=q))
        while x==0:
            x=mpz(hash_elems(Ls[i],Rs[i],g,P, q=q))
            rand+=1
        x_1=pow(x,q-2,q)
        for j in range(n_prime):            
            g_prime.append((pow(int(g_r[j]),x,p)*pow(int(g_l[j]),x_1,p))% p)
        
        x2=pow(x,2,q)
        x_2=pow(x,q-3,q)
        P_prime=(pow(Ls[i],x2,p)*P_prime)%p * pow(Rs[i],x_2,p) %p
        i+=1
    if P_prime==pow(g_prime[0],a_prime[0],p):
        return True
    return False

def test():
    g=[mpz(3),mpz(4)]
    P=1 #3^2%23=9, 4^3%23=18, 9*18%23=1
    a=[mpz(2),mpz(3)]
    p=mpz(23)
    q=mpz(11)
    a,Ls,Rs=extended_schnorr_proof(g,P,np.array(a),p,q)
    assert extended_schnorr_verification(g,P,a,Ls,Rs,p,q)
#test()