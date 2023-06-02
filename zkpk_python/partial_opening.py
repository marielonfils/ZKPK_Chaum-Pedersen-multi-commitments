import numpy as np
from utils.group import random_exp, compute_size3
from utils.hash2 import hash_elems
from extended_schnorr import extended_schnorr_proof, extended_schnorr_verification
from utils.constants import generate_constants,h,u,p,q 
from gmpy2 import powmod as pow,mpz
import time
import cProfile
import pstats
start_time = time.time()

# groups size
l=2+1#128

# generators
t0=time.time()
gs = generate_constants(l)


# votes and commitment on votes
index_j=np.random.randint(l)
vs=np.random.randint(2,size=(l))
bit=vs[index_j]
V=mpz(1)
gamma=random_exp()
V = V*pow(h,gamma,p)%p
for j in range(l):    
    V = V *pow(gs[j],mpz(vs[j]),p)%p

t00 = time.time()

#Prover's function
def generate_proof(gs,h,u,V,index_j,bit,vs,gamma):

    #step 0#
    g_bolds=np.delete(gs,index_j)

    #step 1#
    #1.1
    vs_tilde=np.delete(vs,index_j)
    #1.2
    alpha=random_exp()
    ss=np.zeros(l-1,dtype=object)
    for i in range(l-1):
        ss[i]=random_exp()
    #1.3
    S=pow(h,alpha,p)
    for i in range(l-1):
        S = S * pow(g_bolds[i],ss[i],p)%p
    
    #step 2 - 3 - 4#
    x=mpz(hash_elems(S,gs,h,u,V,index_j,bit,q=q))
    rand=0
    while x==0:
        x=mpz(hash_elems(S,gs,h,u,V,index_j,bit,rand,q=q))
        rand+=1
    
    #step 5
    #5.1
    mu = (gamma+alpha*x%q)%q
    #5.2
    bs=np.zeros(l-1,dtype=object)
    for i in range(l-1):
        bs[i] = (vs_tilde[i]+ss[i]*x%q)%q
   
    #step 7#
    P = pow(h,q-mu,p)*pow(S,x,p)%p*V%p*pow(gs[index_j],q-bit,p)%p
     
    #step 8#
    e=extended_schnorr_proof(g_bolds,P,bs,p,q)  
    return S,mu,e
    

#Verifier's function
def verify_proof(gs,h,u,V,index_j,bit,S,mu,e):

    x=mpz(hash_elems(S,gs,h,u,V,index_j,bit,q=q))
    rand=0
    while x==0:
        x=mpz(hash_elems(S,gs,h,u,V,index_j,bit,rand,q=q))
        rand+=1

    #step 0#
    g_bolds = np.delete(gs,index_j)

    #step 7#
    P=pow(h,q-mu,p)*pow(S,x,p)%p*V%p*pow(gs[index_j],q-bit,p)%p  

    #step 8#
    a_e,Ls_e,Rs_e=e
    if not extended_schnorr_verification(g_bolds,P,a_e,Ls_e,Rs_e,p,q):
        return 1
    return 0


n_false=0
print("generation time:", t00-t0)
for i in range(1):
    t1=time.time()
    #profiler = cProfile.Profile()
    #profiler.enable()
    S,mu,e=generate_proof(gs,h,u,V,index_j,bit,vs,gamma)
    #profiler.disable()
    #profiler.dump_stats("p_1_1")
    t2=time.time()
    print("prover time:", t2-t1)
    print(compute_size3(S,mu,e,l-1,1))
    t3=time.time()
    #profiler2 = cProfile.Profile()
    #profiler2.enable()
    t=verify_proof(gs,h,u,V,index_j,bit,S,mu,e)
    #profiler2.disable()
    #profiler2.dump_stats("v_1_1")
    t4=time.time()
    print("verifier time:", t4-t3)
    if t!=0:
        n_false+=1
        print(t)
print(n_false)
print("--- %s seconds ---" % (time.time() - start_time))
#stats=pstats.Stats("p_1_1")
#stats.print_callers("powmod")
#stats=pstats.Stats("v_1_1")
#stats.print_callers("powmod")