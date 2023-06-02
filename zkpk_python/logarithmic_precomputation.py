import numpy as np
import sys
from utils.group import random_exp, compute_size, PowRadix
from utils.hash2 import hash_elems
from bulletproof import bullet_proof, bullet_verification
from extended_schnorr import extended_schnorr_proof, extended_schnorr_verification
from utils.constants import generate_constants,h,u,p,q 
from gmpy2 import powmod as pow,mpz
import time
#import cProfile
#import pstats
start_time = time.time()


def delta(y,z):
    y_m=mpz(0)
    y_m_i=mpz(1)
    z_2=pow(z,2,q)
    z_m_i=z_2*z % q
    z_m=0
    for _ in range(m):
        y_m=(y_m +y_m_i)%q
        y_m_i= y_m_i* y % q
        z_m= (z_m+z_m_i)%q
        z_m_i= z_m_i * z %q
    return ((y_m*(z-z_2)%q )%q- z_m)  %q

#Prover's function
def generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas):
    #step 1#
    #1.1
    a_Ls=vs.copy()
    #verify if modulo 2 TODO
    #1.2
    a_Rs=(a_Ls-np.ones((l,m)))
    #1.3
    alpha=random_exp()
    rho=random_exp()
    s_Ls=np.zeros((l,m),dtype=object)
    s_Rs=np.zeros((l,m),dtype=object)
    for j in range(l):
        for i in range(m):
            s_Ls[j][i]=random_exp()
            s_Rs[j][i]=random_exp()
    #1.4 - 1.5
    A= hpow(alpha)
    S= hpow(rho) 
    for i in range(l):
        for j in range(m):
            A= A*gsboldspow[i][j](mpz(a_Ls[i][j]))%p 
            A= A*hsboldspow[i][j](mpz(a_Rs[i][j]))%p 
            S= S*gsboldspow[i][j](s_Ls[i][j])%p 
            S= S*hsboldspow[i][j](s_Rs[i][j])%p 
    #step 2 - 3 - 4#
    y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while y==0:
        y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1
    z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,q=q))
    rand =0
    while z==0:
        z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,rand,q=q))
    #step 4
    zs_1=np.ones(m,dtype=object)
    zs_m=np.zeros(m,dtype=object)
    zs_m_i= pow(z,2,q)
    
    ys_m=np.zeros(m,dtype=object)
    ys_m_i=mpz(1)    
    ys_m_sum=mpz(0)
    for i in range(m):
        zs_1[i]=z
        zs_m[i]=zs_m_i
        ys_m[i]=ys_m_i
        ys_m_sum = (ys_m_sum + ys_m_i) %q
        zs_m_i= zs_m_i *z % q
        ys_m_i= ys_m_i *y %q

    t1s=np.zeros((l),dtype=object)
    t2s=np.zeros((l),dtype=object)
    for j in range(l):
        sum_t1=mpz(0)
        sum_t2=mpz(0)
        for i in range(m):
            s_Rs_ys_m=ys_m[i]*s_Rs[j][i]%q
            sum_t1= (sum_t1 +((mpz(a_Ls[j][i])-z)%q)*s_Rs_ys_m % q)%q
            sum_t1 = (sum_t1+s_Ls[j][i]*zs_m[i]%q+(((mpz(a_Rs[j][i])+z)%q*ys_m[i])%q)*s_Ls[j][i]%q)%q
            sum_t2= (sum_t2+ s_Ls[j][i]*s_Rs_ys_m%q)%q
        t1s[j]=sum_t1
        t2s[j]=sum_t2

    #step 5#
    #5.1
    tau1 = random_exp()
    tau2 = random_exp()
    #5.2
    T1=hpow(tau1)
    T2=hpow(tau2)
    
    for i in range(l):
        T1 = T1*gspow[i](mpz(t1s[i]))%p 
        T2 = T2*gspow[i](mpz(t2s[i]))%p 

    #step 6 - 7 -8#
    x=mpz(hash_elems(T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while x==0:
        x=mpz(hash_elems(T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1
    #step 9#
    #9.1 - 9.2
    ls=np.zeros((l,m),dtype=object)
    rs=np.zeros((l,m),dtype=object)
    t_hats=np.zeros(l,dtype=object)
    x_2=pow(x,2,q)
    for j in range(l):        
        rs_i=np.zeros(m,dtype=object)
        sum=mpz(0)
        for i in range(m):
            ls_i_i=((mpz(a_Ls[j][i])-z)%q+s_Ls[j][i]*x %q)%q
            ls[j][i]=ls_i_i
            rhs_i=((mpz(a_Rs[j][i])+z)%q+s_Rs[j][i]*x%q)%q
            rs_i_i= ((rhs_i*ys_m[i])%q + zs_m[i]) %q        
            rs_i[i]=rs_i_i
            sum = (sum +ls_i_i*rs_i_i %q )%q
        rs[j,:]=rs_i        
        t_hats[j]=sum

    #9.4
    tau_x= ((tau2*x_2 %q + tau1*x %q)%q)%q
    for i in range(m):
        tau_x= (tau_x+ zs_m[i]*gammas[i]) % q
    #9.5
    mu=(alpha+rho*x % q)%q
    #step 10 - 11 - 12#
    phi=mpz(hash_elems(tau_x,mu,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while phi==0:
        phi=mpz(hash_elems(tau_x,mu,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1

    #step 15#
    #15.1 -15.4
    h_bolds_prime=np.zeros((l,m),dtype=object)
    g_bolds_prime=np.zeros((l,m),dtype=object)    
    phi_i=mpz(1)
    y_1=pow(y,q-2,q)
    phi_1=pow(phi,q-2,q)
    t_bar = mpz(0)
    phis_l=np.zeros(l,dtype=object)
    phi_j=mpz(1) 
    delta= ys_m_sum* (z-pow(z,2,q)) % q
    P_bar = hpow(q-tau_x) *pow(T1,x,p) % p * pow(T2,x_2,p) % p
    phis_l_1m=np.zeros((l,m),dtype=object)
    for i in range(m-1):
        delta = (delta - zs_m[i+1]) %q
        P_bar = P_bar * pow(Vs[i],zs_m[i],p) % p
    delta=(delta-zs_m[m-1]*z %q)%q
    P_bar = P_bar * pow(Vs[m-1],zs_m[m-1],p) % p
    P_ext=P_bar 
    P=A*pow(S,x,p) % p        
    for j in range(l):
        y_i=mpz(1)
        t_bar= (t_bar+t_hats[j]*phi_j%q) % q
        phis_l[j]=phi_j
        term= gspow[j](delta)  
        P_bar=  P_bar*term %p
        P_ext= P_ext*term%p
        P_bar= P_bar* hspow[j](phi_j)%p    
        for i in range(m):
            h_ij= hsboldspow[j][i](y_i) 
            g_ij= gsboldspow[j][i](phi_i) 
            h_bolds_prime[j][i]= h_ij
            g_bolds_prime[j][i]= g_ij
            P= P* pow(g_ij,mpz(q-z*phis_l[j] % q),p) % p
            P= P* pow(h_ij, mpz((z*ys_m[i]%q+zs_m[i]) %q),p) % p 
            phis_l_1m[j][i]=phi_j*ls[j][i]%q 
            y_i=y_i*y_1 % q
        phi_i=phi_i*phi_1 % q  
        phi_j= phi_j *phi % q        
            
    #step 16#
    a_1=phis_l_1m
    b_1=np.array(rs)
    g_bolds=np.array(g_bolds)
    h_bolds=np.array(h_bolds)
    t_hats=np.array(t_hats)
    bp1=bullet_proof(g_bolds_prime.flatten(),h_bolds_prime.flatten(), P*hpow(q-mu)%p,u,a_1.flatten(),b_1.flatten(),t_bar,p,q)
    bp2=bullet_proof(gs,hs,P_bar,u,t_hats.flatten(),phis_l.flatten(),t_bar,p,q)
    #step 17#
    e1=extended_schnorr_proof(gs,P_ext,t_hats,p,q)  
    return A,S,T1,T2,tau_x,mu,phi,t_bar,bp1,bp2,e1
    

#Verifier's function
def verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,T1,T2,tau_x,mu,phi,t_bar,bp_1,bp_2,e1):

    y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while y==0:
        y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1
    z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,q=q))
    rand =0
    while z==0:
        z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,rand,q=q))
    x=mpz(hash_elems(T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while x==0:
        x=mpz(hash_elems(T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1

    zs_m=np.zeros(m,dtype=object)
    zs_m_i= pow(z,2,q)
    ys_m=np.zeros(m,dtype=object)
    ys_m_i=mpz(1)
    x_2=pow(x,2,q)
    for i in range(m):
        zs_m[i]=zs_m_i
        ys_m[i]=ys_m_i
        zs_m_i=zs_m_i*z % q
        ys_m_i=ys_m_i*y %q


    #step 15#
    #15.1 -15.2
    h_bolds_prime=np.zeros((l,m),dtype=object)
    g_bolds_prime=np.zeros((l,m),dtype=object)
    y_1=pow(y,q-2,q)
    phi_i=mpz(1)
    phi_1=pow(phi,q-2,q)
    phis_l=np.zeros(l,dtype=object)
    phi_i_l=mpz(1)
    for j in range(l):
        phis_l[j]=phi_i_l
        phi_i_l= phi_i_l*phi %q
        y_i=mpz(1)
        for i in range(m):
            h_bolds_prime[j][i]= hsboldspow[j][i](y_i)
            g_bolds_prime[j][i]= gsboldspow[j][i](phi_i)
            y_i=y_i*y_1 % q
        phi_i=phi_i*phi_1 % q
    #15.3 -15.4
    delt= delta(y,z)
    P=A*pow(S,x,p) % p
    P_bar = hpow(q-tau_x)* pow(T1,x,p) % p * pow(T2,x_2,p) % p
    P_bar_e1=P_bar
    for i in range(m):
        P_bar = P_bar*pow(Vs[i],zs_m[i],p) % p
        P_bar_e1=P_bar
    for j in range(l):    
        term= gspow[j](delt)     
        P_bar= P_bar*term% p
        P_bar_e1=P_bar_e1*term %p
        P_bar= P_bar*pow(hs[j],phis_l[j],p) % p
        for i in range(m):
            P=  P*pow(g_bolds_prime[j][i], q-z*phis_l[j] % q, p) %p
            P=  P*pow(h_bolds_prime[j][i], (z*ys_m[i]%q+zs_m[i]) % q,p)%p    
    
    #step 16#
    a_bp1,b_bp1,Ls_bp1,Rs_bp1=bp_1
    if not bullet_verification(g_bolds_prime.flatten(),h_bolds_prime.flatten(), P*hpow(q-mu)%p,
                               u,a_bp1,b_bp1,t_bar,Ls_bp1,Rs_bp1,p,q):
        return 1
    
    a_bp2,b_bp2,Ls_bp2,Rs_bp2=bp_2
    if not  bullet_verification(gs,hs,P_bar,u,a_bp2,b_bp2,t_bar,Ls_bp2,Rs_bp2,p,q):
        return 2    

    #step 17#
    a_e1,Ls_e1,Rs_e1=e1
    if not extended_schnorr_verification(gs,P_bar_e1,a_e1,Ls_e1,Rs_e1,p,q):
        return 3
    return 0




n_false=0
lm=[1]#1,2,4,8,16,32,64,128]
ks=[1]#2,4,6,8,10,1]#[1,2,4,6,8,10]
k=sys.argv[1]
l=sys.argv[2]
lm.append(int(l))
ks.append(int(k))
f1="1_1_p.stats"
f2="1_1_v.stats"
generation_t=[]
prover_t=[]
verifier_t=[]
t=[]
size=[]

for k in ks:
    generation_t=[]
    prover_t=[]
    verifier_t=[]
    t=[]
    size=[]
    for l in lm:
    
    
        # groups size
        m=l#128

        # generators
        t0=time.time()
        gs,hs,g_bolds,h_bolds = generate_constants(l,m)
        #precomputation
        gs_radix = [PowRadix(gs[i], k) for i in range(l)]
        gspow = [gs_radix[i].pow for i in range(l)]
        hs_radix = [PowRadix(hs[i], k) for i in range(l)]
        hspow = [hs_radix[i].pow for i in range(l)]
        gsbolds_radix = [[PowRadix(g_bolds[i][j], k) for j in range(m)]for i in range(l)]
        gsboldspow = [[gsbolds_radix[i][j].pow for j in range(m)] for i in range(l)]
        hsbolds_radix = [[PowRadix(h_bolds[i][j], k) for j in range(m)]for i in range(l)]
        hsboldspow = [[hsbolds_radix[i][j].pow for j in range(m)] for i in range(l)]
        h_radix = PowRadix(h, k) 
        hpow = h_radix.pow
        t01=time.time()

        # votes and commitment on votes
        Vs=np.zeros(m,dtype=object)
        vs=np.random.randint(2,size=(l,m))
        t02=time.time()
        gammas=np.zeros(m,dtype=object)
        for i in range(m):
            gammas_i=random_exp()
            gammas[i]=gammas_i
            product=hpow(gammas_i)
            for j in range(l):
                product= product*gspow[j](mpz(vs[j][i]))%p
            Vs[i]=product
        t00 = time.time()
        n_false=0
        generation_t.append(t00-t02+t01-t0)
        p_t=[]
        v_t=[]
        for i in range(1):
            t1=time.time()
            #profiler = cProfile.Profile()
            #profiler.enable()
            A,S,T1,T2,tau_x,mu,phi,t_bar,bp1,bp2,e1=generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas)
            #profiler.disable()
            #profiler.dump_stats("p_1_1")
            t2=time.time()
            size.append(compute_size(A,S,T1,T2,tau_x,mu,phi,t_bar,bp1,bp2,e1,l,m))

            p_t.append( t2-t1)
            #profiler2 = cProfile.Profile()
            #profiler2.enable()
            t=verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,T1,T2,tau_x,mu,phi,t_bar,bp1,bp2,e1)
            #profiler2.disable()
            #profiler2.dump_stats("v_1_1")
            t3=time.time()
            v_t.append( t3-t2)
            if t!=0:
                n_false+=1
                print(t)
        print(k,l,m,n_false,t00-t02+t01-t0,p_t,v_t)
        prover_t.append(np.mean(p_t))
        verifier_t.append(np.mean(v_t))
    print("k {} l {}".format(k,l))
    print("generation_t=",generation_t)
    print("prover_t=",prover_t)
    print("verifier_t=",verifier_t)
    print("size=",size)
    print("--- %s seconds ---" % (time.time() - start_time))
#stats=pstats.Stats("p_1_1")
#stats.print_callers("powmod")
#stats=pstats.Stats("v_1_1")
#stats.print_callers("powmod")