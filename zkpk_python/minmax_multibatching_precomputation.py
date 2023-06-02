import numpy as np
from utils.group import random_exp,compute_size2, PowRadix
from utils.hash2 import hash_elems
import json
from bulletproof import bullet_proof, bullet_verification
from extended_schnorr import extended_schnorr_proof, extended_schnorr_verification
from utils.constants import generate_constants,g,h,u,p,q 
from gmpy2 import powmod as pow,mpz
import time
#import cProfile
#import pstats
start_time = time.time()


#Prover's function
#@profile
def generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas):
    #step 1#
    #1.1 - 1.2
    a_Ls=np.zeros((n,m),dtype=object)
    a_Rs=np.zeros((n,m),dtype=object)    
    for i in range(m):
        v_j =0 
        for k in range(l):
            v_j = (v_j+vs[k][i])%q

        shifted = v_j
        for j in range(n):
            a_i=shifted&1
            a_Ls[j][i]=a_i
            a_Rs[j][i]=a_i-1
            shifted = shifted>>1    
    #1.3
    alpha=random_exp()
    rho=random_exp()
    s_Ls=np.zeros((n,m),dtype=object)
    s_Rs=np.zeros((n,m),dtype=object)
    for j in range(n):
        for i in range(m):
            s_Ls[j][i]=random_exp()
            s_Rs[j][i]=random_exp()
    #1.4 - 1.5
    A=hpow(alpha)
    S=hpow(rho)
    twos_n=np.zeros(n,dtype=object)
    twos_i=mpz(1)
    twos_sum=mpz(0)           
    for i in range(n):
        twos_n[i]=twos_i
        twos_sum = (twos_sum+twos_i)%q
        twos_i = twos_i*mpz(2)%q
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
    
    #step 5#
    #5.1
    tau0 = random_exp()
    tau1 = random_exp()
    tau2 = random_exp()

    #5.2
    zs_m=np.zeros(m,dtype=object)
    zs_m_sum=mpz(0)
    zs_m_i= pow(z,2,q)
    vs_bar=np.zeros(l,dtype=object)
    c_bar = mpz(0)      
    for i in range(m):
        zs_m[i]=zs_m_i
        zs_m_sum = (zs_m_sum+zs_m_i)%q       
        zs_m_i= zs_m_i *z % q
    for k in range(l):
        c_sum=mpz(0)
        for i in range(m):
            c_sum= (c_sum+vs[k][i]*zs_m[i]%q)%q
            c_bar=(c_bar+vs[k][i]*zs_m[i]%q)%q
        vs_bar[k]=c_sum
    #5.3 - 5.4
    t1=mpz(0)
    t2=mpz(0)

    h_bolds_prime=np.zeros((n,m),dtype=object)
    y_1=pow(y,q-2,q)
    y_i=mpz(1)  
    P=A   
    exp_z = q-z 
    t0=mpz(0)
    ys_m=np.zeros(n*m,dtype=object)
    ys_m_i=mpz(1)    
    ys_m_sum=mpz(0) 
    for j in range(n):        
        for i in range(m):
            ind=j*m+i
            ys_m[j*m+i]=ys_m_i
            ys_m_sum = (ys_m_sum+ys_m_i)%q
            s_Rs_ys_m=ys_m_i*s_Rs[j][i]%q
            t1= (t1 +((mpz(a_Ls[j][i])-z)%q)*s_Rs_ys_m % q)%q
            t1 = (t1+s_Ls[j][i]*zs_m[i]%q*twos_n[j]%q+(((mpz(a_Rs[j][i])+z)%q*ys_m_i)%q)*s_Ls[j][i]%q)%q
            t2= (t2+ s_Ls[j][i]*s_Rs_ys_m%q)%q  
            

            ind=j*m+i
            h_ij= hsboldspow[j][i](y_i) 
            h_bolds_prime[j][i]= h_ij
            P= P*gsboldspow[j][i](exp_z)%p 
            P= P* pow(h_bolds_prime[j][i],mpz((z*ys_m_i%q+zs_m[i]*twos_n[j]%q) %q),p) % p 
            y_i=y_i*y_1 % q  
            l1=(a_Ls[j][i]-z)%q
            r1=(a_Rs[j][i]+z)%q*ys_m_i%q+zs_m[i]*twos_n[j]
            t0 = (t0+l1*r1%q)%q 
            ys_m_i= ys_m_i *y %q
    #5.5
    T0=hpow(tau0)*gpow(c_bar)%p
    T1=hpow(tau1)*gpow(t1)%p
    T2=hpow(tau2)*gpow(t2)%p

    #step 6 - 7 -8#
    x=mpz(hash_elems(T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while x==0:
        x=mpz(hash_elems(T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1    

    #step 9#
    #9.1 - 9.3
    ls=np.zeros((n,m),dtype=object)
    rs=np.zeros((n,m),dtype=object)
    t_hat=mpz(0)
    x_2=pow(x,2,q)
    for j in range(n):        
        for i in range(m):
            ind = j*m+i
            ls_i_i=((mpz(a_Ls[j][i])-z)%q+s_Ls[j][i]*x %q)%q
            ls[j][i]=ls_i_i
            rhs_i=((mpz(a_Rs[j][i])+z)%q+s_Rs[j][i]*x%q)%q
            rs_i_i= ((rhs_i*ys_m[ind])%q + zs_m[i]*twos_n[j]) %q        
            rs[j][i]=rs_i_i
            t_hat = (t_hat +ls_i_i*rs_i_i %q )%q
    #9.4
    tau_x= ((tau2*x_2 %q + tau1*x %q)%q + tau0)%q
    #9.5
    mu=(alpha+rho*x % q)%q
    P= P*hpow(q-mu)%p*pow(S,x,p)%p

   

    #step 13#
    bp1=bullet_proof(g_bolds.flatten(),h_bolds_prime.flatten(),P,u,ls.flatten(),rs.flatten(),t_hat,p,q)

    #step 14 - step 15#
    #14.1 -14.3
    rho_prime = random_exp()
    ss=np.zeros(l,dtype=object)
    S_prime=hpow(rho_prime)#pow(h,rho_prime,p)
    
    #step 16#
    #16.1 
    tau1_prime = random_exp()
    #16.2
    t1_prime=mpz(0)
    for k in range(l):
        s_i = random_exp()
        ss[k]= s_i
        S_prime = S_prime * gspow[k](s_i)%p
        t1_prime = (t1_prime+ss[k])%q    
    #16.3
    T1_prime = hpow(tau1_prime)*gpow(t1_prime)%p

    #step 17 - 18 - 19#
    x_prime=mpz(hash_elems(T1_prime,S_prime,tau_x,mu,t_hat,T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while x_prime==0:
        x_prime=mpz(hash_elems(T1_prime,S_prime,tau_x,mu,t_hat,T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1
    
    #step 20#
    #20.1 - 20.2
    ls_prime = np.zeros(l,dtype=object)
    rs_prime = np.zeros(l,dtype=object)
    for k in range(l):
        ls_prime[k]=(vs_bar[k]+ss[k]*x_prime%q)%q
        rs_prime[k]=mpz(1)
    #20.3
    t_prime = (c_bar + t1_prime*x_prime%q)%q 
    #20.4
    tau_x_prime = (tau0 + tau1_prime*x_prime%q)%q
    #20.5
    mu_prime = rho_prime*x_prime%q
    for i in range(m):
        mu_prime = (mu_prime+gammas[i]*zs_m[i]%q)%q
    
    #step 23#
    P2=hpow(q-mu_prime)*pow(S_prime,x_prime,p)%p
    for i in range(m):
        P2 = P2*pow(Vs[i],zs_m[i],p)%p
    
    P3=P2
    for k in range(l):
        P2 = P2 * hs[k]%p  
   
    bp2=bullet_proof(gs,hs,P2,u,ls_prime,rs_prime,t_prime,p,q)

    #step 24#
    e1=extended_schnorr_proof(gs,P3,ls_prime,p,q)  
    return A,S,T0,T1,T2,tau_x,mu,t_hat,bp1,S_prime,T1_prime,mu_prime,tau_x_prime,t_prime,bp2,e1
    

#Verifier's function
#@profile
def verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,T0,T1,T2,tau_x,mu,t_hat,bp_1,S_prime,T1_prime,mu_prime,tau_x_prime,t_prime,bp_2,e1):

    y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while y==0:
        y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1
    z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,q=q))
    rand =0
    while z==0:
        z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,rand,q=q))
        x=mpz(hash_elems(T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    x=mpz(hash_elems(T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    while x==0:
        x=mpz(hash_elems(T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1 
    x_prime=mpz(hash_elems(T1_prime,S_prime,tau_x,mu,t_hat,T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=q))
    rand=0
    while x_prime==0:
        x_prime=mpz(hash_elems(T1_prime,S_prime,tau_x,mu,t_hat,T0,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=q))
        rand+=1

    zs_m=np.zeros(m,dtype=object)
    zs_m_sum=mpz(0)
    zs_m_i= pow(z,2,q)
    twos_n=np.zeros(n,dtype=object)
    twos_i=mpz(1)
    twos_sum=mpz(0)
    ys_m=np.zeros(n*m,dtype=object)
    ys_m_i=mpz(1)    
    ys_m_sum=mpz(0)

    #5.2
    for i in range(n*m):
        if i < n:
            twos_n[i]=twos_i
            twos_sum = (twos_sum+twos_i)%q
            twos_i = twos_i*mpz(2)%q
        ys_m[i]=ys_m_i
        ys_m_sum = (ys_m_sum+ys_m_i)%q
        ys_m_i= ys_m_i *y %q
    for i in range(m):
        zs_m[i]=zs_m_i
        zs_m_sum = (zs_m_sum+zs_m_i)%q       
        zs_m_i= zs_m_i *z % q       
    
    #step 12
    delta= ys_m_sum* (z-pow(z,2,q)) % q
    delta = (delta - twos_sum*zs_m_sum%q*z%q)%q
    lhs = gpow(t_hat)*hpow(tau_x)%p 
    rhs= T0*pow(T1,x,p)%p*pow(T2,pow(x,2,q),p)%p*gpow(delta)%p 
    if lhs != rhs:
        return 1
    
    #step 11#
    #11.1 
    h_bolds_prime=np.zeros((n,m),dtype=object)
    y_1=pow(y,q-2,q) 
    #11.1 -11.2
    P=A*pow(S,x,p) % p   
    y_i=mpz(1)  
    exp_z = mpz(q-z)   
    for j in range(n):                
        for i in range(m):
            ind=j*m+i
            h_ij=hsboldspow[j][i](y_i) 
            h_bolds_prime[j][i]= h_ij
            P= P*gsboldspow[j][i](exp_z)%p 
            P= P* pow(h_bolds_prime[j][i],mpz((z*ys_m[ind]%q+zs_m[i]*twos_n[j]%q) %q),p) % p             
            y_i=y_i*y_1 % q
    P= P*hpow(q-mu)%p 

    #step 13#
    a_bp1,b_bp1,Ls_bp1,Rs_bp1=bp_1
    if not bullet_verification(g_bolds.flatten(),h_bolds_prime.flatten(),
                               P,u,a_bp1,b_bp1,t_hat,Ls_bp1,Rs_bp1,p,q):
        return 2
    
    #step 22#
    lhs2 = gpow(t_prime)*hpow(tau_x_prime)%p 
    rhs2 = T0 * pow(T1_prime,x_prime,p)%p
    if rhs2 != lhs2:
        return 3
    
    #step 23#
    P2= hpow(q-mu_prime)*pow(S_prime,x_prime,p)%p 
    for i in range(m):
        P2 = P2*pow(Vs[i],zs_m[i],p)%p
    P3=P2
    for k in range(l):
        P2 = P2 * hs[k]%p
    a_bp2,b_bp2,Ls_bp2,Rs_bp2=bp_2
    if not  bullet_verification(gs,hs,P2,u,a_bp2,b_bp2,t_prime,Ls_bp2,Rs_bp2,p,q):
        return 4   

    #step 24#
    a_e1,Ls_e1,Rs_e1=e1
    if not extended_schnorr_verification(gs,P3,a_e1,Ls_e1,Rs_e1,p,q):
        return 5
    return 0


n_false=0
lm=[1]#,2,4,8,16,32,64,128]
ks=[10]#,2,4,6,8,10]
f1="1_1_p.stats"
f2="1_1_v.stats"
generation_t=[]
prover_t=[]
verifier_t=[]
t=[]
size=[]
for l in lm:
    for k in ks:
    
        # groups size
        #l=128
        n=2
        m=l    
        # generators
        time0=time.time()
        gs,hs,g_bolds,h_bolds = generate_constants(l,m,n)
        #precomputation
        gs_radix = [PowRadix(gs[i], k) for i in range(l)]
        gspow = [gs_radix[i].pow for i in range(l)]
        gsbolds_radix = [[PowRadix(g_bolds[i][j], k) for j in range(m)]for i in range(n)]
        gsboldspow = [[gsbolds_radix[i][j].pow for j in range(m)] for i in range(n)]
        hsbolds_radix = [[PowRadix(h_bolds[i][j], k) for j in range(m)]for i in range(n)]
        hsboldspow = [[hsbolds_radix[i][j].pow for j in range(m)] for i in range(n)]
        g_radix = PowRadix(h, k) 
        gpow = g_radix.pow
        h_radix = PowRadix(g, k) 
        hpow = h_radix.pow
        time01=time.time()
        # votes and commitment on votes
        Vs=np.zeros(m,dtype=object)
        vs=np.random.randint(2,size=(l,m))
        limit=np.zeros(m)
        for i in range(m):
            for j in range(l):
                if limit[i]>=2**(n-1):
                    vs[j][i]=0
                else:
                    if vs[j][i]==1:
                        limit[i]+=1
        time02=time.time()
        gammas=np.zeros(m,dtype=object)
        for i in range(m):
            gammas_i=random_exp()
            gammas[i]=gammas_i
            product=hpow(gammas_i)
            for j in range(l):
                product= product*gspow[j](mpz(vs[j][i]))%p 
            Vs[i]=product
        t00 = time.time()
        generation_t.append(t00-time02+time01-time0)
        v_t=[]
        p_t=[]
        print("#generated")
        for i in range(3):
            t1=time.time()
            #profiler = cProfile.Profile()
            #profiler.enable()
            A,S,T0,T1,T2,tau_x,mu,t_hat,bp1,S_prime,T1_prime,mu_prime,tau_x_prime,t_prime,bp2,e1=generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas)
            #profiler.disable()
            #profiler.dump_stats(f1)
            t2=time.time()
            #print("prover time:", t2-t1)
            p_t.append(t2-t1)
            size.append(compute_size2(A,S,T0,T1,T2,tau_x,mu,t_hat,bp1,S_prime,T1_prime,mu_prime,tau_x_prime,t_prime,bp2,e1,l,m,n))
            t3=time.time()
            #profiler2 = cProfile.Profile()
            #profiler2.enable()
            t=verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,T0,T1,T2,tau_x,mu,t_hat,bp1,S_prime,T1_prime,mu_prime,tau_x_prime,t_prime,bp2,e1)
            #profiler2.disable()
            #profiler2.dump_stats(f2)
            t4=time.time()
            #print("verifier time:", t4-t3)
            v_t.append(t4-t3)
            if t!=0:
                n_false+=1
                print("error",t)
        print(k,l,m,n,n_false,p_t,v_t)
        prover_t.append(np.mean(p_t))
        verifier_t.append(np.mean(v_t))
print("generation_t=",generation_t)
print("prover_t=",prover_t)
print("verifier_t=",verifier_t)
print("size=",size)
print("--- %s seconds ---" % (time.time() - start_time))
#stats=pstats.Stats(f1)
#stats.print_callers("powmod")
#stats=pstats.Stats(f2)
#stats.print_callers("powmod")