import numpy as np
import utils.elgamal as eg
from utils.hash import hashg 
from utils.hash2 import hash_elems
import json
from bulletproof import bullet_proof, bullet_verification
from extended_schnorr import extended_schnorr_proof, extended_schnorr_verification
from utils.constants import generate_constants,h,u,p,q 
from gmpy2 import powmod as pow, invert,mpz
import time
from sys import getsizeof
start_time = time.time()
#change QNbits in elgamal.py


# groups size
l=128
m=128

# generators
t0=time.time()
G=eg.ElgamalGroup(p,q)
"""h=eg.random_generator(G.p,G.q)
u=eg.random_generator(G.p,G.q)
gs=np.zeros(l)
hs=np.zeros(l)
g_bolds=np.zeros((l,m))
h_bolds=np.zeros((l,m))
for j in range(l):
    gs[j]=eg.random_generator(G.p, G.q)
    hs[j]=eg.random_generator(G.p, G.q)
    g_bolds_i=np.zeros(m)
    h_bolds_i=np.zeros(m)
    for i in range(m):
        g_bolds_i[i]=eg.random_generator(G.p, G.q)
        h_bolds_i[i]=eg.random_generator(G.p, G.q)
    g_bolds[j,:]=g_bolds_i
    h_bolds[j,:]=h_bolds_i"""
gs,hs,g_bolds,h_bolds = generate_constants(l,m)
"""gs=[mpz(2),mpz(2)]
hs=[mpz(2),mpz(2)]
g_bolds=[[mpz(3),mpz(3)],[mpz(3),mpz(3)]]
h_bolds=[[mpz(3),mpz(3)],[mpz(3),mpz(3)]]
h=mpz(3)
g=mpz(3)
u=mpz(2)"""

"""gs_list=gs.astype(int).tolist()
hs_list=hs.astype(int).tolist()
g_bolds_list=g_bolds.astype(int).tolist()
h_bolds_list=h_bolds.astype(int).tolist()"""

# votes and commitment on votes
Vs=np.zeros(m,dtype=object)
vs=np.random.randint(2,size=(l,m))
#TODO!!!
gammas=np.zeros(m,dtype=object)
for i in range(m):
    gammas_i=G.random_exp()
    gammas[i]=gammas_i
    product=pow(h,gammas_i,G.p)
    for j in range(l):
        product= product *pow(gs[j],mpz(vs[j][i]),G.p)%G.p
    Vs[i]=product
#Vs_list=Vs.astype(int).tolist()
t00 = time.time()

def delta(y,z):
    # delta(y,z) := (z-z^2) <1^m,y^m> - \sum_{k=1}^m z^{k+2}
    #TODO :verify if y^m correct
    y_m=mpz(0)
    y_m_i=mpz(1)
    z_2=pow(z,2,G.q)
    z_m_i=z_2*z % G.q
    z_m=0
    for _ in range(m):
        y_m=(y_m +y_m_i)%G.q
        y_m_i= y_m_i* y % G.q
        z_m= (z_m+z_m_i)%G.q
        z_m_i= z_m_i * z %G.q
    return ((y_m*(z-z_2)%G.q )%G.q- z_m)  %G.q

#Prover's function
def generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas):
    #step 1#
    #1.1
    a_Ls=vs.copy()
    #verify if modulo 2 TODO
    #1.2
    a_Rs=(a_Ls-np.ones((l,m)))
    #1.3
    alpha=G.random_exp()
    rho=G.random_exp()
    s_Ls=np.zeros((l,m),dtype=object)
    s_Rs=np.zeros((l,m),dtype=object)
    for j in range(l):
        for i in range(m):
            s_Ls[j][i]=G.random_exp()
            s_Rs[j][i]=G.random_exp()
    #1.4 - 1.5
    A=pow(h,alpha,G.p)
    S=pow(h,rho,G.p)
    for i in range(l):
        for j in range(m):
            A=A*pow(g_bolds[i][j],mpz(a_Ls[i][j]),G.p) % G.p
            A=A*pow(h_bolds[i][j],mpz(a_Rs[i][j]),G.p) % G.p
            S=S*pow(g_bolds[i][j],s_Ls[i][j],G.p) % G.p
            S=S*pow(h_bolds[i][j],s_Rs[i][j],G.p) % G.p
    #step 2 - 3 - 4#
    y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=G.q))
    #y = hashg(json.dumps({"A": A,"S":S, "gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list}),G.q)
    rand=0
    while y==0:
        y=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=G.q))
        #y = hashg(json.dumps({"A": A,"S":S, "gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list, "rand":rand}),G.q)
        rand+=1
    z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,q=G.q))
    #z = hashg(json.dumps({"A": A,"S": S,"gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list, "r":42}),G.q)
    rand =0
    while z==0:
        z=mpz(hash_elems(A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,42,rand,q=G.q))
        #z = hashg(json.dumps({"A": A,"S": S,"gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list, "r":42,"rand":rand}),G.q)
    #step 4
    zs_1=np.ones(m,dtype=object)
    zs_m=np.zeros(m,dtype=object)
    zs_m_i= pow(z,2,G.q)
    
    ys_m=np.zeros(m,dtype=object)
    ys_m_i=mpz(1)    
    ys_m_sum=mpz(0)
    for i in range(m):
        zs_1[i]=z
        zs_m[i]=zs_m_i
        ys_m[i]=ys_m_i
        ys_m_sum = (ys_m_sum + ys_m_i) %G.q
        zs_m_i= zs_m_i *z % G.q
        ys_m_i= ys_m_i *y %G.q
    """zs_1_l=np.broadcast_to(zs_1,(l,m))
    ys_m_l=np.broadcast_to(ys_m,(l,m))
    zs_m_l=np.broadcast_to(zs_m,(l,m)) """
    #s_Rs_ys_m=np.zeros((l,m),dtype=object)#np.multiply(ys_m_l,s_Rs)%G.q

    #t1s=(np.einsum('ij,ij->i',(a_Ls-zs_1_l)%G.q,s_Rs_ys_m)+np.einsum('ij,ij->i',s_Ls,zs_m_l)%G.q+np.einsum('ij,ij->i',s_Ls,np.multiply(ys_m_l,(a_Rs+zs_1_l)%G.q)%G.q)%G.q)%G.q
    #t2s=np.einsum('ij,ij->i',s_Ls,s_Rs_ys_m)%G.q
    t1s=np.zeros((l),dtype=object)
    t2s=np.zeros((l),dtype=object)
    for j in range(l):
        sum_t1=mpz(0)
        sum_t2=mpz(0)
        for i in range(m):
            s_Rs_ys_m=ys_m[i]*s_Rs[j][i]%G.q
            sum_t1= (sum_t1 +((mpz(a_Ls[j][i])-z)%G.q)*s_Rs_ys_m % G.q)%G.q
            sum_t1 = (sum_t1+s_Ls[j][i]*zs_m[i]%G.q+(((mpz(a_Rs[j][i])+z)%G.q*ys_m[i])%G.q)*s_Ls[j][i]%G.q)%G.q
            sum_t2= (sum_t2+ s_Ls[j][i]*s_Rs_ys_m%G.q)%G.q
        t1s[j]=sum_t1
        t2s[j]=sum_t2

    #step 5#
    #5.1
    tau1 = mpz(2)#G.random_exp()
    tau2 = mpz(2)#G.random_exp()
    #5.2
    T1=pow(h,tau1,G.p)
    T2=pow(h,tau2,G.p)
    
    for i in range(l):
        T1 = T1* pow(gs[i],mpz(t1s[i]),G.p) % G.p
        T2 = T2* pow(gs[i],mpz(t2s[i]),G.p) % G.p

    #step 6 - 7 -8#
    x=mpz(hash_elems(T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=G.q))
    #x = hashg(json.dumps({"T1": T1,"T2":T2,"A": A,"S":S,"gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list}),G.q)
    rand=0
    while x==0:
        x=mpz(hash_elems(T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=G.q))
        #x = hashg(json.dumps({"T1": T1,"T2":T2,"A": A,"S":S,"gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list,"rand":rand}),G.q)
        rand+=1
    #step 9#
    #9.1 - 9.2
    ls=np.zeros((l,m),dtype=object)
    rs=np.zeros((l,m),dtype=object)
    t_hats=np.zeros(l,dtype=object)
    x_2=pow(x,2,G.q)
    for j in range(l):        
        #ls[j,:]=(a_Ls[j]-zs_1)%G.q+s_Ls[j]*x %G.q%G.q
        rs_i=np.zeros(m,dtype=object)
        sum=mpz(0)
        for i in range(m):
            ls_i_i=((mpz(a_Ls[j][i])-z)%G.q+s_Ls[j][i]*x %G.q)%G.q
            ls[j][i]=ls_i_i
            rhs_i=((mpz(a_Rs[j][i])+z)%G.q+s_Rs[j][i]*x%G.q)%G.q
            rs_i_i= ((rhs_i*ys_m[i])%G.q + zs_m[i]) %G.q        
            rs_i[i]=rs_i_i
            sum = (sum +ls_i_i*rs_i_i %G.q )%G.q
        rs[j,:]=rs_i        
        t_hats[j]=sum#np.dot(ls[j],rs[j].T)%G.q
    #9.3
    #t_hats= [np.dot(ls[j],rs[j].T) % G.q for j in range(l)] #TODO try without transposing
    #9.4
    tau_x= tau2*x_2 %G.q + tau1*x %G.q 
    for i in range(m):
        tau_x= (tau_x+ zs_m[i]*gammas[i]) % G.q
    #9.5
    mu=(alpha+rho*x % G.q)%G.q
    #step 10 - 11 - 12#
    phi=mpz(hash_elems(tau_x,mu,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,q=G.q))
    #phi = hashg(json.dumps({"tau_x": tau_x, "mu":mu,"T1": T1,"T2":T2,"A": A,"S":S, "gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list}),G.q)
    rand=0
    while phi==0:
        phi=mpz(hash_elems(tau_x,mu,T1,T2,A,S,gs,hs, h,u,g_bolds,h_bolds, Vs,rand,q=G.q))
        #phi = hashg(json.dumps({"tau_x": tau_x, "mu":mu,"T1": T1,"T2":T2,"A": A,"S":S, "gs":gs_list, "hs" :hs_list, "h":h, "u":u, "g_bolds":g_bolds_list, "h_bolds":h_bolds_list, "Vs":Vs_list,"rand":rand}),G.q)
        rand+=1

    #step 13#
    """t_bar =0
    phis_l=np.zeros(l)
    phi_j=1
    for j in range(l):
        t_bar= (t_bar+t_hats[j]*phi_j%G.q) % G.q
        phis_l[j]=phi_j
        phi_j= phi_j *phi % G.q"""

    #step 15#
    #15.1 -15.4
    h_bolds_prime=np.zeros((l,m),dtype=object)
    g_bolds_prime=np.zeros((l,m),dtype=object)    
    phi_i=mpz(1)
    y_1=pow(y,G.q-2,G.q)
    phi_1=pow(phi,G.q-2,G.q)
    t_bar = mpz(0)
    phis_l=np.zeros(l,dtype=object)
    phi_j=mpz(1) 
    delta= ys_m_sum* (z-pow(z,2,G.q)) % G.q#np.sum(ys_m)%G.q* (z-pow(z,2,G.q)) % G.q
    P_bar = pow(h,G.q-tau_x,G.p) * pow(T1,x,G.p) % G.p * pow(T2,x_2,G.p) % G.p
    phis_l_1m=np.zeros((l,m),dtype=object)
    for i in range(m-1):
        delta = (delta - zs_m[i+1]) %G.q
        P_bar = P_bar * pow(Vs[i],zs_m[i],G.p) % G.p
    delta=(delta-zs_m[m-1]*z %G.q)%G.q
    P_bar = P_bar * pow(Vs[m-1],zs_m[m-1],G.p) % G.p
    P_ext=P_bar 
    P=A*pow(S,x,G.p) % G.p        
    for j in range(l):
        y_i=mpz(1)
        t_bar= (t_bar+t_hats[j]*phi_j%G.q) % G.q
        phis_l[j]=phi_j
        term= pow(gs[j],delta,G.p)
        P_bar=  P_bar*term %G.p
        P_ext= P_ext*term%G.p
        P_bar=  P_bar*pow(hs[j],phi_j,G.p) % G.p        
        for i in range(m):
            h_ij=pow(h_bolds[j][i],y_i,G.p)
            g_ij= pow(g_bolds[j][i],phi_i,G.p)
            h_bolds_prime[j][i]= h_ij
            g_bolds_prime[j][i]= g_ij
            P= P* pow(g_ij,mpz(G.q-z*phis_l[j] % G.q),G.p) % G.p
            P= P* pow(h_ij, mpz((z*ys_m[i]%G.q+zs_m[i]) %G.q),G.p) % G.p 
            phis_l_1m[j][i]=phi_j*ls[j][i]%G.q 
            y_i=y_i*y_1 % G.q
        phi_i=phi_i*phi_1 % G.q  
        phi_j= phi_j *phi % G.q        
            
    #step 16#
    #TODO : Bulletproof u?    
    #phis_l_1m = np.repeat(phis_l,m).reshape((l,m))
    a_1=phis_l_1m#np.multiply(np.array(ls),phis_l_1m)%G.q
    b_1=np.array(rs)
    g_bolds=np.array(g_bolds)
    h_bolds=np.array(h_bolds)
    t_hats=np.array(t_hats)
    bp1=bullet_proof(g_bolds_prime.flatten(),h_bolds_prime.flatten(),P*pow(h,G.q-mu,G.p)%G.p,u,a_1.flatten(),b_1.flatten(),t_bar,G.p,G.q)
    bp2=bullet_proof(gs,hs,P_bar,u,t_hats.flatten(),phis_l.flatten(),t_bar,G.p,G.q)
    #step 17#
    e1=extended_schnorr_proof(gs,P_ext,t_hats,G.p,G.q)  
    print("size",eg.compute_size(A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp1,bp2,e1,l,m))
    return A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp1,bp2,e1
    

#Verifier's function
def verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp_1,bp_2,e1):


    zs_m=np.zeros(m,dtype=object)
    zs_m_i= pow(z,2,G.q)
    ys_m=np.zeros(m,dtype=object)
    ys_m_i=mpz(1)
    x_2=pow(x,2,G.q)
    for i in range(m):
        zs_m[i]=zs_m_i
        ys_m[i]=ys_m_i
        zs_m_i=zs_m_i*z % G.q
        ys_m_i=ys_m_i*y %G.q


    #step 15#
    #15.1 -15.2
    h_bolds_prime=np.zeros((l,m),dtype=object)
    g_bolds_prime=np.zeros((l,m),dtype=object)
    y_1=pow(y,G.q-2,G.q)
    phi_i=mpz(1)
    phi_1=pow(phi,G.q-2,G.q)
    phis_l=np.zeros(l,dtype=object)
    phi_i_l=mpz(1)
    for j in range(l):
        phis_l[j]=phi_i_l
        phi_i_l= phi_i_l*phi %G.q
        y_i=mpz(1)
        for i in range(m):
            h_bolds_prime[j][i]= pow(h_bolds[j][i],y_i,G.p)
            g_bolds_prime[j][i]= pow(g_bolds[j][i],phi_i,G.p)
            y_i=y_i*y_1 % G.q
        phi_i=phi_i*phi_1 % G.q
    #15.3 -15.4
    delt= delta(y,z)
    #TODO verify how to compute/receive A, S, tau_x, T1, T2 
    P=A*pow(S,x,G.p) % G.p
    P_bar = pow(h,G.q-tau_x,G.p) * pow(T1,x,G.p) % G.p * pow(T2,x_2,G.p) % G.p
    P_bar_e1=P_bar
    #TODO fusionner avec boucles précédentes !
    for i in range(m):
        P_bar = P_bar*pow(Vs[i],zs_m[i],G.p) % G.p
        P_bar_e1=P_bar
    for j in range(l):    
        term=pow(gs[j],delt,G.p)     
        P_bar= P_bar*term% G.p
        P_bar_e1=P_bar_e1*term %G.p
        P_bar= P_bar*pow(hs[j],phis_l[j],G.p) % G.p
        for i in range(m):
            P= P*pow(g_bolds_prime[j][i], G.q-z*phis_l[j] % G.q, G.p) %G.p
            P= P*pow(h_bolds_prime[j][i], (z*ys_m[i]%G.q+zs_m[i]) % G.q,G.p)%G.p    
    
    #step 16#
    (a_bp1,b_bp1,Ls_bp1,Rs_bp1,xs_bp1),x_bp1=bp_1
    if not bullet_verification(g_bolds_prime.flatten(),h_bolds_prime.flatten(),
                               P*pow(h,G.q-mu,G.p)%G.p,u,a_bp1,b_bp1,t_bar,xs_bp1,x_bp1,Ls_bp1,Rs_bp1,G.p,G.q):
        return 1
    
    (a_bp2,b_bp2,Ls_bp2,Rs_bp2,xs_bp2),x_bp2=bp_2
    if not  bullet_verification(gs,hs,P_bar,u,a_bp2,b_bp2,t_bar,xs_bp2,x_bp2,Ls_bp2,Rs_bp2,G.p,G.q):
        return 2    

    #step 17#
    a_e1,xs_e1,Ls_e1,Rs_e1=e1
    if not extended_schnorr_verification(gs,P_bar_e1,a_e1,xs_e1,Ls_e1,Rs_e1,G.p,G.q):
        return 3
    return 0


n_false=0
print("generation time:", t00-t0)
for i in range(1):
    t1=time.time()
    A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp1,bp2,e1=generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas)
    t2=time.time()
    print("prover time:", t2-t1)
    t=verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp1,bp2,e1)
    t3=time.time()
    print("verifier time:", t3-t2)
    if t!=0:
        n_false+=1
        print(t)
print(n_false)
print("--- %s seconds ---" % (time.time() - start_time))