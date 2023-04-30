import numpy as np
import utils.elgamal as eg
from utils.hash import hashg 
import json
from bulletproof import bullet_proof, bullet_verification
from extended_schnorr import extended_schnorr_proof, extended_schnorr_verification
#change QNbits in elgamal.py


# groups size
l=4
m=4
#p=7

# generators
G=eg.gen_group()
h=eg.random_generator(G.p,G.q)
u=eg.random_generator(G.p,G.q)
gs=[]
hs=[]
g_bolds=[]
h_bolds=[]
for _ in range(l):
    gs.append(eg.random_generator(G.p, G.q))
    hs.append(eg.random_generator(G.p, G.q))
    g_bolds_i=[]
    h_bolds_i=[]
    for _ in range(m):
        g_bolds_i.append(eg.random_generator(G.p, G.q))
        h_bolds_i.append(eg.random_generator(G.p, G.q))
    g_bolds.append(g_bolds_i)
    h_bolds.append(h_bolds_i)

# votes and commitment on votes
Vs=[]
vs=np.random.randint(2,size=(l,m))
#TODO!!!

gammas=[]
for i in range(m):
    gammas_i=G.random_exp()
    gammas.append(gammas_i)
    product=pow(h,gammas_i,G.p)
    for j in range(l):
        vs[j,i]=(i+j)%2
        product= product *pow(gs[j],int(vs[j][i]),G.p)%G.p
    Vs.append(int(product))


def delta(y,z):
    # delta(y,z) := (z-z^2) <1^m,y^m> - \sum_{k=1}^m z^{k+2}
    #TODO :verify if y^m correct
    y_m=0
    y_m_i=1
    z_2=pow(z,2,G.q)
    z_m_i=z_2*z % G.q
    z_m=0
    for i in range(m):
        y_m=(y_m +y_m_i)%G.q
        y_m_i*=y % G.q
        z_m= (z_m+z_m_i)%G.q
        z_m_i*=z %G.q
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
    s_Ls=np.zeros((l,m))
    s_Rs=np.zeros((l,m))
    for j in range(l):
        for i in range(m):
            s_Ls[j][i]=G.random_exp()
            s_Rs[j][i]=G.random_exp()
    #1.4 - 1.5
    A=pow(h,alpha,G.p)
    S=pow(h,rho,G.p)
    for i in range(l):
        for j in range(m):
            A=A*pow(g_bolds[i][j],a_Ls[i][j].item(),G.p) % G.p
            A=A*pow(h_bolds[i][j],int(a_Rs[i][j].item()),G.p) % G.p
            S=S*pow(g_bolds[i][j],int(s_Ls[i][j]),G.p) % G.p
            S=S*pow(h_bolds[i][j],int(s_Rs[i][j]),G.p) % G.p
    
    #step 2 - 3 - 4#
    y = hashg(json.dumps({"A": A,"S":S, "gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs}),G.q)
    rand=0
    while y==0:
        y = hashg(json.dumps({"A": A,"S":S, "gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs, "rand":rand}),G.q)
        rand+=1
    z = hashg(json.dumps({"A": A,"S": S,"gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs, "r":42}),G.q)
    rand =0
    while z==0:
        z = hashg(json.dumps({"A": A,"S": S,"gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs, "r":42,"rand":rand}),G.q)
    #step 4
    zs_1=z*np.ones(m)
    zs_m=np.ones(m)
    zs_m_i= pow(z,2,G.q)
    ys_m=np.ones(m)
    ys_m_i=1    
    for i in range(m):
        zs_m[i]=zs_m_i
        ys_m[i]=ys_m_i
        zs_m_i= zs_m_i *z % G.q
        ys_m_i= ys_m_i *y %G.q
    zs_1_l=np.broadcast_to(zs_1,(l,m))
    ys_m_l=np.broadcast_to(ys_m,(l,m))
    zs_m_l=np.broadcast_to(zs_m,(l,m)) 
    s_Rs_ys_m=np.multiply(ys_m_l,s_Rs)%G.q

    t1s=(np.einsum('ij,ij->i',(a_Ls-zs_1_l)%G.q,s_Rs_ys_m)+np.einsum('ij,ij->i',s_Ls,zs_m_l)%G.q+np.einsum('ij,ij->i',s_Ls,np.multiply(ys_m_l,(a_Rs+zs_1_l)%G.q)%G.q)%G.q)%G.q
    t2s=np.einsum('ij,ij->i',s_Ls,s_Rs_ys_m)%G.q

    #step 5#
    #5.1
    tau1 = G.random_exp()
    tau2 = G.random_exp()
    #5.2
    T1=pow(h,tau1,G.p)
    T2=pow(h,tau2,G.p)
    
    for i in range(l):
        T1 = T1* pow(gs[i],int(t1s[i]),G.p) % G.p
        T2 = T2* pow(gs[i],int(t2s[i]),G.p) % G.p

    #step 6 - 7 -8#
    x = hashg(json.dumps({"T1": T1,"T2":T2,"A": A,"S":S,"gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs}),G.q)
    rand=0
    while x==0:
        x = hashg(json.dumps({"T1": T1,"T2":T2,"A": A,"S":S,"gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs,"rand":rand}),G.q)
        rand+=1
    #step 9#
    #9.1 - 9.2
    ls=[]
    rs=[]
    x_2=pow(x,2,G.q)
    for j in range(l):
        ls.append((a_Ls[j]-zs_1)%G.q+s_Ls[j]*x %G.q%G.q)
        rs_i=[]
        rhs=np.array((a_Rs[j]+zs_1)%G.q+s_Rs[j]*x%G.q)%G.q
        #TODO verify this expression below
        for i in range(m):            
            rs_i.append(((rhs[i]*ys_m[i])%G.q + zs_m[i]) %G.q)
            #zs_m*=z%G.q
        rs.append(rs_i)
    #9.3
    t_hats= [np.dot(np.array(ls[j]),np.array(rs[j]).T) % G.q for j in range(l)] #TODO try without transposing
    #9.4
    tau_x= tau2*x_2 %G.q + tau1*x %G.q 
    for j in range(m):
        tau_x= (tau_x+ zs_m[j]*gammas[j]) % G.q
    #tau_x+= (zs_m[m-1]*z) %G.q *gammas[m-1] % G.q
    tau_x=int(tau_x)
    #9.5
    mu=alpha+rho*x % G.q

    #step 10 - 11 - 12#
    phi = hashg(json.dumps({"tau_x": tau_x, "mu":mu,"T1": T1,"T2":T2,"A": A,"S":S, "gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs}),G.q)
    rand=0
    while phi==0:
        phi = hashg(json.dumps({"tau_x": tau_x, "mu":mu,"T1": T1,"T2":T2,"A": A,"S":S, "gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs,"rand":rand}),G.q)
        rand+=1
    #step 13#
    t_bar =0
    phis_l=np.zeros(l)
    phi_j=1
    for j in range(l):
        t_bar= (t_bar+t_hats[j]*phi_j%G.q) % G.q
        phis_l[j]=phi_j
        phi_j= phi_j *phi % G.q

    #step 15#
    #15.1 -15.2
    h_bolds_prime=np.zeros((l,m))
    g_bolds_prime=np.zeros((l,m))    
    phi_i=1
    y_1=pow(y,G.q-2,G.q)
    phi_1=pow(phi,G.q-2,G.q)
    for j in range(l):
        y_i=1
        for i in range(m):
            # TODO : more efficient inverse computation
            h_bolds_prime[j][i]= pow(h_bolds[j][i],y_i,G.p)
            g_bolds_prime[j][i]= pow(g_bolds[j][i],phi_i,G.p)
            y_i=y_i*y_1 % G.q
        phi_i=phi_i*phi_1 % G.q
    #15.3 -15.4
    delta= np.sum(ys_m)%G.q* (z-pow(z,2,G.q)) % G.q 
    P_bar = pow(h,G.q-tau_x,G.p) * pow(T1,x,G.p) % G.p * pow(T2,x_2,G.p) % G.p
    for i in range(m-1):
        delta = (delta - zs_m[i+1]) %G.q
        P_bar = P_bar * pow(Vs[i],int(zs_m[i]),G.p) % G.p
    #TDO verify delta
    delta=(delta-zs_m[m-1]*z %G.q)%G.q
    P_bar = P_bar * pow(Vs[m-1],int(zs_m[m-1]),G.p) % G.p
    P_ext=P_bar
    P=A*pow(S,x,G.p) % G.p  
    for j in range(l):
        term= pow(gs[j],int(delta),G.p)
        P_bar=  P_bar*term %G.p
        P_ext= P_ext*term%G.p
        P_bar=  P_bar*pow(hs[j],int(phis_l[j]),G.p) % G.p     
        for i in range(m):
            P= P* pow(int(g_bolds_prime[j][i]),G.q-int(z*phis_l[j] % G.q),G.p) % G.p #%G.p#pow(g_bolds_prime[j][i], -z*phis_l[j] % G.q, G.p) %G.p
            P= P* pow(int(h_bolds_prime[j][i]), int((z*ys_m[i]%G.q+zs_m[i]) %G.q),G.p) % G.p  
    """for i in range(m):
        term=pow(Vs[i],int(zs_m[i]),G.p)
        P_bar = P_bar*term % G.p
        P_ext=P_ext*term % G.p
        print(term)"""
    #step 16#
    #TODO : Bulletproof u?, c?
    
    phis_l_1m = np.repeat(phis_l,m).reshape((l,m))
    a_1=np.multiply(np.array(ls),phis_l_1m)%G.q
    b_1=np.array(rs)
    c_1=a_1@b_1.T
    c_2=np.dot(t_hats,phis_l)
    g_bolds=np.array(g_bolds)
    h_bolds=np.array(h_bolds)
    t_hats=np.array(t_hats)
    bp1=bullet_proof(g_bolds_prime.flatten().astype(int).tolist(),h_bolds_prime.flatten().astype(int).tolist(),P*pow(h,G.q-mu,G.p)%G.p,u,a_1.flatten(),b_1.flatten(),int(t_bar),G.p,G.q)
    bp2=bullet_proof(gs,hs,P_bar,u,t_hats.flatten(),phis_l.flatten(),int(t_bar),G.p,G.q)

    #step 17#
    #TODO : extended schnorr
    e1=extended_schnorr_proof(gs,P_ext,t_hats,G.p,G.q)

    

    return A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp1,bp2,e1

#Verifier's function
def verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp_1,bp_2,e1):

    zs_m=np.zeros(m)
    zs_m_i= pow(z,2,G.q)
    ys_m=np.zeros(m)
    ys_m_i=1
    x_2=pow(x,2,G.q)
    for i in range(m):
        zs_m[i]=zs_m_i
        ys_m[i]=ys_m_i
        zs_m_i=zs_m_i*z % G.q
        ys_m_i=ys_m_i*y %G.q


    #step 15#
    #15.1 -15.2
    h_bolds_prime=np.zeros((l,m))
    g_bolds_prime=np.zeros((l,m))
    y_1=pow(y,G.q-2,G.q)
    phi_i=1
    phi_1=pow(phi,G.q-2,G.q)
    phis_l=np.zeros(l)
    phi_i_l=1
    for j in range(l):
        phis_l[j]=phi_i_l
        phi_i_l= phi_i_l*phi %G.q
        y_i=1
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
        P_bar = P_bar*pow(Vs[i],int(zs_m[i]),G.p) % G.p
        P_bar_e1=P_bar
    for j in range(l):    
        term=pow(gs[j],delt,G.p)     
        P_bar= P_bar*term% G.p
        P_bar_e1=P_bar_e1*term %G.p
        P_bar= P_bar*pow(hs[j],int(phis_l[j]),G.p) % G.p
        for i in range(m):
            P= P*pow(int(g_bolds_prime[j][i]), int(G.q-z*phis_l[j] % G.q), G.p) %G.p
            P= P*pow(int(h_bolds_prime[j][i]), int(z*ys_m[i]%G.q+zs_m[i] % G.q),G.p)%G.p    
    
    #step 16#
    (a_bp1,b_bp1,Ls_bp1,Rs_bp1,xs_bp1),x_bp1=bp_1
    if not bullet_verification(g_bolds_prime.flatten().astype(int).tolist(),h_bolds_prime.flatten().astype(int).tolist(),P*pow(h,G.q-mu,G.p)%G.p,u,a_bp1,b_bp1,int(t_bar),xs_bp1,x_bp1,Ls_bp1,Rs_bp1,G.p,G.q):
        return 1

    (a_bp2,b_bp2,Ls_bp2,Rs_bp2,xs_bp2),x_bp2=bp_2
    if not bullet_verification(gs,hs,P_bar,u,a_bp2,b_bp2,int(t_bar),xs_bp2,x_bp2,Ls_bp2,Rs_bp2,G.p,G.q):
        return 2
    

    #step 17#
    a_e1,xs_e1,Ls_e1,Rs_e1=e1
    if not extended_schnorr_verification(gs,P_bar_e1,a_e1,xs_e1,Ls_e1,Rs_e1,G.p,G.q):
        return 3
    return 0
n_false=0
for i in range(100):
    A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp1,bp2,e1=generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas)
    t=verify_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,A,S,y,z,T1,T2,x,tau_x,mu,phi,t_bar,bp1,bp2,e1)
    if t!=0:
        n_false+=1
        print(t)
print(n_false)