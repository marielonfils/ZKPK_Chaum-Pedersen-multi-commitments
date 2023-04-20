import numpy as np
import utils.elgamal as eg
from utils.hash import hashg 
import json
#change QNbits in elgamal.py
#TODO modulo in for *= %G.p

# groups size
l=2
m=3
p=7

# generators
G=eg.gen_group()
h=eg.random_generator(G.p,G.q)
u=eg.random_generator(G.p, G.q)
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
gammas=[]
for i in range(m):
    gammas_i=G.random_exp()
    gammas.append(gammas_i)
    product=pow(h,gammas_i,G.p)
    for j in range(l):
        product*=pow(gs[j],vs[j][i])
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
        y_m+=y_m_i
        y_m_i*=y % G.q
        z_m+=z_m_i
        z_m_i*=z %G.q
    return ((y_m*(z-z_2)%G.q )%G.q- z_m)  %G.q

#Prover's function
def generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas):
    #step 1#
    a_Ls=vs.copy()
    #verify if modulo 2 TODO
    a_Rs=a_Ls-np.ones((l,m))
    alpha=G.random_exp()
    rho=G.random_exp()
    s_Ls=[]
    s_Rs=[]
    for _ in range(l):
        s_Ls_i=[]
        s_Rs_i=[]
        for _ in range(m):
            s_Ls_i.append(G.random_exp())
            s_Rs_i.append(G.random_exp())
        s_Ls.append(s_Ls_i)
        s_Rs.append(s_Rs_i)
    A=pow(h,alpha,G.p)
    S=pow(h,rho,G.p)
    for i in range(l):
        for j in range(m):
            A*=pow(g_bolds[i][j],a_Ls[i][j].item(),G.p) % G.p
            A*=pow(h_bolds[i][j],int(a_Rs[i][j].item()),G.p) % G.p
            S*=pow(g_bolds[i][j],s_Ls[i][j],G.p) % G.p
            S*=pow(h_bolds[i][j],s_Rs[i][j],G.p) % G.p
    
    #step 2 - 3 - 4#
    y = hashg(json.dumps({"A": A,"gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs}),G.q)
    z = hashg(json.dumps({"S": S,"gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs}),G.q)
    #TODO : verify hashes ?
    

    #step 4#
    #TODO : verify that we have nothing todo here

    #step 5#
    tau1 = G.random_exp()
    tau2 = G.random_exp()
    T1=pow(h,tau1,G.p)
    T2=pow(h,tau2,G.p)
    #TODO: compute (t0s),t1s and t2s
    #t0s=np.ones(l)
    #TODO check that t1s and t2s contain int for pow function below
    t1s=np.ones(l)
    t2s=np.ones(l)
    for i in range(l):
        T1*= pow(gs[i],int(t1s[i]),G.p) % G.p
        T2*= pow(gs[i],int(t2s[i]),G.p) % G.p

    #step 6 - 7 -8#
    #TODO : est-ce qu'il faut ajouter A,S,y,z,tau1,tau2,t0s,t1s,t2s dans le hash ?
    x = hashg(json.dumps({"T1": T1,"T2":T2,"gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs}),G.q)
    
    #step 9#
    ls=[]
    rs=[]
    zs_1=z*np.ones(m)
    zs_m=np.ones(m)
    zs_m_i= pow(z,2,G.q)
    #TODO : \mathbf{y}^m c'est bien un vecteur (1,y,...y^m) ?
    ys_m=np.ones(m)
    ys_m_i=1
    x_2=pow(x,2,G.q)
    for i in range(m):
        zs_m[i]=zs_m_i
        ys_m[i]=ys_m_i
        zs_m_i*=z % G.q
        ys_m_i*=y %G.q
    for j in range(l):
        ls.append(np.array(a_Ls[j]-zs_1+np.array(s_Ls[j])*x %G.q)%G.q)
        rs_i=[]
        rhs=np.array(a_Rs[j]+zs_1+np.array(s_Rs[j])*x%G.q)%G.q
        #TODO verify this expression below
        for i in range(m):            
            rs_i.append(((rhs[i]*ys_m[i])%G.q + zs_m[i]) %G.q)
            zs_m*=z%G.q
        rs.append(rs_i)
    t_hats= [np.dot(np.array(ls[j]),np.array(rs[j]).T) % G.q for j in range(l)]
    tau_x= tau2*x_2 %G.q + tau1*x %G.q 
    for j in range(1,m-1):
        tau_x+= zs_m[j]*gammas[j-1] % G.q
    tau_x+= (zs_m[m-1]*z) %G.q *gammas[m-1] % G.q
    tau_x=int(tau_x)
    mu=alpha+rho*x % G.q

    #step 10 - 11 - 12#
    #TODO : verify if we need to add other things in the hash
    phi = hashg(json.dumps({"tau_x": tau_x, "mu":mu, "gs":gs, "hs" :hs, "h":h, "u":u, "g_bolds":g_bolds, "h_bolds":h_bolds, "Vs":Vs}),G.q)

    #step 13#
    t_bar =0
    phi_j=1
    for j in range(l):
        t_bar+=t_hats[j]*phi_j % G.q
        phi_j*=phi
    
    #step 14#
    #TODO send t_bar to V!!

    #step 15#
    h_bolds_prime=np.zeros((l,m))
    g_bolds_prime=np.zeros((l,m))
    y_i=1
    phi_i=1
    phis_l=np.zeros(l)
    phi_i_l=1
    y_1=pow(y,G.q-1,G.p)
    phi_1=pow(phi,G.q-1,G.p)
    for j in range(l):
        phis_l[j]=phi_i_l
        phi_i_l*=phi %G.q
        for i in range(m):
            # TODO : more efficient inverse computation
            h_bolds_prime[j][i]= pow(h_bolds[j][i],y_i,G.p)
            g_bolds_prime[j][i]= pow(g_bolds[j][i],phi_i,G.p)
            y_i*=y_1 % G.q
        phi_i*=phi_1 % G.q

    delta= np.sum(ys_m)* (z-pow(z,2,G.q)) % G.q 
    den=pow(h,tau_x,G.p)
    P_bar = pow(den,G.q-1,G.p) * pow(T1,x,G.p) % G.p * pow(T2,x_2,G.p) * G.p
    for i in range(m-1):
        delta-= zs_m[i] %G.q
        P_bar *= pow(Vs[i],int(zs_m[i]),G.p) % G.p
    delta-=zs_m[m-1]*z %G.q

    P=A*S % G.p    
    for j in range(l):
        P_bar*= pow(hs[j],int(phis_l[j]),G.p) % G.p
        P_bar*= pow(gs[j],int(delta),G.p) %G.p
        for i in range(m):
            den = pow(int(g_bolds_prime[j][i]),int(z*phis_l[j] % G.q), G.p )
            P*= pow(den,G.q-1,G.p) % G.p #%G.p#pow(g_bolds_prime[j][i], -z*phis_l[j] % G.q, G.p) %G.p
            #TODO verify the expression below
            P*= pow(int(h_bolds_prime[j][i]), int((z*ys_m[i]+zs_m[i]) %G.q),G.p) % G.p
            P=1    
    
    #step 16#
    #TODO : Bulletproof

    #step 17#
    #TODO : extended schnorr

    return 1

#Verifier's function
def verify_proof(gs,hs,hu,g_bolds,h_bolds,Vs,y,z,x,phi,t_bar):

    zs_m=np.ones(m)
    zs_m_i= pow(z,2,G.q)
    ys_m=np.ones(m)
    ys_m_i=1
    x_2=pow(x,2,G.q)
    for i in range(m):
        zs_m[i]=zs_m
        ys_m[i]=ys_m_i
        zs_m_i*=z % G.q
        ys_m_i*=y %G.q

    #step 14#
    #TODO send t_bar to V!!

    #step 15#
    h_bolds_prime=[]
    g_bolds_prime=[]
    y_i=1
    phi_i=1
    phis_l=[]
    phi_i_l=1
    for j in range(l):
        phis_l[j]=phi_i_l
        phi_i_l*=phi %G.q
        for i in range(m):
            h_bolds_prime[j][i]= pow(h_bolds[j][i],y_i,G.p)
            g_bolds_prime[j][i]= pow(g_bolds[j][i],phi_i,G.p)
            y_i/=y % G.q
        phi_i/=phi % G.q
    
    #TODO compute delta
    delt= delta(y,z)
    #TODO verify how to compute/receive A, S, tau_x, T1, T2 
    #P=A*S % G.p
    #P_bar = pow(h,-tau_x,G.p) * pow(T1,x,G.p) % G.p * pow(T2,x_2,G.p) * G.p
    for j in range(l):
        P_bar*= pow(h[j],phis_l[j]) % G.p
        P_bar*= pow(g[j],delt,G.p)
        for i in range(m):
            P*= pow(g_bolds_prime[j][i], -z*phis_l[j] % G.q, G.p) %G.p
            #TODO verify the expression below
            P*= pow(h_bolds_prime[j][i], z*ys_m[i]+zs_m[i]) % G.p
    
    for i in range(m):
        P_bar *= pow(Vs[i],zs_m[i],G.p) % G.p
    
    #step 16#
    #TODO : Bulletproof

    #step 17#
    #TODO : extended schnorr
    return 0

generate_proof(gs,hs,h,u,g_bolds,h_bolds,Vs,vs,gammas)