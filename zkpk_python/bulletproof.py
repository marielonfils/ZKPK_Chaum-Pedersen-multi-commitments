import numpy as np
import json
from utils.hash import hashg



def bullet_proof(gs,hs,P,u,a_s,bs,c,p,q):
    Ls=[]
    Rs=[]
    xs=[]
    def proof(gs,hs,P,u,a_s,bs,n):
        if n== 1:
            return a_s,bs,Ls,Rs,xs
        
        #step 1#
        n_prime=int(n/2)
        a_l=a_s[:n_prime]
        a_r=a_s[n_prime:]
        b_l=bs[:n_prime]
        b_r=bs[n_prime:]
        g_l=gs[:n_prime]
        g_r=gs[n_prime:]
        h_l=hs[:n_prime]
        h_r=hs[n_prime:]

        c_L=int(np.dot(a_l,b_r)%q)
        c_R=int(np.dot(a_r,b_l)%q)
        L=pow(u,c_L,p)
        R=pow(u,c_R,p)
        for i in range(n_prime):
            L= L*pow(int(g_r[i]),int(a_l[i]),p)%p
            L= L*pow(int(h_l[i]),int(b_r[i]),p)%p
            R= R*pow(int(g_l[i]),int(a_r[i]),p)%p
            R= R*pow(int(h_r[i]),int(b_l[i]),p)%p
        Ls.append(L)
        Rs.append(R)

        #step 2-3-4 #
        x=hashg(json.dumps({"L": L,"R":R,"g" :gs,"h":hs,"a":a_s.astype(int).tolist(),"b":bs.astype(int).tolist(),"c":c,"u":u,"P":P,"p":p,"q":q}),q)
        
        rand=0
        while x==0:
            x=hashg(json.dumps({"L": L,"R":R,"g" :gs,"h":hs,"a":a_s.astype(int).tolist(),"b":bs.astype(int).tolist(),"c":c,"u":u,"P":P,"p":p,"q":q, "rand":rand}),q)
            rand+=1
        x_1= pow(x, q-2,q)
        xs.append(x)

        #step 5#
        g_prime=[]
        h_prime=[]
        for i in range(n_prime):
            g_prime.append((pow(int(g_r[i]),x,p)*pow(int(g_l[i]),x_1,p))% p)
            h_prime.append((pow(int(h_r[i]),x_1,p)*pow(int(h_l[i]),x,p))% p)
        x2=pow(x,2,q)
        x_2=pow(x,q-3,q)
        P_prime=(pow(L,x2,p)*P)%p * pow(R,x_2,p) %p
        a_prime=((a_r*(x_1))%q+ (a_l*x)%q)%q 
        b_prime=((b_r*x)%q+(b_l*(x_1))%q)%q  


        return proof(g_prime,h_prime,P_prime,u,a_prime,b_prime,n_prime)#!

    if len(gs)%2 !=0:
        ValueError()
    x=hashg(json.dumps({"g" :gs,"h":hs,"a":a_s.astype(int).tolist(),"b":bs.astype(int).tolist(),"c":c,"u":u,"P":P,"p":p,"q":q}),q)
    rand=0
    while x==0:
        x=hashg(json.dumps({"g" :gs,"h":hs,"a":a_s.astype(int).tolist(),"b":bs.astype(int).tolist(),"c":c,"u":u,"P":P,"p":p,"q":q,"rand":rand}),q)
        rand+=1
    u_xc=pow(u,x*c%q,p)
    u_x=pow(u,x,p)
    return proof(gs,hs,P*u_xc%p,u_x,a_s,bs,len(gs)),x
    #return proof(gs,hs,P,u,a_s,bs,len(gs)),1

def bullet_verification(gs,hs,P,u,a,b,c,xs,x,Ls,Rs,p,q):
    n_prime=len(gs)
    g_prime=gs
    h_prime=hs
    P_prime=P*pow(u,x*c%q,p)%p
    u_prime=pow(u,x,p)
    i=0
    while(n_prime >1):
        #step 1#
        n_prime=int(n_prime/2)

        #step 5#
        g_l=g_prime[:n_prime]
        g_r=g_prime[n_prime:]
        h_l=h_prime[:n_prime]
        h_r=h_prime[n_prime:]

        g_prime=[]
        h_prime=[]
        x=xs[i]
        x_1=pow(x,q-2,q)
        for j in range(n_prime):            
            g_prime.append((pow(g_r[j],x,p)*pow(g_l[j],x_1,p))% p)
            h_prime.append((pow(h_r[j],x_1,p)*pow(h_l[j],x,p))% p)
        
        x2=pow(x,2,q)
        x_2=pow(x,q-3,q)
        P_prime=(pow(Ls[i],x2,p)*P_prime)%p * pow(Rs[i],x_2,p) %p
        i+=1

    c_prime=int(a[0]*b[0] %q)
    u_c=pow(u_prime,c_prime,p)
    if P_prime==((pow(g_prime[0],int(a[0]),p) * pow(h_prime[0],int(b[0]),p) % p )*u_c %p):
        return True
    return False

def test():
    g=[3,4]
    h=[3,4]
    P=1 #3^2%23=9, 4^3%23=18, 9*18%23=1
    a=[2,3]
    b=[2,3] #<a,b> = 13=2mod 11
    c=2
    u = 2  #u^(<a,b>)=4
    p=23
    q=11
    #x=2
    (a,b,Ls,Rs,xs),x=bullet_proof(g,h,P,u,np.array(a),np.array(b),c,p,q)
    assert bullet_verification(g,h,P,u,a,b,c,xs,x,Ls,Rs,p,q)
#test()
