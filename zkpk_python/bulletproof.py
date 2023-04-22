import numpy as np
import json
from utils.hash import hashg

Ls=[]
Rs=[]
xs=[]

def bullet_proof(gs,hs,u,P,a_s,bs,p,q):

    def proof(gs,hs,P,a_s,bs,n):
        if n== 1:
            return (a_s,bs,Ls,Rs,xs)
        
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
            L= L*pow(g_r[i],int(a_l[i]),p)%p
            L= L*pow(h_l[i],int(b_r[i]),p)%p
            R= R*pow(g_l[i],int(a_r[i]),p)%p
            R= R*pow(h_r[i],int(b_l[i]),p)%p
        Ls.append(L)
        Rs.append(R)

        #step 2-3-4 #
        x = hashg(json.dumps({"L": L,"R":R, "g" :gs, "P":P}),q)
        rand=0
        while x==0:
            x = hashg(json.dumps({"L": L,"R":R, "g" :gs, "P":P, "rand":rand}),q)
            rand+=1
        x_1= pow(x, q-2,q)
        xs.append(x)

        #step 5#
        g_prime=[]
        h_prime=[]
        for i in range(n_prime):
            g_prime.append((pow(g_r[i],x,p)*pow(g_l[i],x_1,p))% p)
            h_prime.append((pow(h_r[i],x_1,p)*pow(h_l[i],x,p))% p)
        x2=pow(x,2,q)
        x_2=pow(x,q-3,q)
        P_prime=(pow(L,x2,p)*P)%p * pow(R,x_2,p) %p
        a_prime=((a_r*(x_1))%q+ (a_l*x)%q)%q 
        b_prime=((b_r*x)%q+(b_l*(x_1))%q)%q  

        return proof(g_prime,h_prime,P_prime,a_prime,b_prime,n_prime)

    if len(gs)%2 !=0:
        ValueError()
    return proof(gs,hs,P,a_s,bs,len(gs))

def bullet_verification(gs,hs,u,P,a,b,xs,Ls,Rs,p,q):
    n_prime=len(gs)
    g_prime=gs
    h_prime=hs
    P_prime=P
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

    c=int(a[0]*b[0] %q)
    u_c=pow(u,c,p)
    if P_prime==((pow(g_prime[0],int(a[0]),p) * pow(h_prime[0],int(b[0]),p) % p )*u_c %p):
        return True
    return False

def test():
    g=[3,4]
    h=[3,4]
    P=4 #3^2%23=9, 4^3%23=18, 9*18%23=1
    a=[2,3]
    b=[2,3] #<a,b> = 13=2mod 11
    u = 2  #u^(<a,b>)=4
    p=23
    q=11
    a,b,Ls,Rs,xs=bullet_proof(g,h,u,P,np.array(a),np.array(b),p,q)
    print(a,b,xs,Ls,Rs)
    assert bullet_verification(g,h,u,P,a,b,xs,Ls,Rs,p,q)
test()
