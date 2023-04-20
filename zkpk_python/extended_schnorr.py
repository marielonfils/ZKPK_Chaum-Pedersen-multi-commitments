import numpy as np
import json
from utils.hash import hashg

xs=[]
Ls=[]
Rs=[]

def extended_schnorr_proof(g,P,a,p,q):
    def proof(g,P,a,n):

        if n==1:
            return (a,xs,Ls,Rs)   

        #step 1#
        n_prime=int(n/2)
        a_l=a[:n_prime]
        a_r=a[n_prime:]
        g_l=g[:n_prime]
        g_r=g[n_prime:]

        L=1
        R=1
        for i in range(n_prime):
            L*=pow(g_r[i],int(a_l[i]),p)%p
            R*=pow(g_l[i],int(a_r[i]),p)%p
        Ls.append(L)
        Rs.append(R)

        #step 2-3-4 #
        x = hashg(json.dumps({"L": L,"R":R, "g" :g, "P":P}),q)
        x_1= pow(x, q-1,p)
        xs.append(x)

        #step 5#
        g_prime=[]
        for i in range(n_prime):
            g_prime.append((pow(g_r[i],x,p)*pow(g_l[i],x_1,p))% p)
        x2=pow(x,2,p)
        x_2=pow(x,-2,p)
        P_prime=(pow(L,x2,p)*P)%p * pow(R,x_2,p) %p
        a_prime=((a_r*x_1)%q+ (a_l*x)%q)%q   
        return proof(g_prime,P_prime,a_prime,n_prime)
    if len(g)%2 !=0:
        ValueError()
    return proof(g,P,a,len(g))


def extended_schnorr_verification(g,P,a,xs,Ls,Rs,p,q):
    n_prime=len(g)
    a_prime=a
    g_prime=g
    P_prime=p
    i=0
    while(n_prime >1):
        #step 1#
        n_prime=int(n_prime/2)

        #step 5#
        g_l=g_prime[:n_prime]
        g_r=g_prime[n_prime:]

        g_prime=[]
        x=xs[i]
        x_1=pow(x,q-1,p)
        for j in range(n_prime):            
            g_prime.append((pow(g_r[j],x,p)*pow(g_l[j],x_1,p))% p)
        
        x2=pow(x,2,p)
        x_2=pow(x,-2,p)
        P_prime=(pow(Ls[i],x2,p)*P_prime)%p * pow(Rs[i],x_2,p) %p
        i+=1
    if P==pow(g_prime[0],int(a_prime[0]),p):
        return True
    return False

def test():
    g=[3,5]
    P=2 #3^2%7=2, 5^3%7=1, 2*1%5=2
    a=[2,3]
    p=7
    q=6
    a,xs,Ls,Rs=extended_schnorr_proof(g,P,np.array(a),p,q)
    assert extended_schnorr_verification(g,P,a,xs,Ls,Rs,p,q)
test()