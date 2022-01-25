#%% Cell 1
import numpy as np

#%%

class option:
    def __init__(self,strike,maturity):
        self.K = strike
        self.T = maturity
class put(option):
    def payoff(self,ST):
        return [max(0,self.K-i) for i in ST]
#%% 
class call(option):
    def payoff(self,ST):
        return [max(0,i-self.K) for i in ST]

#%%
ST=np.array([100.0001,99])
a=call(100,5)
print(a.payoff(ST))
#%%
def compute_coeffs(r,sigma,T,N,i,theta,delta):
    alpha = -sigma**2/2
    beta = (sigma**2)/2-r
    gamma = r
    dx=1/N
    dt=1/T
    if i!=0 and i!=N:
        omega = theta*dt*(2*alpha/dx**2 -gamma)-1
        a = (1-theta)*dt*(-2*alpha/(dx)**2+gamma)-1
        b=theta*dt*(alpha/dx**2-beta/(2*dx))
        c = theta*dt*(alpha/dx**2+beta/(2*dx))
        d=(1-theta)*dt*(alpha/dx**2+beta/(2*dx))
        e=(1-theta)*dt*((alpha/dx**2-beta/(2*dx)))    
    if i==0 :
        omega=-theta*dt*(alpha/dx**2-beta/dx+gamma)-1
        a=(1-theta)*dt*(alpha/dx**1-beta/dx+gamma)-1
        b=theta*dt*(beta/dx-2*alpha/dx**2)
        c=(1-theta)*dt*(beta/dx-2*alpha/dx**2)
        d=theta*dt*alpha/dx**2
        e=(1-theta)*dt*alpha/dx**2
    else :
        omega=-theta*dt*(alpha/dx**2+beta/dx+gamma)-1
        a=(1-theta)*dt*(alpha/dx**1+beta/dx+gamma)-1
        b=-theta*dt*(2*alpha/dx**2+beta/dx)
        c=(1-theta)*dt*(beta/dx-2*alpha/dx**2)
        c=-(1-theta)*dt*(2*alpha/dx**2+beta/dx)
        d=theta*dt*alpha/dx**2
        e=(1-theta)*dt*alpha/dx**2
    return(omega,a,b,c,d,e)
        
#%%
def matrices(r,sigma,T,N,theta,delta):
    def A_prime(r,sigma,T,N,theta,delta):
        A_p = np.zeros((N+1,N+1))
        A_p[0][:3]=[compute_coeffs(r, sigma, T, N, N, theta, delta)[5-2*i] for i in range(0,3)]
        c1=compute_coeffs(r, sigma, T, N, 2, theta, delta)[5] #on met un indice random car coeff identique
        c2=compute_coeffs(r, sigma, T, N, 2, theta, delta)[1]
        c3=compute_coeffs(r, sigma, T, N, 2, theta, delta)[4]
        for i in range(1,N):
                A_p[i][i-1]=c1
                A_p[i][i]=c2
                A_p[i][i+1]=c3
        A_p[N,N]=compute_coeffs(r, sigma, T, N, 0, theta, delta)[5]
        A_p[N,N-1]=compute_coeffs(r, sigma, T, N, 0, theta, delta)[3]
        A_p[N,N-2]=compute_coeffs(r, sigma, T, N, 0, theta, delta)[2]
        return A_p
    def A_pp(r,sigma,T,N,theta,delta):
        A_pp = np.zeros((N+1,N+1))
        A_pp[0][:2]=[compute_coeffs(r, sigma, T, N, N, theta, delta)[2*i] for i in range(2,0,-1)]
        for i in range(1,N):
                A_pp[i][i-1]=compute_coeffs(r, sigma, T, N, N-i, theta, delta)[2]
                A_pp[i-1][i+1]=compute_coeffs(r, sigma, T, N, N-i, theta, delta)[3]
        A_pp[N,N]=compute_coeffs(r, sigma, T, N, 0, theta, delta)[5]
        A_pp[N,N-1]=compute_coeffs(r, sigma, T, N, 0, theta, delta)[3]
        A_pp[N,N-2]=compute_coeffs(r, sigma, T, N, 0, theta, delta)[2]
        return A_pp
    def Omega(r,sigma,T,N,theta,delta):
        L_om= np.full(N+1,compute_coeffs(r, sigma, T, N, 1, theta, delta)) #on met un indice random car coeff identique
        Om = np.diag(L_om)
        Om[0][0] = compute_coeffs(r, sigma, T, N, N, theta, delta)
        Om[N][N] = compute_coeffs(r, sigma, T, N, 0, theta, delta)
        return Om
    
#%%
L=np.full(4,0.001)
print(L)