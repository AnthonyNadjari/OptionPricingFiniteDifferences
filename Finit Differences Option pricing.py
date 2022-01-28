#%% Cell 1
import numpy as np
import pandas as pd
#%%


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
      diag_low = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[5] for i in range(1,N)]
      diag_low.append(compute_coeffs(r, sigma, T, N, 0, theta, delta)[3])
      A_low = np.diag(diag_low,-1)

      diag = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[1] for i in range(N)]
      diag[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[5]
      diag.append(compute_coeffs(r, sigma, T, N, 0, theta, delta)[5])
      A_diag = np.diag(diag)
      
      diag_up = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[4] for i in range(0,N)]
      diag_up[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[3]
      A_up = np.diag(diag_up,1)
      A_p=A_low+A_diag+A_up

      return A_p
  def A_pp(r,sigma,T,N,theta,delta):
    A_pp = np.zeros((N+1,N+1))
    diag_low = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[2] for i in range(1,N+1)]
    
    diag=np.zeros(N-1)
    diag.append(compute_coeffs(r, sigma, T, N, 0, theta, delta)[4])
    diag[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[4]

    diag_up = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[3] for i in range(0,N)]
    diag_up[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[2]

    np.fill_diagonal(A_pp[1:,:-1],diag_low)
    np.fill_diagonal(A_pp,diag)
    np.fill_diagonal(A_pp[:-1,1:],diag_up)
  def Omega(r,sigma,T,N,theta,delta):
    L_om= np.full(N+1,compute_coeffs(r, sigma, T, N, 1, theta, delta)) #on met un indice random car coeff identique
    Om = np.diag(L_om)
    Om[0][0] = compute_coeffs(r, sigma, T, N, N, theta, delta)
    Om[N][N] = compute_coeffs(r, sigma, T, N, 0, theta, delta)
    return Om
  def system(r,sigma,T,N,theta,delta):
    A_p = A_prime(r,sigma,T,N,theta,delta)
    A_pprime = A_pp(r,sigma,T,N,theta,delta)
    om=Omega(r,sigma,T,N,theta,delta)
    A=np.dot(np.linalg.inv(om-A_pprime),A_p)
    b=np.transpose([1/T for i in range(T)])
    b=np.dot(np.linalg.inv(om-A_pprime),b)
    return A,b
  return system(r,sigma,T,N,theta,delta)   

class product():
  def __init__(self,nom):
    self.name_file=nom
    self.import_payoff()
  def import_payoff(self):
    fichier=pd.read_csv(self.name_file)
    setattr(self,"r",len(fichier["Interest rate"]))    
    setattr(self,"payoff",fichier["Payoff"])
    setattr(self,"sigma",fichier["Sigma"])
    setattr(self,"delta",fichier["Delta"])#check if we can put only a number i.o an array
    setattr(self,"theta",fichier["Theta"]) #check if we can put only a number i.o an array
    setattr(self,"N",len(fichier["Payoff"]))

fichier = "Data.xlsx"
