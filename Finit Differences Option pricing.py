#%% Cell 1
import numpy as np
import pandas as pd
import os
os.chdir(r"C:\Users\nadja\OneDrive\Bureau\code")
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
        a=(1-theta)*dt*(alpha/dx**2-beta/dx+gamma)-1
        b=theta*dt*(beta/dx-2*alpha/dx**2)
        c=(1-theta)*dt*(beta/dx-2*alpha/dx**2)
        d=theta*dt*alpha/dx**2
        e=(1-theta)*dt*alpha/dx**2
    else :
        omega=-theta*dt*(alpha/dx**2+beta/dx+gamma)-1
        a=(1-theta)*dt*(alpha/dx**2+beta/dx+gamma)-1
        b=-theta*dt*(2*alpha/dx**2+beta/dx)
        c=-(1-theta)*dt*(2*alpha/dx**2+beta/dx)
        d=theta*dt*alpha/dx**2
        e=(1-theta)*dt*alpha/dx**2
    return(omega,a,b,c,d,e)
         
#%%
def matrices(r,sigma,T,N,theta,delta,X_N):
  def A_prime(r,sigma,T,N,theta,delta):
      A_p = np.zeros((N+1,N+1))
      diag_low_first = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[5] for i in range(1,N)]
      diag_low=np.append(diag_low_first,[compute_coeffs(r, sigma, T, N, 0, theta, delta)[3]])
      
      diag_first = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[1] for i in range(N)]
      diag_first[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[5]
      diag=np.append(diag_first,[compute_coeffs(r, sigma, T, N, 0, theta, delta)[5]])
      
      
      diag_up = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[4] for i in range(0,N)]
      diag_up[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[3]
      
      np.fill_diagonal(A_p[1:,:-1],diag_low)
      np.fill_diagonal(A_p,diag)
      np.fill_diagonal(A_p[:-1,1:],diag_up)

      return A_p
  def A_pp(r,sigma,T,N,theta,delta):
    A_2p = np.zeros((N+1,N+1))
    diag_low = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[2] for i in range(1,N+1)]
    
    diag_first=np.zeros(N-1)
    diag=np.append(diag_first,[compute_coeffs(r, sigma, T, N, 0, theta, delta)[4]])
    #print(compute_coeffs(r, sigma, T, N, N, theta, delta)[4])
    diag[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[4]

    diag_up = [compute_coeffs(r, sigma, T, N, N-i, theta, delta)[3] for i in range(0,N)]
    diag_up[0]=compute_coeffs(r, sigma, T, N, N, theta, delta)[2]

    np.fill_diagonal(A_2p[1:,:-1],diag_low)
    np.fill_diagonal(A_2p,diag)
    np.fill_diagonal(A_2p[:-1,1:],diag_up)
    return A_2p
  def Omega(r,sigma,T,N,theta,delta):
    L_om= np.full(N+1,compute_coeffs(r, sigma, T, N, 1, theta, delta)[0]) #on met un indice random car coeff identique
    Om = np.diag(L_om)
    Om[0][0] = compute_coeffs(r, sigma, T, N, N, theta, delta)[0]
    Om[N][N] = compute_coeffs(r, sigma, T, N, 0, theta, delta)[0]
    return Om
  def system(r,sigma,T,N,theta,delta):
    A_p = A_prime(r,sigma,T,N,theta,delta)
    
    A_pprime = A_pp(r,sigma,T,N,theta,delta)
    
    om=Omega(r,sigma,T,N,theta,delta)
    A=np.dot(np.linalg.inv(om-A_pprime),A_p)
    b=np.transpose([1/T for i in range(T+1)])
    b=np.dot(np.linalg.inv(om-A_pprime),b)
    
    return A,b
  def solve_system(A,b,X,t,T):
      for i in range(T+1):
          X=np.dot(A,X)+b
          print(X)
      return X
  
  A=system(r,sigma,T,N,theta,delta)[0]
  b=system(r,sigma,T,N,theta,delta)[1]
  #print(solve_system(A,b,X_N,0,T)[0])
  return solve_system(A,b,X_N,0,T)   
    
class product():
  def __init__(self,nom):
    self.name_file=nom
    self.import_payoff()
  def import_payoff(self):
    fichier=pd.read_csv(self.name_file,sep=";")
    setattr(self,"r",fichier["Interest rate"][0])    
    setattr(self,"payoff",fichier["Payoff"])
    setattr(self,"sigma",fichier["Sigma"][0])
    setattr(self,"delta",fichier["Delta"][0])
    #check if we can put only a number i.o an array
    setattr(self,"theta",fichier["Theta"][0]) #check if we can put only a number i.o an array
    setattr(self,"N",np.shape(fichier["Payoff"])[0])
    setattr(self,"T",self.N)
    setattr(self,"result",matrices(self.r,self.sigma,self.T,self.N,self.theta,self.delta,self.payoff))

#%% Cellule Test
nom_fichier = "Data_pde_pricer.csv"
a=product(nom_fichier)
fichier = pd.read_csv(nom_fichier,sep=";")

#%% Cellule Test
nom_fichier = "Data_pde_pricer.csv"
fichier = pd.read_csv(nom_fichier,sep=";")
a=product(nom_fichier)
mat=a.result

#%% Cellule Test
b=np.full(11,1)
a=np.ones((11,11))
print(a)
print(b)
print(np.dot(a,b))
