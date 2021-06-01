import numpy as np
import matplotlib.pyplot as plt
#import pandas as pd
#import sympy as sy

k = 237 # coef de transf de calor # material
h = 70 # coeficiente de convectivo na extremidade direita
L = 0.2 # comprimento total da barra
nos = 5 # qnt. de nós
D = 0.019 # diâmetro
Te = 120 # temperatura na extremidade esquerda (°C)
Td = 25 # temperatura do ar (°C)
A = (((np.pi)*(D**2))/4)
P = np.pi*D

constantes = {
    'coefiência de condutividade térmica (k)' : k,
    'coeficiênte de convecção (h)' : h,
    'comprimento da barra': L,
    'quantidade de nós': nos,
    'diâmetro da colher': D,
    'temperatura da concha': Te,
    'temperatura da cabo': Td,
    'Área da seção transversal':A
    }

class temperatura_barra():
    def resultado(self, h, k, nos, L, A, P, Te, Td):
        e = nos - 1 # n de elementos
        self.barra = np.linspace(0,L,nos) # discretização da barra
        #self.barra = np.arange(0,L,0.025) # discretização da barra
        pts_d = np.zeros(len(self.barra)-2)
        pts_d = np.append(pts_d, 1)
        l = L/e #comprimento de cada elemento
        
        M1 = np.matrix([
            [1, -1],
            [-1, 1]
            ])
        
        M2 = np.matrix([
            [2, 1],
            [1, 2]
            ])
        
        vet1 = np.array([1,1])
        
        
        self.ks = []
        self.fs = []
        
        # montagem das matrizes locais
        
        for i in pts_d:
            var1 = k*A/l
            
            var2 = h*P*l/6
            
            var3 = h*A
            
            var4 = h*P*l*Td/2
            
            var5 = h*A*Td
            
            mat_convec = np.matrix([[0, 0],[0, var3]])
            
            vet_convec = np.array([0,var5])
            
            k_e = var1*M1 + var2*M2 + mat_convec*i
            
            f_e = var4*vet1 + vet_convec*i
            
            self.ks.append(k_e)
            self.fs.append(f_e)
            
        # montagem da matriz global K
        
        [n] = set(k_e.shape)
        
        N = len(self.ks) #quantidade de matrizes k_e
        
        self.K = np.zeros((N+1, N+1)) #matriz de zeros
        
        for i, k in enumerate(self.ks):
            self.K[i:i+n, i:i+n] += k
        
        # montagem da vetor global F
        
        [n] = set(tuple((f_e.shape[0],f_e.shape[0])))
        
        N = len(self.fs) # quantidade de matrizes f_e
        
        self.F = np.zeros((N + 1)) #array de zeros
        
        for i, k in enumerate(self.fs):
            self.F[i:i+n] += k

        #Q = np.array([Te, 0, 0, 0, 0])

        self.K = np.delete(self.K, 0, 0)
        self.K = np.delete(self.K, 0, 1)
        self.F = np.delete(self.F, 0, 0)
        #Q = np.delete(Q, 0, 0)
        
        self.F[0] = self.F[0] + Te*self.K[1,0]*(-1)

        R = np.linalg.inv(self.K) @ (self.F)
        
        R = np.round(R, 1)
        
        R = np.insert(R, 0, Te)
        
        return R
    
    
temp = temperatura_barra()

R = temp.resultado(h, k, nos, L, A, P, Te, Td)

fig = plt.figure()
ax = plt.axes([0.1,0.1,0.9,0.9])

def plot_h(start, end):
    for i in np.linspace(start,end,5):
        R = temp.resultado(i, k, nos, L, A, P, Te, Td)
        ax.plot(temp.barra, R, label= 'h = ' + str(i))

def plot_L(start, end):
    for i in np.round(np.linspace(start,end,5),2):
        R = temp.resultado(h, k, nos, i, A, P, Te, Td)
        ax.plot(temp.barra, R, label= 'L = ' + str(i))

def plot_k(start, end):
    for i in np.linspace(start,end,5):
        R = temp.resultado(h, i, nos, L, A, P, Te, Td)
        ax.plot(temp.barra, R, label= 'k = ' + str(i))

materiais = {
     'alumínio': 237,
     'aço inox': 15.1,
     'madeira': 0.16
     }


#plot_k(15.1,237)
#plot_h(10,70)
plot_L(0.1,0.2)
ax.set_xlabel('comprimento (m)')
ax.set_ylabel('temperatura (°C)')
ax.legend()
plt.show()

'''
data = {
        'barra': temp.barra,
        'temperaturas': R
        }

data = pd.DataFrame(data)

print(data)
'''

