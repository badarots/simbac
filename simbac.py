import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.image as im

class Sim:    
    # inicializa a simulação usando uma matriz com as posições iniciais das células
    def __init__(self, m, **kwargs):
        
        self.cmax = kwargs.get('cmax', 3.5)
        self.cdif = kwargs.get('cdif', 0.25)
        self.a = kwargs.get('a', 0.181)
        self.q = kwargs.get('q', 4.35)
        self.k = kwargs.get('k', 1/self.a)
        self.r = kwargs.get('r', 1)
            
        self.m = np.array(m)
        self.c = kwargs.get('cinit', self.cmax)
        self.c = np.zeros(self.m.shape) + self.c
        
        self.xm, self.ym = m.shape
    
        self.pop = np.count_nonzero(self.m)
        np.random.seed()
    
    # Roda um passo da simulação
    def att(self):
        pop_x, pop_y = np.where(self.m != 0)   #encontra onde estão as celulas ocupadas
        self.pop = len(pop_x)         
        itera = np.random.permutation(np.arange(self.pop))   # gera um vetor com ordenação aleatória
        
        for i in itera: # percorre as posições ocupadas aleatóriamente
            self.vida_morte(pop_x[i], pop_y[i])
        
        self.c += self.cdif*(self.cmax - self.c)
        
    # decide o destino da célula na posição (i, j)
    def vida_morte(self, i, j):
        xi, xf, yi, yf = self.fronteira(i, j) # encontra a fronteira da vizinhança, serve para o caso da célula estar na borda da matriz
        
        vizinhos = self.m[i+xi:i+xf, j+yi:j+yf]  # estrai da matriz princial os vizinhos da célula
        
        zeros = np.where(vizinhos == 0)  # encontra os espaços vazios na vizinhança
        vagas = len(zeros[0])            # encontra o número de espaços vazinhos
        
        comida = np.sum(self.c[i+xi:i+xf, j+yi:j+yf])
                
        p = self.p(comida)
        destino = np.random.random()
        
        if self.r*p <= destino < self.r*p + self.a*(1-p):
            self.m[i, j] = 0
            self.pop -= 1
        else:
            fome = self.a
            reproduziu = False 
            if destino < self.r*p:
                fome = 1 + self.a
                if vagas != 0:
                    r = np.random.randint(len(zeros[0]))
                    vizinhos[zeros[0][r], zeros[1][r]] = self.m[i, j]
                    reproduziu = True
                else:
                    Δi = Δj = 0
                    while Δi == 0 and Δj == 0:
                        Δi = np.random.randint(xi, xf)
                        Δj = np.random.randint(yi, yf)
                        
                    x = i + Δi
                    y = j + Δj
                    while 0 <= x < self.xm and 0 <= y < self.ym:                    
                        if self.m[x, y] == 0:
                            self.m[x, y] = self.m[i, j]
                            reproduziu = True
                            break
                        x += Δi
                        y += Δj
                
            if reproduziu:
                self.pop += 1
            
            fome *= self.q
            if comida < fome:
                self.c[i+xi:i+xf, j+yi:j+yf] = 0
            else:                
                x = np.arange(xi, xf) + i
                y = np.arange(yi, yf) + j
                size = (xf - xi)*(yf - yi)
                pos = np.random.permutation(np.array(np.meshgrid(x, y)).T.reshape((size,2)))

                for x in pos:
                    if self.c[x[0], x[1]] >= fome:
                        self.c[x[0], x[1]] -= fome
                        break
                    else:
                        fome -= self.c[x[0], x[1]]
                        self.c[x[0], x[1]] = 0                
            
                
    # encontra a fronteira da vizinhaça, serve para o caso da célula estar na borda da matriz    
    def fronteira(self, i, j):
        xi = yi = -1
        xf = yf = 2

        if i == 0:
            xi = 0
        elif i + 1 == self.xm:
            xf = 1

        if j == 0:
            yi = 0
        elif j + 1 == self.ym:
            yf = 1

        return (xi, xf, yi, yf)
    
    def p(self, comida):
        #comida /= self.fmin
        return comida/(comida + self.k)
# Configuração das cores do gráfico de matrizes
cmap = colors.ListedColormap(['white', 'green','red','blue'])
bounds=[0,1,2,3]
norm = colors.BoundaryNorm(bounds, cmap.N)

# Função para plotar gráficos
def grafico(pop, sim, eixo = '', passos = -1):
    cmax = np.max(sim.cmax)
    # Plota os gráficos
    plt.close('all')
    fig = plt.figure(figsize = (7,7))     # cria a imagem

    ax = fig.add_subplot(221)       # gera o gráfico da matriz
    a1 = plt.imshow(sim.m, interpolation='nearest', origin='lower',
                    cmap=cmap, norm=norm)               
    plt.title('Pop = ' + str(sim.pop) + ', t = ' + str(len(pop) -1))

    fig.add_subplot(222)
    a2 = plt.imshow(sim.c, origin='lower', vmax = cmax, vmin = 0)
    plt.title('Comida')

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")

    cbar = plt.colorbar(cax=cax)
    ticks = int(np.max(sim.cmax)+0.5)
    cbar.set_ticks(np.linspace(0, cmax, 5))

    fig.add_subplot(212)      # gera o gráfico da população
    if eixo == 'log':
        a3, = plt.plot(np.log(pop))
        plt.ylabel('log(N)')
        plt.ylim(0,np.log(sim.m.shape[0]**2))
    else:
        a3, = plt.plot(pop)
        plt.ylabel('N')
        plt.ylim(0, sim.m.shape[0]**2)
    plt.xlabel('t')
    if passos != 0:
        plt.xlim(0,passos)

    plt.tight_layout()
    return (fig, ax, a1, a2, a3)
    
# Configuração das cores do gráfico de matrizes
cmap = colors.ListedColormap(['white', 'green','red','blue'])
bounds=[0,1,2,3]
norm = colors.BoundaryNorm(bounds, cmap.N)

# Função para plotar gráficos
def grafico(pop, sim, eixo = '', passos = -1):
    cmax = np.max(sim.cmax)
    # Plota os gráficos
    plt.close('all')
    fig = plt.figure(figsize = (7,7))     # cria a imagem

    ax = fig.add_subplot(221)       # gera o gráfico da matriz
    a1 = plt.imshow(sim.m, interpolation='nearest', origin='lower',
                    cmap=cmap, norm=norm)               
    plt.title('Pop = ' + str(sim.pop) + ', t = ' + str(len(pop) -1))

    fig.add_subplot(222)
    a2 = plt.imshow(sim.c, origin='lower', vmax = cmax, vmin = 0)
    plt.title('Comida')

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")

    cbar = plt.colorbar(cax=cax)
    ticks = int(np.max(sim.cmax)+0.5)
    cbar.set_ticks(np.linspace(0, cmax, 5))

    fig.add_subplot(212)      # gera o gráfico da população
    if eixo == 'log':
        a3, = plt.plot(np.log(pop))
        plt.ylabel('log(N)')
        plt.ylim(0,np.log(sim.m.shape[0]**2))
    else:
        a3, = plt.plot(pop)
        plt.ylabel('N')
        plt.ylim(0, sim.m.shape[0]**2)
    plt.xlabel('t')
    if passos != 0:
        plt.xlim(0,passos)

    plt.tight_layout()
    return (fig, ax, a1, a2, a3)
    
# Matrizes desenhadas
    
def circulo(n):
    a=np.zeros((n,n))
    m = n//2
    r = n*3/10

    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            if np.round(np.sqrt((i-m)**2 + (j-m)**2)) <= r:
                a[i, j] = 1
    return a

def cor1(n):
    a = np.zeros((n,n))
    mx = n//2
    my = round(n*0.75)

    b = 0.18
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            y = i - my
            x = j - mx
            if x != 0 or y != 0:
                sen = y/np.sqrt(x**2 + y**2)
                cos = x/np.sqrt(x**2 + y**2)
                r = np.sqrt(x**2 + y**2)
                if 2 - 2*sen + sen*np.sqrt(abs(cos))/(sen + 1.4) -r/(n*b) > 0:
                    a[i, j] = 1
    return a

def cor2(n):
    a=np.zeros((n,n))
    mx = n//2
    my = round(n*0.4)

    b = 2.8
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            y = (i - my)/(n/b)
            x = (j - mx)/(n/b)
            if (x**2 + y**2 -1)**3 - 2 * x**2 * y**3 < 0:
                a[i, j] = 1
    return a

def cos(n):
    a=np.zeros((n,n))
    mx = n//2
    my = n//2
    r0 = np.sqrt(2)*n/2
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            y = (i - my)
            x = (j - mx)
            r = np.sqrt(x**2 + y**2)
            a[i, j] = (np.cos(np.pi*r/r0)+1)/2
    return a
def cos2(n):
    a=np.zeros((n,n))
    mx = n//2
    my = n//2
    r0 = n/2
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            y = (i - my)
            x = (j - mx)
            r = np.sqrt(x**2 + y**2)
            if r <= r0:
                a[i, j] = (np.cos(np.pi*r/r0)+1)/2
    return a
def cos3(n):
    a=np.zeros((n,n))
    mx = n//2
    my = n//2
    r0 = np.sqrt(2)*n/2
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            y = (i - my)
            x = (j - mx)
            r = np.sqrt(x**2 + y**2)
            a[i, j] = np.cos(np.pi/2*r/r0)
    return a
def cos4(n):
    a=np.zeros((n,n))
    mx = n//2
    my = n//2
    r0 = 0.9*n/2
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            y = (i - my)
            x = (j - mx)
            r = np.sqrt(x**2 + y**2)
            if r <= r0:
                a[i, j] = np.cos(np.pi/2*r/r0)
    return a

def tri(n):
    a=np.zeros((n,n))
    mx = 0
    my = n//2
    r0 = n/2
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            y = (i - my)
            x = (j - mx)
            r = abs(r0 - x)
            if -r <= y <= r:
                a[i, j] = 1
    return a

def sig(x, ri, rf, a, b):
    return ri + (rf-ri)/(1+exp(b(a-x)))

# matriz da imagem 'caminho<n>.png'
def caminho(n):
    caminho = im.imread('caminho{}.png'.format(n))
    caminho =  1 -np.sum(caminho, axis = 2)/3
    return caminho