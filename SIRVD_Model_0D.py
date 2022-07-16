import numpy as np
import matplotlib.pyplot as plt

t_tot = 2000
dt = 0.1
t = np.arange(0, t_tot+0.01, dt)
n = len(t)

S = np.empty(n)
I = np.empty(n)
R = np.empty(n)
V = np.empty(n)
D = np.empty(n)

'''rate of infection'''
beta = np.empty(n)
beta[:] = 0.4417

'''recovery rate'''
gamma = np.empty(n)
gamma[:] = 0.1508

'''death rate'''
delta = np.empty(n)
delta[:] = 0.01

''' vaccination rate'''
alpha = np.empty(n)
alpha[:] = 0.0

'''susceptibility rate - recovered to susceptible'''
sigma = np.empty(n)
sigma[:] = 0.003

'''total Population'''
N = 500
S[0] = N

'''Phasenraumplot Visualisierung'''
'''def show_phasenraumplot(Cx, Cy, Cx0, Cy0, A, B, t):
    #plt.clf()
    plt.title("Phasenraumplot mit A=" + str(A) + ", B=" + str(B) + ", [X0]="
              + str(Cx0) + ", [Y0]=" + str(Cy0))
    plt.plot(Cx, Cy,linewidth=0.6,label=["A=" + str(A) + " B=" + str(B)])
    #plt.grid(linestyle='dotted', linewidth=0.5)
    plt.xlabel('[X]')
    plt.ylabel('[Y]')
    plt.tight_layout()
    #plt.draw()
    #plt.savefig("Phasenraumplot_A=" + str(A) + "_B=" + str(B) + "_[X0]=" + str(Cx0) + "_[Y0]=" + str(Cy0) + ".pdf", dpi=1000)
    plt.show()
    return'''

'''Konzentrationsverlauf Visualisierung'''
def show_konzentrationsverlauf(S, R, V, I, D, t):
    plt.clf()
    plt.figure(figsize=(18, 8))

    #plt.title("Konzentrationsverlauf mit A=" + str(A) + ", B=" + str(B) +", [X0]=" + str(Cx0) + ", [Y0]=" + str(Cy0))
    plt.plot(t, S, label='susceptible',linewidth=1.5)
    plt.plot(t, I, label='infected',linewidth=1.5)
    plt.plot(t, R, label='recovered',linewidth=1.5)
    plt.plot(t, V, label='vaccinated',linewidth=1.5)
    plt.plot(t, D, label='deceased',linewidth=1.5)
    #plt.margins(x=-0.2,y=0.05)
    plt.legend()
    plt.grid(linestyle='dotted', linewidth=0.25)
    plt.ylabel('n')
    plt.xlabel('t')

    plt.tight_layout()
    #plt.draw()
    #plt.savefig("render/gesamtplot.png", dpi=1000)
    plt.show()
    return

'''Iterationsschleife'''
for i in range(1, n):
    S[i] = S[i-1] + (- ((beta[i-1] * I[i-1] * S[i-1]) / (N)) + sigma[i-1] * R[i-1] - alpha[i-1] * S[i-1]) * dt
    I[i] = I[i-1] + (((beta[i-1] * I[i-1] * S[i-1]) / (N)) - gamma[i-1] * I[i-1] - delta[i-1] * I[i-1]) * dt
    R[i] = R[i-1] + (gamma[i-1] * I[i-1] - sigma[i-1] * R[i-1]) * dt
    V[i] = V[i-1] + (alpha[i-1] * S[i-1]) * dt
    D[i] = D[i-1] + (delta[i-1] * I[i-1]) * dt

    if i == 150:
        S[i]-=1
        I[i]+=1


show_konzentrationsverlauf(S, R, V, I, D, t)






