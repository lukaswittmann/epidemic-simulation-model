import numpy as np
import matplotlib.pyplot as plt



t_tot = 1000
dt = 0.5
t = np.arange(0, t_tot+0.01, dt)
n = len(t)

S = np.empty(n)
I = np.empty(n)
R = np.empty(n)
V = np.empty(n)
D = np.empty(n)


I_v = np.empty(n)
D_v = np.empty(n)

'''rate of infection'''
beta = np.empty(n)
beta[:] = 0.5

'''recovery rate'''
gamma = np.empty(n)
gamma[:] = 0.1

'''death rate'''
delta = np.empty(n)
delta[:] = 0.001

''' vaccination rate'''
alpha = np.empty(n)
alpha[:] = 0.00

'''susceptibility rate - recovered to susceptible'''
sigma = np.empty(n)
sigma[:] = 0.003




'''rate of infection for vaccinated'''
beta_v = np.empty(n)
beta_v[:] = beta[:] * 0.5

'''recovery rate for vaccinated'''
gamma_v = np.empty(n)
gamma_v[:] = gamma[:] * 10

'''death rate for vaccinated'''
delta_v = np.empty(n)
delta_v[:] = delta[:] * 0.1

'''susceptibility rate - vaccinated to susceptible'''
sigma_v = np.empty(n)
sigma_v[:] = sigma[:] * 1.2

'''total Population'''
N = 100
V[0] = 0.0*N
S[0] = N - V[0]


'''Konzentrationsverlauf Visualisierung'''
def show_konzentrationsverlauf(S, R, V, I, D, t, I_v, D_v):
    plt.figure(figsize=(24, 10))

    plt.title("Parameters used: " + "beta=" + str(beta[1]) + ", gamma=" +  str(gamma[1]) + ", delta=" +  str(delta[1]) + ", alpha=" +  str(alpha[1]) + ", sigma=" +  str(sigma[1]) + ", beta_v=" +  str(beta_v[1]) + ", gamma_v=" +  str(gamma_v[1]) + ", delta_v=" +  str(delta_v[1]) + ", sigma_v=" +  str(sigma_v[1]))

    plt.plot(t, S, label='susceptible',linewidth=1.5)
    plt.plot(t, I, label='infected unvaccinated',linewidth=1, linestyle='dotted')
    plt.plot(t, I_v, label='infected vaccinated', linewidth=1, linestyle='dotted')
    plt.plot(t, I_v + I, label='infected', linewidth=1.5, linestyle='dotted')
    plt.plot(t, R, label='recovered',linewidth=1.5)
    plt.plot(t, V, label='vaccinated',linewidth=1.5)
    plt.plot(t, D, label='deceased unvaccinated',linewidth=0.5, linestyle='dashed')
    plt.plot(t, D_v, label='deceased vaccinated', linewidth=0.5, linestyle='dashed')
    plt.plot(t, D_v + D, label='deceased', linewidth=1.5, linestyle='dashed')
    plt.grid(linestyle='dotted', linewidth=0.25)
    plt.ylabel('n [%]')
    plt.xlabel('t [a.u.]')
    plt.legend()
    #plt.tight_layout()
    plt.show()
    return

'''Iterationsschleife'''
for i in range(1, n):
    S[i] = S[i-1] + (- ((beta[i-1] * I[i-1] * S[i-1]) / (N)) + sigma[i-1] * R[i-1] - alpha[i-1] * S[i-1] + (sigma_v[i-1] * V[i-1])) * dt
    I[i] = I[i-1] + (((beta[i-1] * I[i-1] * S[i-1]) / (N)) - gamma[i-1] * I[i-1] - delta[i-1] * I[i-1]) * dt
    R[i] = R[i-1] + (gamma[i-1] * I[i-1] - sigma[i-1] * R[i-1]) * dt
    V[i] = V[i-1] + (alpha[i-1] * S[i-1] - (beta_v[i-1] * V[i-1]) + gamma_v[i-1] * I_v[i-1] - (sigma_v[i-1] * V[i-1])) * dt
    D[i] = D[i-1] + (delta[i-1] * I[i-1]) * dt

    I_v[i] = I_v[i - 1] + ((beta_v[i-1] * V[i-1]) - delta_v[i-1] * I_v[i-1] - gamma_v[i-1] * I_v[i-1]) * dt
    D_v[i] = D_v[i - 1] + (delta_v[i-1] * I_v[i-1]) * dt

    if i == 100:
        S[i] -= 0.01
        I[i] += 0.01

print("beta=" + str(beta[1]) + " gamma=" +  str(gamma[1]) + " delta=" +  str(delta[1]) + " alpha=" +  str(alpha[1]) + " sigma=" +  str(sigma[1]) + " beta_v=" +  str(beta_v[1]) + " gamma_v=" +  str(gamma_v[1]) + " delta_v=" +  str(delta_v[1]) + " sigma_v=" +  str(sigma_v[1]))

show_konzentrationsverlauf(S, R, V, I, D, t, I_v, D_v)






