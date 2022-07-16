import numpy as np
import matplotlib.pyplot as plt
import time
import random
#import cv2



'''Definition der Laufvariablen'''
t_tot = 3000
dt = 1
l_tot = 500
t = np.arange(0, t_tot+0.01, dt)
n = len(t)


'''rate of infection'''
beta = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
beta[:,:] = 0.4417

'''recovery rate'''
gamma = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
gamma[:,:] = 0.1508

'''death rate'''
delta = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
delta[:,:] = 0.01

''' vaccination rate'''
alpha = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
alpha[:,:] = 0.0

'''susceptibility rate - recovered to susceptible'''
sigma = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
sigma[:,:] = 0.003

'''total Population'''
N = 1000


'''Arrays fuer die beiden Konzentrationen: (Zeit, y-Laenge, x-Laenge)'''
S = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
I = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
R = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
V = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
D = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
Sn = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
In = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
Rn = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
Vn = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)
Dn = np.zeros(shape=(l_tot+1, l_tot+1), dtype=float)

'''Verschobene Konzentrationsarrays fuer die Diffusion'''
I_oben = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
I_unten = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
I_links = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)
I_rechts = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)

I_diff_gesamt = np.empty(shape=(l_tot+1, l_tot+1), dtype=float)

S[:,:] = N

DI = 0.15

betrachtungsintervall = 1
h = 0

plt.ion()


def update_figures(S, R, I , D, dt, i):
    figure1.set_data(S[:,:])
    figure2.set_data(R[:, :])
    figure3.set_data(I[:, :])
    figure4.set_data(D[:, :])
    fig.suptitle("t=" + str(np.round((dt * i), decimals=1)) + ", Iteration " + str(i), fontsize=14)
    #plt.tight_layout()
    plt.draw()
    #plt.savefig("render/image_gut_" + str(h) + ".png", dpi=250)
    plt.pause(0.001)
    return

# create figure
fig = plt.figure(figsize=(15, 14))

# setting values to rows and column variables
rows = 2
columns = 2

# Adds a subplot at the 1st position
fig.add_subplot(rows, columns, 1)

# showing image
figure1 = plt.imshow(S[:,:], cmap='plasma',interpolation="bilinear") #Greys bicubic bilinear
plt.axis('off')
plt.clim(0, N)
fig.colorbar(figure1)
plt.title("susceptible", fontsize=10)

# Adds a subplot at the 2nd position
fig.add_subplot(rows, columns, 2)

# showing image
figure2 = plt.imshow(R[:,:], cmap='plasma',interpolation="bilinear")
plt.axis('off')
plt.clim(0, N)
fig.colorbar(figure2)
plt.title("recovered", fontsize=10)

# Adds a subplot at the 3rd position
fig.add_subplot(rows, columns, 3)

# showing image
figure3 = plt.imshow(I[:,:], cmap='plasma',interpolation="bilinear")
plt.axis('off')
fig.colorbar(figure3)
plt.clim(0, int(N/7))
plt.title("infected", fontsize=10)

# Adds a subplot at the 4th position
fig.add_subplot(rows, columns, 4)

# showing image
figure4 = plt.imshow(D[:,:], cmap='plasma',interpolation="bilinear")
plt.axis('off')
plt.clim(0, int(N/6))
fig.colorbar(figure4)
plt.title("dead", fontsize=10)





fig.suptitle("t=0, Iteration 0", fontsize=12)

plt.tight_layout()

start = time.time()



print("Start der Berechnung...")
for i in range(1, n):
    I_oben[:, :] = np.roll(I[:, :], -1, axis=0)
    I_oben[-1, :] = I[-1, :]

    I_unten[:, :] = np.roll(I[:, :], 1, axis=0)
    I_unten[0, :] = I[0, :]

    I_links[:, :] = np.roll(I[:, :], -1, axis=1)
    I_links[:, -1] = I[:, -1]

    I_rechts[:, :] = np.roll(I[:, :], 1, axis=1)
    I_rechts[:, 0] = I[:, 0]

    I_diff_gesamt[:, :] = DI * (- 4 * I[:, :] + I_oben[:, :] + I_unten[:, :] + I_links[:, :] + I_rechts[:, :])

    Sn[:, :] = S[:, :] + (-((beta[:, :] * I[:, :] * S[:, :]) / N) + sigma[:, :] * R[:, :] - alpha[:, :] * S[:, :]) * dt
    In[:, :] = I[:, :] + (((beta[:, :] * I[:, :] * S[:, :]) / N) - gamma[:, :] * I[:, :] - delta[:, :] * I[:, :] + I_diff_gesamt[:, :]) * dt
    Rn[:, :] = R[:, :] + (gamma[:, :] * I[:, :] - sigma[:, :] * R[:, :]) * dt
    Vn[:, :] = V[:, :] + (alpha[:, :] * S[:, :]) * dt
    Dn[:, :] = D[:, :] + (delta[:, :] * I[:, :]) * dt

    S[:, :] = Sn[:, :]
    I[:, :] = In[:, :]
    R[:, :] = Rn[:, :]
    V[:, :] = Vn[:, :]
    D[:, :] = Dn[:, :]

    end = time.time()

    #if i == 50:
        #I[50,50] = 1
        #S[50, 50] = N-1
        #I[30, 80] = 1
        #S[30, 80] = 499
        #I[90, 10] = 1
        #S[90, 10] = 499

    if i % int(random.random() * 100 + 1) == 0 and i < 500:
        f1 = 1
        f2 = 1
        x1 = int(random.random() * l_tot)
        y1 = int(random.random() * l_tot)
        x2 = int(random.random() * l_tot)
        y2 = int(random.random() * l_tot)
        I[x1,y1] += f1
        S[x1,y1] -= f1
        I[x2,y2] += f2
        S[x2,y2] -= f2

    if i % ((t_tot/dt)/1000) == 0:
        end = time.time()
        print("Berechnung: " + str(np.round((i/n)*100, decimals=1)) + "% erledigt, Iteration " + str(i) + " mit " + str(np.round(i / (end - start), decimals=1)) + " Iters/s Restdauer: " + str(np.round((((t_tot/dt)-i)/(i / (end - start)))/60, decimals=1)) + " min")

    if i == 1:
        update_figures(S, R, I, D, dt, i)
        h += 1

    if i % betrachtungsintervall == 0:
        update_figures(S, R, I, D, dt, i)
        h += 1
