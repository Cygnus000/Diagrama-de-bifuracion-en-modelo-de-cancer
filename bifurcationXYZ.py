# codigo adaptado de https://stackoverflow.com/questions/62842261/in-python-how-to-make-a-bifurcation-diagram-of-the-lorenz-system-under-a-varyin
import numpy as np
import matplotlib.pyplot as plt
    # Definicion del sistema de cancer
    # aqui el parametro que vamos a variar es 'a' y los demas fijos
def cancer_system(x, y, z, a, b=2.5, c=0.6, d=1.5, e=4.5, f=1, g=0.2, h=0.5):
    x_dot = x * (1-x)- a * x * y - b * x * z
    y_dot = c * y * (1 - y) - d * x * y
    z_dot = (e * x * z) / (x + f) - g * x * z - h * z
    return x_dot, y_dot, z_dot

da = 0.0001#0.001  # tamaño de paso
a = np.arange(0.85, 1.05, da)  # rango de variacion  0.85 max 1.05
dt = 0.001  # paso de tiempo
t = np.arange(0, 250, dt)# 100, dt)  # rango temporal

# arreglos de solucion
xs = np.empty(len(t) + 1)
ys = np.empty(len(t) + 1)
zs = np.empty(len(t) + 1)

# valores iniciales x0,y0,z0 para el sistema
xs[0], ys[0], zs[0] = (1,1,1)

# Guardamos las coordenadas de los puntos y graficamos con plt.plot
ax_maxes = []
ax_mins = []
x_maxes = []
x_mins = []

ay_maxes = []
ay_mins = []
y_maxes = []
y_mins = []

az_maxes = []
az_mins = []
z_maxes = []
z_mins = []

for A in a:
    print(f"{A=:.4f}") # Imprimimos el avance en el parametro a
    for i in range(len(t)):
        # Runge-kutta para resolver el sistema a lo largo de un tiempo t
        x_dot, y_dot, z_dot = cancer_system(xs[i], ys[i], zs[i], A)
        xs[i + 1] = xs[i] + (x_dot * dt)
        xs[i + 1] = 0.5 *  (xs[i] + (x_dot * dt) + xs[i + 1])
        ys[i + 1] = ys[i] + (y_dot * dt)
        ys[i + 1] = 0.5 *  (ys[i] + (y_dot * dt) + ys[i + 1])
        zs[i + 1] = zs[i] + (z_dot * dt)
        zs[i + 1] = 0.5 *  (zs[i] + (z_dot * dt) + zs[i + 1])
    # calculamos y guardamos los puntos de inflexion
    for i in range(1, len(xs) - 1):
        # guardamos los maximos locales
        if xs[i - 1] < xs[i] and xs[i] > xs[i + 1]:
            ax_maxes.append(A)
            x_maxes.append(xs[i])
        # guardamos los minimos locales
        elif xs[i - 1] > xs[i] and xs[i] < xs[i + 1]:
            ax_mins.append(A)
            x_mins.append(xs[i])
 
    for i in range(1, len(ys) - 1):
        # guardamos los maximos locales
        if ys[i - 1] < ys[i] and ys[i] > ys[i + 1]:
            ay_maxes.append(A)
            y_maxes.append(ys[i])
        # guardamos los minimos locales
        elif ys[i - 1] > ys[i] and ys[i] < ys[i + 1]:
            ay_mins.append(A)
            y_mins.append(ys[i])
            
    for i in range(1, len(zs) - 1):
        # guardamos los maximos locales
        if zs[i - 1] < zs[i] and zs[i] > zs[i + 1]:
            az_maxes.append(A)
            z_maxes.append(zs[i])
        # guardamos los minimos locales
        elif zs[i - 1] > zs[i] and zs[i] < zs[i + 1]:
            az_mins.append(A)
            z_mins.append(zs[i])

    # "usamos los valores finales para usarlos como valores iniciales
    # para permanecer cerca del atractor durante la siguiente iteracion"
    xs[0], ys[0], zs[0] = xs[i], ys[i], zs[i]

# Eliminando los primeros puntos que contienen ruido
del ax_maxes[:10]
del ax_mins[:10]
del x_maxes[:10]
del x_mins[:10]
del ay_maxes[:10]
del ay_mins[:10]
del y_maxes[:10]
del y_mins[:10]
del az_maxes[:10]
del az_mins[:10]
del z_maxes[:10]
del z_mins[:10]
# salvando los datos en archivo
xpuntos = []
ypuntos = []
zpuntos = []
for valor in range(len(ax_maxes)):
    xpuntos.append((ax_maxes[valor],x_maxes[valor]))
for valor in range(len(ax_mins)):
    xpuntos.append((ax_mins[valor],x_mins[valor]))
np.savetxt('Bif_X_vs_a.dat', xpuntos)

for valor in range(len(ay_maxes)):
    ypuntos.append((ay_maxes[valor],y_maxes[valor]))
for valor in range(len(ay_mins)):
    ypuntos.append((ay_mins[valor],y_mins[valor]))
np.savetxt('Bif_Y_vs_a.dat', ypuntos)

for valor in range(len(az_maxes)):
    zpuntos.append((az_maxes[valor],z_maxes[valor]))
for valor in range(len(az_mins)):
    zpuntos.append((az_mins[valor],z_mins[valor]))
np.savetxt('Bif_Z_vs_a.dat', zpuntos)

# Graficamos
    # tipo de letra latex
plt.rcParams.update({
"text.usetex" : True,
"font.family" : 'serif',
"font.serif"  : 'cm'
})
# X
plt.scatter(ax_maxes, x_maxes, color="black", s=0.8, alpha=0.2)
plt.scatter(ax_mins, x_mins, color="red", s=0.8, alpha=0.2)
plt.legend(['$max(x(t))$','$min(x(t))$'], loc="best") 
#plt.xlim(0.85, 1.1)
#plt.ylim(0, 1)
plt.title("Diagrama de bifurcación del modelo de Cáncer")
plt.xlabel("parámetro $a$")
plt.ylabel("$x(t)$")
plt.savefig('bifurcation_X_vs_a.eps', format='eps')
plt.show()
plt.clf()
# Y
plt.scatter(ay_maxes, y_maxes, color="black", s=0.8, alpha=0.2)
plt.scatter(ay_mins, y_mins, color="red", s=0.8, alpha=0.2)
plt.legend(['$max(y(t))$','$min(y(t))$'], loc="best") 
#plt.xlim(0.85, 1.1)
#plt.ylim(0, 1)
plt.title("Diagrama de bifurcación del modelo de Cáncer")
plt.xlabel("parámetro $a$")
plt.ylabel("$y(t)$")
plt.savefig('bifurcation_Y_vs_a.eps', format='eps')
plt.show()
plt.clf()
# Z
plt.scatter(az_maxes, z_maxes, color="black", s=0.8, alpha=0.2)
plt.scatter(az_mins, z_mins, color="red", s=0.8, alpha=0.2)
plt.legend(['$max(z(t))$','$min(z(t))$'], loc="best") 
#plt.xlim(0.85, 1.1)
#plt.ylim(0, 1)
plt.title("Diagrama de bifurcación del modelo de Cáncer")
plt.xlabel("parámetro $a$")
plt.ylabel("$z(t)$")
plt.savefig('bifurcation_Z_vs_a.eps', format='eps')
plt.show()
