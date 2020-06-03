import math
import matplotlib.pyplot as plt
import numpy as np
import sympy

coefficients_ch_eq = [1, 125, 5100, 65000, 1]
poles = np.roots(coefficients_ch_eq)


def get_eq(coefficients):
    eq = ''
    size = len(coefficients)
    for i in range(size):
        if i != 0:
            if coefficients[i] > 0:
                eq += '+'
            else:
                eq += '-'
        if i < size-1:
            if coefficients[i] != 1:
                eq += str(coefficients[i]) + '*s**'
            else:
                eq += 's**'
            eq += str(size-i-1)
        else:
            eq += str(coefficients[i])
    return eq


def derivative_eq(coefficients):
    size = len(coefficients)-1
    der = []
    for i in range(size):
        der.append(-coefficients[i]*(size-i))
    return der


def break_away():
    break_point = np.roots(derivative_eq(coefficients_ch_eq))
    eq = get_eq(coefficients_ch_eq)
    expresion = []
    for s in break_point:
        temp = eval(eq)
        temp -= 1
        temp *= -1
        if temp > 0:
            expresion.append(s)
    return expresion


def routh():
    x = []
    temp = []
    for i in range(0,len(coefficients_ch_eq),2):
        if i == len(coefficients_ch_eq) - 1:
            temp.append('k')
        else:
            temp.append(str(coefficients_ch_eq[i]))
    x.append(temp.copy())
    temp.clear()
    for i in range(1,len(coefficients_ch_eq),2):
        if i == len(coefficients_ch_eq) -1:
            temp.append('k')
        else:
            temp.append(str(coefficients_ch_eq[i]))
    x.append(temp.copy())
    while len(x[0]) != len(x[1]):
        x[1].append('0')
    temp.clear()
    for i in range(len(coefficients_ch_eq) - 2):
        for j in range(len(x[i])-1):
            strr = '(' + x[1+i][0] + '*' + x[i][j+1] + '-' + x[i][0] + '*' + x[1+i][j+1]+')/' + x[1+i][0]
            try:
                temp.append(str(eval(strr)))
            except NameError:
                temp.append(strr)
        temp.append('0')
        x.append(temp.copy())
        temp.clear()
    k = sympy.solve(x[len(x)-2][0])
    return math.sqrt(k[0]/eval(x[len(x)-3][0]))


def get_centroid():
    centroid = 0
    for i in poles:
        print(i)
        centroid += i.real
    return centroid/len(poles)


def angles_asymptotes():
    angles = []
    for i in range(len(coefficients_ch_eq) - 1):
        angles.append((180*(2*i+1)/(len(coefficients_ch_eq) - 1))%360)
        angles.append(-angles[len(angles)-1])
    res = []
    for i in angles:
        if res.count(i) == 0 and res.count(360-abs(i)) == 0:
            res.append(i)
    return res


def angles_departure():
    angles = []
    counter = 0
    for i in poles:
        if i.imag != 0:
            angle = 0
            for j in poles:
                if i == j:
                    continue
                angle += math.pi/2 + math.atan2(abs(i.real-j.real),abs(i.imag-j.imag))
            angles.append([(math.pi - angle)*180/math.pi,counter])
        counter += 1
    return angles


def draw_poles():
    x = []
    y = []
    plt.grid()
    for pole in poles:
        x.append(pole.real)
        y.append(pole.imag)
    plt.title('root locus')
    plt.plot([-90,40],[0,0], color='black', linewidth=.25)
    plt.plot([0,0], [-40, 40], color='black', linewidth=.25)
    plt.scatter(x, y, color='red', marker='x')
    plt.xlim(-100,70)
    plt.ylim(-100,100)


def draw_asymptotes():
    draw_poles()
    centroid = get_centroid()
    angles = angles_asymptotes()
    for i in angles:
        plt.plot([centroid, centroid + 40], [0, math.tan(i*math.pi/180)*-40],linestyle='dotted', linewidth=1, color='black')
        plt.plot([centroid, centroid - 40], [0, math.tan(i * math.pi / 180) * 40], linestyle='dotted', linewidth=1,color='black')


def draw_root_locus():
    draw_asymptotes()
    x = []
    y = []
    for i in poles:
        if i.imag == 0:
            x.append(i.real)
            y.append(0)
    plt.plot(x,y,color='blue',linewidth=2)
    a = np.array([[1.0,-routh(),routh()**2],[1.0,routh(),routh()**2], [1.0,0,0]])
    b = np.array([0,0,break_away()[0]])
    z = np.linalg.solve(a,b)
    expre = str(z[0])+'+'+str(z[2])+'*y**2'
    y = np.linspace(-routh()-25,routh()+25,100)
    x = eval(expre)
    plt.plot(x,y,color='blue',linewidth=2)
    xax = []
    yax = []
    for pole in poles:
        if pole.imag != 0:
            xax.append(pole.real)
            yax.append(pole.imag)
    xax.append(xax[0]+math.tan((angles_departure()[0][0]+270)*math.pi/180)*yax[0]-4.75)
    yax.append(0)
    a = np.array([[1,0,0], [1,yax[0],yax[0]**2], [1,yax[1],yax[1]**2]])
    b = np.array([xax[2],xax[0],xax[1]])
    z = np.linalg.solve(a, b)
    expre = str(z[0]) + '+'+ str(z[1]) + '*y+' + str(z[2]) + '*y**2'
    y = np.linspace(yax[0], yax[0] + 25, 50)
    x = eval(expre)
    plt.plot(x, y, color='blue', linewidth=2)
    y = np.linspace(yax[1] - 25, yax[1], 50)
    x = eval(expre)
    plt.plot(x, y, color='blue', linewidth=2)


def k_100value():
    k = np.linspace(0,10000000,10000)
    for i in k:
        coefficients_ch_eq[len(coefficients_ch_eq)-1] = i
        roots = np.roots(coefficients_ch_eq)
        plt.plot(roots[0].real, roots[0].imag, color='green', linewidth=1, marker='*')
        plt.plot(roots[1].real, roots[1].imag, color='yellow', linewidth=1, marker='*')
        plt.plot(roots[2].real, roots[2].imag, color='gray', linewidth=1, marker='*')
        plt.plot(roots[3].real, roots[3].imag, color='brown', linewidth=1, marker='*')


draw_root_locus()
k_100value()
plt.show()