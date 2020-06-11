import os
from random import random
from math import log
from math import exp
from math import pi
from numpy import linspace
from scipy import stats
import matplotlib.pyplot as plt
 
def uniformZeroOne():
    return random()
 
def uniform(size, a, b):
    xList = []
    for i in range(size):
        r = uniformZeroOne()
        x = a+(b-a)*r
        xList.append(x)
    return xList
 
def exponential(size, alpha):
    xList = []
    for i in range(size):
        r = uniformZeroOne()
        x = -(1/alpha)*log(r)
        xList.append(x)
    return xList
 
def gamma(size, k, alpha):
    xList = []
    for i in range(size):
        tr = 1
        for j in range(k):
            r = uniformZeroOne()
            tr = tr * r
        x = -(1/alpha)*log(tr)
        xList.append(x)
    return xList
 
# def pascal(size, k, p):
#   xList = []
#   qr = log(1 - p)
#   for i in range(size):
#       tr = 1
#       for j in range(k):
#           r = uniformZeroOne()
#           tr = tr * r
#       xList.append(int(log(tr)//qr))
#   return xList
 
def geometric(size, p):
    xList = []
    qr = log(1 - p)
    for i in range(size):
        r = uniformZeroOne()
        xList.append(int(log(r)//qr))
    return xList
 
def pascal(size, k, p):
    xList = []
    for i in range(size):
        x = sum(geometric(k,p))
        xList.append(x)
    return xList
 
def binom(size, n, p):
    xList = []
    for i in range(size):
        x = 0
        for j in range(n):
            r = uniformZeroOne()
            if (r <= p):
                x += 1
        xList.append(x)
    return xList
 
def poisson(size, lambda_):
    xList = []
    b = exp(-lambda_)
    for i in range(size):
        x = 0
        tr = 1
        while (True):
            r = uniformZeroOne()
            tr = tr * r
            if (tr >= b):
                x += 1
            else:
                break
        xList.append(x)
    return xList
 
def getAcumulada(array):
    acumulada = [array[0]]
    for i in range(1,len(array)):
        acumulada.append(acumulada[i-1]+array[i])
    return acumulada
 
def empirical(size, probabilities):
    acumulada = getAcumulada(probabilities)
    xList = []
    for i in range(size):
        r = uniformZeroOne()
        x = 1
        while (x <= len(acumulada)):
            if(r > acumulada[x-1]):
                x += 1
            else:
                break
        xList.append(x)
    return xList
 
def normal(size, mu, sigma):
    k = 12
    xList = []
    for i in range(size):
        sum = 0
        for j in range(k):
            sum += uniformZeroOne()
        x = sigma*(sum-6) + mu
        xList.append(x)
    return xList
 
def hypergeometric(size, N, n, Np):
    xList = []
    p = Np/N
    for i in range(size):
        x = 0
        Ni = N
        Pi = p
        for j in range(n):
            r = uniformZeroOne()
            if (r <= Pi):
                s = 1
                x += 1
            else:
                s = 0
            Pi = (Ni*Pi - s)/(Ni-1)
            Ni = Ni - 1
        xList.append(x)
    return xList
 
def testChiCuad(f_obs, f_esp, alpha):
    chi2 = 0
    k = len(f_obs)
    for i in (range(k)):
        chi2 += ((f_obs[i] - f_esp[i])**2)/f_esp[i]
    chiTabla = stats.chi2.ppf(1-alpha,k-1)
    print ("=== Chi Calculado ===")
    print(chi2)
    print("")
    print("=== Chi Tabla ===")
    print(chiTabla)
     
    if chi2<chiTabla:
        print ("Al ser el chi2 calculado menor al valor de tabla, la hipótesis nula de que no existe diferencia entre la distribución de la muestra y la distribución supuesta se acepta")  
    else:
        print ("Al ser el chi2 calculado mayor al valor de tabla, la hipótesis nula de que no existe diferencia entre la distribución de la muestra y la distribución supuesta se rechaza")
 
def kstest(lista, alpha, parameters, type):
    print("Realizando test de Kolmogorov-Smirnov, esto puede tardar unos segundos si la muestra es grande")
    lista.sort()
    dmax = 0
    n = len(lista)
    if (type == 'Exponential'):
        for j in range(n):
            d = (j+1)/n - (1-exp(-parameters*lista[j]))
            if (d > dmax):
                dmax = d
            d = 1-exp(-parameters*lista[j]) - j/n
            if (d > dmax):
                dmax = d
    elif (type == 'Normal'):
        for j in range(n):
            d = (j+1)/n - stats.norm.cdf(lista[j],parameters[0],parameters[1])
            if (d > dmax):
                dmax = d
            d = stats.norm.cdf(lista[j],parameters[0],parameters[1]) - j/n
            if (d > dmax):
                dmax = d
    elif (type == 'Gamma'):
        loc = 0
        scale = 1/parameters[0]
        for j in range(n):
            d = (j+1)/n - stats.gamma.cdf(lista[j],parameters[1],loc,scale)
            if (d > dmax):
                dmax = d
            d = stats.gamma.cdf(lista[j],parameters[1],loc,scale) - j/n
            if (d > dmax):
                dmax = d
    elif (type == "Uniform"):
        for j in range(n):
            d = (j+1)/n - (lista[j]-parameters[0])/(parameters[1]-parameters[0])
            if (d > dmax):
                dmax = d
            d = (lista[j]-parameters[0])/(parameters[1]-parameters[0]) - j/n
            if (d > dmax):
                dmax = d
    dtabla = stats.ksone.ppf((1-alpha/2),n)
    print ("=== Desviación máxima ===")
    print(dmax)
    print("")
    print ("=== Desviación tabla ===")
    print(dtabla)
    if dmax<dtabla:
        print ("Al ser el desvío calculado menor al valor de tabla, la hipótesis nula de que no existe diferencia entre la distribución de la muestra y la distribución uniforme se acepta")  
    else:
        print ("Al ser el desvío calculado mayor al valor de tabla, la hipótesis nula de que no existe diferencia entre la distribución de la muestra y la distribución uniforme se rechaza")
 
def calculateIntervals(array,numOfInterv):
    minimo = array[0]
    maximo = array[0]
    for i in range(1, len(array)):
        if array[i] < minimo:
            minimo = array[i]
        if array[i] > maximo:
            maximo = array[i]
    intervalWidth= (maximo-minimo)/numOfInterv
    intervalList=[minimo]
    for i in range(1, numOfInterv):
        intervalList.append(intervalList[i-1]+intervalWidth)
    return intervalList
 
def calculateExpectedFreq(array, type, parameters, totalSize):
    freq = []
    if (type == 'Exponential'):
        for i in range(1,len(array)):
            expRelativeFreq = (1-exp(-parameters*array[i])) - (1-exp(-parameters*array[i-1]))
            expAbsoluteFreq = expRelativeFreq*totalSize
            freq.append(expAbsoluteFreq)
    elif (type == 'Normal'):
        for i in range(1,len(array)):
            expRelativeFreq = (stats.norm.cdf(array[i],parameters[0],parameters[1])) - (stats.norm.cdf(array[i-1],parameters[0],parameters[1]))
            expAbsoluteFreq = expRelativeFreq*totalSize
            freq.append(expAbsoluteFreq)
    elif (type == 'Gamma'):
        loc = 0
        scale = 1/parameters[0]
        for i in range(1,len(array)):
            expRelativeFreq = (stats.gamma.cdf(array[i],parameters[1],loc,scale)) - (stats.gamma.cdf(array[i-1],parameters[1],loc,scale))
            expAbsoluteFreq = expRelativeFreq*totalSize
            freq.append(expAbsoluteFreq)
    rest = totalSize-sum(freq)
    freq.append(rest)
    return freq
 
def calcularIndice(numero, listIntervals):
    indice = 0
    for i in range(1,len(listIntervals)):
        if (numero >= listIntervals[i-1]) & (numero < listIntervals[i]):
            break
        elif (numero == listIntervals[len(listIntervals)-1]):
            indice = len(listIntervals)-2
            break
        else:
            indice += 1
    return indice
 
def getContinuousFreq(array, listIntervals):
    fr = []
    for i in range(len(listIntervals)):
        fr.append(0)
    for i in range(len(array)):
        fr[calcularIndice(array[i], listIntervals)] += 1
    return fr
 
def getDiscretFreq(array, max):
    size = len(array)
    fr = []
    for i in range(max+1):
        fr.append(0)
    for i in range(size):
        fr[array[i]] += 1
    return fr
 
def studyDiscrete(array, parameters, type):
    max = 0
    for i in range(len(array)):
        if (max < array[i]):
            max = array[i]
    obsFreq = getDiscretFreq(array, max)
    if (type == 'Empirical'):
        obsFreq.pop(0)
    else:
        plt.plot(obsFreq, "bo", color = 'g',label='Observado')
    x = range(0, max+1)
    expFreq = []
    if (type == 'Poisson'):
        for i in x:
            expFreq.append(functionPoisson(i, parameters)*size)
        plt.title(f"{type} (lambda= {parameters})") 
    elif (type == 'Empirical'):
        x = range(1, max+1)
        plt.plot(x,obsFreq, "bo", color = 'g',label='Observado')
        for i in parameters:
            expFreq.append(i*size)
        plt.title(f"{type}") 
    elif (type == 'Pascal'):
        for i in x:
            expFreq.append(functionPascal(i, parameters[0], parameters[1])*size)
        plt.title(f"{type} (k= {parameters[0]}  -  p= {parameters[1]})")
    elif (type == 'Hypergeometric'):
        for i in x:
            expFreq.append(functionHypergeometric(i, parameters[0], parameters[1], parameters[2])*size)
        plt.title(f"{type} (N= {parameters[0]}  -  n= {parameters[1]}  -  Np= {parameters[2]})")
    elif (type == 'Binomial'):
        for i in x:
            expFreq.append(functionBinomial(i, parameters[0], parameters[1])*size)
        plt.title(f"{type} (n= {parameters[0]}  -  p= {parameters[1]})")
    plt.plot(x, expFreq, color='k',label='Esperado')    
    plt.legend(loc="upper right")
    plt.show()
    testChiCuad(obsFreq, expFreq, 0.05)
 
def studyContinuous(array, parameters, type):
    size = len(array)
    listIntervals = calculateIntervals(array, numOfInterv)
    intervalWidth = (listIntervals[2]-listIntervals[1])/2
    listMiddle=[]
    for i in range(len(listIntervals)):
        listMiddle.append(listIntervals[i]+intervalWidth)
    frObservated = getContinuousFreq(array, listIntervals)
    plt.hist(array, 20, density=True, histtype="stepfilled", alpha=.7, linewidth=5, color='g', label='Observado')
    xmin, xmax = plt.xlim()
    if (type == "Exponential"):
        frExpected = calculateExpectedFreq(listIntervals, 'Exponential', alpha, size)
        x = linspace(0, xmax, 100)
        y = stats.expon.pdf(x, scale=1/parameters)
        plt.title(f"{type} (alpha= {parameters})") 
    if (type == "Normal"):
        frExpected = calculateExpectedFreq(listIntervals, 'Normal', [parameters[0], parameters[1]], size)
        x = linspace(xmin, xmax, 100)
        y = stats.norm.pdf(x, parameters[0], parameters[1])
        plt.title(f"{type} (mu= {parameters[0]}  -  sigma= {parameters[1]})")
    if (type == "Gamma"):
        frExpected = calculateExpectedFreq(listIntervals, 'Gamma', [parameters[0], parameters[1]], size)
        x = linspace(0, xmax, 100)
        y = stats.gamma.pdf(x, parameters[1], scale = 1/parameters[0])
        plt.title(f"{type} (alpha= {parameters[0]}  -  k= {parameters[1]})")
    plt.plot(x, y, color='k', label='Esperado')    
    plt.legend(loc="upper right")
    plt.show()
    testChiCuad(frObservated, frExpected, 0.05)
    print("")
    kstest(array, 0.05, parameters, type)
 
 
def functionPoisson(x, lambda_):
    return (lambda_**x)*exp(-lambda_)/factorial(x)
 
def functionBinomial(x, n, p):
    return combinatoria(n,x)*(p**x)*(1-p)**(n-x)
 
def functionHypergeometric(x, N, n, Np):
    return combinatoria(Np, x)*combinatoria(N-Np,n-x)/combinatoria(N,n)
 
def functionPascal(x, k, p):
    return combinatoria(k+x-1,x)*(p**k)*((1-p)**(x))
 
def combinatoria(a,b):
    return factorial(a)/(factorial(b)*factorial(a-b))
 
def factorial(num):
    factorial = 1
    if (num > 1):
        for i in range(1,num + 1):
            factorial = factorial*i
    return factorial
 
 
def menu():
    print(" ")    
    print ("MENÚ PRINCIPAL - Selecciona una opción")
    print ("1 - Uniforme")
    print ("2 - Exponencial")
    print ("3 - Gamma")
    print ("4 - Normal")
    print ("5 - Pascal")
    print ("6 - Binomial")
    print ("7 - Hipergeometrica")
    print ("8 - Poisson")
    print ("9 - Empirica Discreta")
    print ("0 - Salir")
 
size = 10000
numOfInterv = 1
r = 0
while r < size:
    numOfInterv += 1
    r = 2**numOfInterv
empiricalProb = [0.273, 0.037, 0.195, 0.009, 0.124, 0.058, 0.062, 0.151, 0.047, 0.044]
 
while True:
    os.system('cls')
    print("Tamaño de la sucesion:", size)
    menu()
   
    opcionMenu = input("Ingrese una opción ")
    print(" ")
    os.system('cls')
 
    if opcionMenu=="1":
        print("Distribucion Uniforme")
        a = float(input("Ingrese a:"))
        b = float(input("Ingrese b:"))
        uniform_list = uniform(size,a,b)
        listIntervals = calculateIntervals(uniform_list, numOfInterv)
        frExpected = []
        frExp = size/len(listIntervals)
        for i in range(len(listIntervals)):
            frExpected.append(frExp)
        frObservated = getContinuousFreq(uniform_list, listIntervals)
        plt.hlines(frExp, a, b, color='k',label='Esperado')
        plt.hist(uniform_list, numOfInterv,color='g',label='Observado',alpha=.7)
        plt.legend(loc="upper right")
        plt.title(f"Uniform (a= {a}  -  b= {b})")
        plt.axis([a-(b-a)/6,b+(b-a)/6, 0, size/(numOfInterv/2)])
        plt.show()
        testChiCuad(frObservated, frExpected, 0.05)
        print("")
        kstest(uniform_list, 0.05, [a,b], "Uniform")
        input("Pulsa una tecla para continuar")
   
    elif opcionMenu=="2":
        print("Distribucion Exponencial")
        alpha = float(input("Ingrese el valor de alpha:"))
        exponential_list = exponential(size, alpha)        
        studyContinuous(exponential_list, alpha, "Exponential")
        input("Pulsa una tecla para continuar")
 
    elif opcionMenu=="3":
        print("Distribucion Gamma")
        alpha = float(input("Ingrese el valor de alpha:"))
        k = int(input("Ingrese el valor de k:"))
        gamma_list = gamma(size, k, alpha)
        studyContinuous(gamma_list, [alpha,k], "Gamma")
        input("Pulsa una tecla para continuar")
 
    elif opcionMenu=="4":
        print("Distribucion Normal")
        mu = float(input("Ingrese el valor de mu:"))
        sigma = float(input("Ingrese el valor de sigma:"))
        normal_list = normal(size, mu, sigma)
        studyContinuous(normal_list, [mu, sigma], "Normal")
        input("Pulsa una tecla para continuar")
   
    elif opcionMenu=="5":
        print("Distribucion Pascal")
        k = int(input("Ingrese el valor de k:"))
        p = float(input("Ingrese el valor de p:"))
        pascal_list = pascal(size, k, p)
        studyDiscrete(pascal_list, [k, p], 'Pascal')
        input("Pulsa una tecla para continuar")
 
    elif opcionMenu=="6":
        print("Distribucion Binomial")
        n = int(input("Ingrese el valor de n:"))
        p = float(input("Ingrese el valor de p:"))
        binom_list = binom(size, n, p)
        studyDiscrete(binom_list, [n, p], 'Binomial')
        input("Pulsa una tecla para continuar")
 
    elif opcionMenu=="7":
        print("Distribucion Hipergeometrica")
        N = int(input("Ingrese el valor de N:"))
        n = int(input("Ingrese el valor de n:"))
        Np = int(input("Ingrese el valor de Np:"))
        hiperg_list = hypergeometric(size, N, n, Np)
        studyDiscrete(hiperg_list, [N, n, Np], 'Hypergeometric')
        input("Pulsa una tecla para continuar")
   
    elif opcionMenu=="8":
        print("Distribucion Poisson")
        lambda_ = float(input("Ingrese el valor de lambda:"))
        poisson_list = poisson(size, lambda_)
        studyDiscrete(poisson_list, lambda_, 'Poisson')
        input("Pulsa una tecla para continuar")
   
    elif opcionMenu=="9":
        print("Distribucion Empirica Discreta")
        empirical_list = empirical(size, empiricalProb)
        studyDiscrete(empirical_list, empiricalProb, 'Empirical')
        input("Pulsa una tecla para continuar")
 
    elif opcionMenu=="0":
        break
   
    else:
        print(" ")
        input("No has pulsado ninguna opción correcta... \n Pulsa una tecla para continuar")
    os.system('cls')