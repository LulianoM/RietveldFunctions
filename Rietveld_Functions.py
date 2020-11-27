import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
from scipy.optimize import curve_fit,minimize
import xlsxwriter


def remove_duplicate_words(string):
    return ' '.join(set(string.split()))

def rename(file):
    file = file.upper().replace('_',' ').replace('/',' ').replace('.',' ').replace('HOME LUCIANOMARTINS DOCUMENTOS NIDF EXPERIMENTO VATERITA-CALCITA DRX','')
    file = file.replace('TXT','')
    file = remove_duplicate_words(file)
    return(file)

def search_position(filename):
    with open(filename) as f:
        header = 1
        for line in f:
            if "[Data]" in line:
                return(header)
            else:
                header+=1
def plot_difratrogram(file_location):
    for file in glob.glob(file_location):
        plt.title('Diffractogram of '+ rename(file))
        position = search_position(file)
        file = pd.read_csv(file, sep=",", header=122)
        df = pd.DataFrame(file, columns = ['     Angle', '       PSD'])

        xdata = file['     Angle'].values
        ydata = file['       PSD'].values
        
        plt.plot(xdata, ydata, color='black')
        plt.xlabel('2$\\theta$')
        plt.ylabel('Intensity (a.u.)')
        fig = plt.figure(figsize=(15/2.54,12/2.54))
        plt.rc('font', size=11)
        plt.show()
def G(x,x0,sigma):
    fac = 1.0/np.sqrt(2*np.pi)/sigma
    return fac*np.exp(-np.power((x-x0)/sigma,2.0))

def L(x,x0,Gamma):
    fac = Gamma/np.pi
    return fac*(1.0/(np.power((x-x0),2.0)+np.power(Gamma,2.0)))

def B(x,b0,b1,b2):
    return b2*x**2+b1*x+b0

def pseudoVoigt(x,x0,eta,Gamma):
    sigmaGaussian = Gamma/np.sqrt(2*np.log(2))
    return ((1-eta)*G(x,x0,sigmaGaussian)+eta*L(x,x0,Gamma))

def plot_rietveld(file_location):
    for filename in glob.glob(file_location):
        
        with open(filename) as f:
            with open("read_file.txt","+w") as new_f:
                header=search_position(filename)
                #print(header)
                header_now=1
                for line in f:
                    if header_now>header:
                        new_f.write(line)
                        header_now+=1
                    else: 
                        header_now+=1   
                        
        data = pd.read_csv("read_file.txt", sep=",", header=0) #122
        #print(data)
        xdata = data['     Angle'].values
        ydata = data['       PSD'].values

        ####### Background Fitting 

        xslicedB = xdata[0:100].tolist()+xdata[1542:1750].tolist()+xdata[-100:-1].tolist()
        yslicedB = ydata[0:100].tolist()+ydata[1542:1750].tolist()+ydata[-100:-1].tolist()

        poptB, pcovB = curve_fit(B, xslicedB, yslicedB, p0=(min(yslicedB),0.0,0.0), ftol=1e-6)
        perrB = np.sqrt(np.diag(pcovB))

        b0 = poptB[0]
        b1 = poptB[1]
        b2 = poptB[2]

        # Fit function
        def DRXFit(x, I, x0, eta, Gamma):
            return (B(x,b0,b1,b2) + I*pseudoVoigt(x,x0,eta,Gamma))


        ############## Calcite peak 

        x0C = 29.5

        xslicedC = xdata[425:524]
        yslicedC = ydata[425:524]

        poptC, pcovC = curve_fit(DRXFit, xslicedC, yslicedC, p0=(0.9*max(yslicedC),x0C,0.5,0.5),bounds=([0.0,x0C-0.5,0.0,0.0],[1.0*max(yslicedC),x0C+0.5,1.0,2.0]), ftol=1e-6)
        perrC = np.sqrt(np.diag(pcovC))

        IC = poptC[0]


        ############## Vaterite peak 
        x0V = 25.0

        xslicedV = xdata[201:300]
        yslicedV = ydata[201:300]

        poptV, pcovV = curve_fit(DRXFit, xslicedV, yslicedV, p0=(0.9*max(yslicedV),x0V,0.5,0.5),bounds=([0.0,x0V-0.5,0.0,0.0],[1.0*max(yslicedV),x0V+0.5,1.0,2.0]), ftol=1e-6)
        perrV = np.sqrt(np.diag(pcovV))

        IV = poptV[0]


        ############## Aragonite peak 
        x0A = 45.9

        xslicedA = xdata[1239:1339]
        yslicedA = ydata[1239:1339]

        poptA, pcovA = curve_fit(DRXFit, xslicedA, yslicedA, p0=(0.9*max(yslicedA),x0A,0.5,0.5), bounds=([0.0,x0A-0.5,0.0,0.0],[1.0*max(yslicedA),x0A+0.5,1.0,2.0]), ftol=1e-6)
        perrA = np.sqrt(np.diag(pcovA))

        IA = poptA[0]


        #### Fraction of each Species

        # aV = 4.0; aA=3.5 (Original from Dickison)

        aV = 9.69
        aA = 0.97

        deno = IC + aV*IV + aA*IA

        xC = IC/deno
        xV = aV*IV/deno
        xA = aA*IA/deno
        
        print('Rietveld of '+ rename(filename))

        fig = plt.figure(figsize=(15/2.54,12/2.54))
        plt.rc('font', size=11)
        plt.plot(xdata,ydata, '.', label='data')
        plt.plot(xdata, B(xdata,*poptB), label='backg.') 
        plt.plot(xslicedC, DRXFit(xslicedC,*poptC), label='calcite') 
        plt.plot(xslicedV, DRXFit(xslicedV,*poptV), label='vaterite') 
        plt.plot(xslicedA, DRXFit(xslicedA,*poptA), label='aragonite') 
        plt.text(60, 0.5*max(ydata), 'xC = %1.2f\nxV = %1.2f\nxA = %1.2f' % (xC*100,xV*100,xA*100), horizontalalignment='center',
             verticalalignment='center')
        plt.legend()
        plt.xlabel('2$\\theta$')
        plt.ylabel('Intensity (a.u.)')
        plt.show()

