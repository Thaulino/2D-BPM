# coding=utf8
import math as ma
import cmath as cma
import pylab as pl
import numpy as np

def PadeApproximation(startAr ,nRef, freq, nMeshAr,xSize,zSize,xStep,zStep):
##      makes BPM with pade (1)
##      startAr : array consists of initial values for approx.
##      nRef : reference index often called nb
##      req : frequency for approx.
##      nMeshAr : array consists of refractive index for each mesh point (x*z big)
##      xSize : amount of mesh point in x direction
##      zSize : amount of mesh points in z direction
##      xStep :stepsize in x direction in m
##      zStep : stepsize in z direction in m
      
      ##### some constants
      #####internal function declaration

      def PadeOneStep(xAr):
            "calculates one step in x - direction of Pade approx."

            ##### some constants

            cVacuum =2.998
            cVacuum =cVacuum*ma.pow(10,8)
            #give wavelenght
            LVacuum = cVacuum/freq
            #k0- wavevector-absolute value
            kVacuum = 2*ma.pi/LVacuum
            #ref wavevector
            kRef=kVacuum*nRef

            d=kVacuum*zStep*nRef
            p=complex(1,d)
            p=p/4

            pC=p.conjugate()

            #alpha
            alpha= p/(ma.pow(kRef,2)*ma.pow(xStep,2))
            #konjugiert komp. alpha
            alphaC = pC/(ma.pow(kRef,2)*ma.pow(xStep,2))

            #Arrays for Constants
            zNowConstantAr=np.zeros((xSize,xSize), np.complex)
            #inverse of zNext... has to be calculated later
            zNextConstantAr=np.zeros((xSize,xSize), np.complex)

            ##### internal function declaration

            def CalcOmega(x,z):
                  res=ma.pow(nMeshAr[x,z],2)-ma.pow(nRef,2)
                  res=res*(ma.pow(kVacuum,2)/ma.pow(kRef,2))
                  return res*p+1-2*alpha

            #konjugiert komp. omega
            def CalcOmegaC(x,z):
                  res=ma.pow(nMeshAr[x,z],2)-ma.pow(nRef,2)
                  res=res*(ma.pow(kVacuum,2)/ma.pow(kRef,2))
                  return res*pC+1-2*alphaC

            

            ##### Fill Constant Arrays

            #Hauptdiagonale Omega/OmegaC und nur links oder rechts alpha/alphaC für x=0 oder x=999
            #x=999
            zNowConstantAr[xSize-1,xSize-1]=CalcOmegaC(xSize-1,zCurrent)+alphaC
            zNowConstantAr[xSize-1,xSize-2]=alphaC
            zNextConstantAr[xSize-1,xSize-1]=CalcOmega(xSize-1,zCurrent)+alpha
            zNextConstantAr[xSize-1,xSize-2]=alpha
            #x=0
            zNowConstantAr[0,0]=CalcOmegaC(0,zCurrent)+alphaC
            zNowConstantAr[0,1]=alphaC
            zNextConstantAr[0,0]=CalcOmega(0,zCurrent)+alpha
            zNextConstantAr[0,1]=alpha
            #Hauptdiagonale mit Omega/OmegaC, links und rechts alpha/alphaC für 0<x<999
            xCurrent=1
            while xCurrent<xSize-1:
                  zNowConstantAr[xCurrent,xCurrent]=CalcOmegaC(xCurrent,zCurrent)
                  zNowConstantAr[xCurrent,xCurrent+1]=alphaC
                  zNowConstantAr[xCurrent,xCurrent-1]=alphaC
                  zNextConstantAr[xCurrent,xCurrent]=CalcOmega(xCurrent,zCurrent)
                  zNextConstantAr[xCurrent,xCurrent+1]=alpha
                  zNextConstantAr[xCurrent,xCurrent-1]=alpha
                  xCurrent=xCurrent+1
                  
            #Matrixmultiplikation um nächste Elemente zu berechnen
            zNextConstantInverseAr=np.linalg.inv(zNextConstantAr)
            multiConstantAr=np.dot(zNextConstantInverseAr,zNowConstantAr)
            resultStepAr=np.dot(multiConstantAr,xAr)
            return np.matrix(resultStepAr[:])


############
      
      #shows content of nMeshAr
      pl.figure(217)
      matrix=np.matrix(nMeshAr)
      pl.imshow(matrix)
      pl.colorbar()
      pl.draw()
      pl.title("Meshed area")
      pl.show()
      
      
      zCurrent=0
      #Array for saving magnetic Field Hz Amplitude
      resultBPM = np.zeros((xSize,zSize), np.complex)
      resultBPM[:,0]=startAr
      
      #loop wird zSize-2 mal durchlaufen, weil eine Ebene als Start gegeben ist
      while zCurrent<zSize-1:
            nextStepResult=PadeOneStep(resultBPM[:,zCurrent])
            i=0
            while i<xSize:
                  resultBPM[i,zCurrent+1]=complex((nextStepResult[0,i].real),(nextStepResult[0,i].imag))
                  i=i+1
            if ma.fmod(zCurrent,50) == 0:
                  print("-----> running")
##                 
            zCurrent=zCurrent+1

            

      #create array with absolute value of BPM result
      resultAbBPM=np.zeros((xSize,zSize))
      i=0
      while i<xSize:
            ii=0
            while ii<zSize:
                  resultAbBPM[i,ii]= abs(resultBPM[i,ii])
                  ii=ii+1
            i=i+1
     

      #  
      pl.figure(1)
      matrix=np.matrix(resultAbBPM[:,:])
      pl.imshow(matrix)
      pl.colorbar()
      pl.title('magnitude distribution')
      pl.draw()
      pl.show()
      

      pl.figure(2)
      x_axis=[]
      y_axis=[]
      xCurrent=0
      while xCurrent<xSize:
            x_axis.append(xCurrent)
            y_axis.append(resultAbBPM[xCurrent,zCurrent])
            xCurrent=xCurrent+1

      pl.title('magnitude distribution at output')
      pl.plot(x_axis,y_axis)
      pl.draw()
      pl.show()

        

      return 0









