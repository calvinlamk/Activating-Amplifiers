from __future__ import division
from PySteppables import *       
import math                      
import numpy as np                   
import CompuCell
import sys
import random
RNG=random.SystemRandom()

#Motility Variables
CtoM=52  
BASAL=100 
SCF=0.5 

#Self-Cutoff
ENDMCS=50000

#Mitosis Variables
RADAVG=3
RADDEV=.5
MTFORCEMIN=-3*10**(-3.88) 
MTFORCEMAX=4*10**(-3.88) 

#Signaling Variables
CONEXPSCF=10000

BETA=1000
BETAAMP=1000
KAPPA=25000

THRESHOLD=7000 

#Sampling and Comp Speed
RESOL=100 
USEDNODES=8 


class ELUGMSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):

        self.pW1 = self.addNewPlotWindow(
            _title='Calibrate',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Phi',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW1.addPlot('BWtoYG', _style='Dots', _color='gray', _size=3) 
        self.pW1.addPlot('BWtoBW', _style='Dots', _color='cyan', _size=3)
        self.pW1.addPlot('BWtoOL', _style='Dots', _color='orange', _size=3)
        self.pW1.addPlot('BWtoRK', _style='Dots', _color='red', _size=3)
        self.pW1.addPlot('RKtoYG', _style='Dots', _color='gray', _size=5)
        self.pW1.addPlot('RKtoBW', _style='Dots', _color='cyan', _size=5)        
        self.pW1.addPlot('RKtoOL', _style='Dots', _color='orange', _size=5) 
        self.pW1.addPlot('RKtoRK', _style='Dots', _color='red', _size=5)
        
        self.pW2 = self.addNewPlotWindow(
            _title='Psi',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Psi',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW2.addPlot('Hamiltonian', _style='Dots', _color='white', _size=3)
 
        self.pW3 = self.addNewPlotWindow(
            _title='Types',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW3.addPlot('Y', _style='Dots', _color='gray', _size=3)
        self.pW3.addPlot('G', _style='Dots', _color='green', _size=3)
        self.pW3.addPlot('B', _style='Dots', _color='blue', _size=3)   
        self.pW3.addPlot('R', _style='Dots', _color='red', _size=3)
        self.pW3.addPlot('O', _style='Dots', _color='orange', _size=3)
        self.pW3.addPlot('W', _style='Dots', _color='cyan', _size=3) 
        self.pW3.addPlot('L', _style='Dots', _color='yellow', _size=3)
        self.pW3.addPlot('K', _style='Dots', _color='pink', _size=3)        
        
        self.pW4 = self.addNewPlotWindow(
            _title='Point System',
            _xAxisTitle='MonteCarlo Step (MCS)',
            _yAxisTitle='Count',
            _xScaleType='linear',
            _yScaleType='linear',
            _grid=True)                
        self.pW4.addPlot('GFP', _style='Dots', _color='green', _size=3)
        self.pW4.addPlot('BFP', _style='Dots', _color='blue', _size=3) 
        self.pW4.addPlot('YFP', _style='Dots', _color='yellow', _size=3)
        self.pW4.addPlot('RFP', _style='Dots', _color='red', _size=3)     
                
        global YtoY,YtoG,GtoY,YtoB,BtoY,YtoR,RtoY,GtoG,GtoB,BtoG,GtoR,RtoG,BtoB,BtoR,RtoB,RtoR,OtoY,YtoO,OtoG,GtoO,OtoB,BtoO,OtoR,RtoO,OtoO,WtoY,YtoW,WtoG,GtoW,WtoB,BtoW,WtoR,RtoW,WtoO,OtoW,WtoW  #the adhesion matrix, call these values and store them for motility code
        global LtoY,YtoL,LtoG,GtoL,LtoB,BtoL,LtoR,RtoL,LtoO,OtoL,LtoW,WtoL,LtoL,KtoY,YtoK,KtoG,GtoK,KtoB,BtoK,KtoR,RtoK,KtoO,OtoK,KtoW,WtoK,KtoL,LtoK,KtoK
        YtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','Y']))
        YtoG=GtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','G']))
        YtoB=BtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','B']))
        YtoR=RtoY=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','Y','Type2','R']))
        GtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','G']))
        GtoB=BtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','B']))
        GtoR=RtoG=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','G','Type2','R']))
        BtoB=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','B','Type2','B']))
        BtoR=RtoB=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','B','Type2','R']))
        RtoR=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','R','Type2','R']))
        OtoY=YtoO=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','O','Type2','Y']))    
        OtoG=GtoO=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','O','Type2','G']))
        OtoB=BtoO=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','O','Type2','B']))
        OtoR=RtoO=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','O','Type2','R']))
        OtoO=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','O','Type2','O']))
        WtoY=YtoW=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','W','Type2','Y']))
        WtoG=GtoW=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','W','Type2','G']))
        WtoB=BtoW=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','W','Type2','B']))
        WtoR=RtoW=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','W','Type2','R']))
        WtoO=OtoW=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','W','Type2','O']))
        WtoW=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','W','Type2','W']))
        
        LtoY=YtoL=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','L','Type2','Y']))
        LtoG=GtoL=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','L','Type2','G']))
        LtoB=BtoL=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','L','Type2','B']))
        LtoR=RtoL=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','L','Type2','R']))
        LtoO=OtoL=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','L','Type2','O']))
        LtoW=WtoL=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','L','Type2','W']))
        LtoL=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','L','Type2','L']))
        
        KtoY=YtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','Y']))
        KtoG=GtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','G']))
        KtoB=BtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','B']))
        KtoR=RtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','R']))
        KtoO=OtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','O']))
        KtoW=WtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','W']))
        KtoL=LtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','L']))
        KtoK=float(self.getXMLElementValue(['Plugin','Name','Contact'],['Energy','Type1','K','Type2','K']))        

        for cell in self.cellList:
            cell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV)
            cell.lambdaSurface=2.2                   
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2 
            cell.lambdaVolume=2.2                    
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3
            cell.dict["PTS"]=[0,0,0,0] #order is green, red, blue, yellow fluroscent protein

    def step(self,mcs):                 
        NUMTY=0 
        NUMTG=0 
        NUMTB=0 
        NUMTR=0 
        NUMTO=0
        NUMTW=0
        NUMTL=0
        NUMTK=0
        
        GFPPTS=0
        RFPPTS=0 
        BFPPTS=0
        YFPPTS=0
        
        SYPSI=0

        CSABWYG=0
        CSABWBW=0
        CSABWOL=0
        CSABWRK=0
        
        CSARKYG=0
        CSARKBW=0
        CSARKOL=0
        CSARKRK=0
        
        BWCYG=0
        BWCBW=0
        BWCOL=0
        BWCRK=0

        RKCYG=0
        RKCBW=0
        RKCOL=0
        RKCRK=0
        
        if mcs==1:
            self.changeNumberOfWorkNodes(USEDNODES) #set to necessary computational nodes                  

        for cell in self.cellList:
            CSAY=0
            CSAG=0
            CSAB=0 
            CSAR=0
            CSAM=0 
            CSAO=0
            CSAW=0
            CSAL=0
            CSAK=0
            
            PTSY=0 
            PTSG=0 
            PTSB=0
            PTSR=0 
            PTSO=0
            PTSW=0
            PTSL=0
            PTSK=0
            DTREPO=0 #delta reporter One
            DTREPT=0 #deta reporter Two
            DTREPH=0 #delta reporter three
            DTREPF=0 #reporter four

            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor is None: #none type refers to medium 
                    continue
                if neighbor.type==1: #gray cells
                    CSAY+=commonSurfaceArea
                if neighbor.type==2: #green cells
                    CSAG+=commonSurfaceArea
                if neighbor.type==3: #blue cells
                    CSAB+=commonSurfaceArea
                    PTSB+=commonSurfaceArea*neighbor.dict["PTS"][2]/(neighbor.surface)
                if neighbor.type==6: #cyan
                    CSAW+=commonSurfaceArea
                    PTSW+=commonSurfaceArea*neighbor.dict["PTS"][2]/(neighbor.surface)                                                 
                if neighbor.type==5: #orange cells
                    CSAO+=commonSurfaceArea
                    PTSO+=commonSurfaceArea*CONEXPSCF/neighbor.surface
                if neighbor.type==7: #yellow cells
                    CSAL+=commonSurfaceArea
                    PTSL+=commonSurfaceArea*CONEXPSCF/neighbor.surface
                if neighbor.type==4: #red cells
                    CSAR+=commonSurfaceArea
                    PTSR+=commonSurfaceArea*CONEXPSCF/neighbor.surface 
                if neighbor.type==8: #pink cells
                    CSAK+=commonSurfaceArea
                    PTSK+=commonSurfaceArea*CONEXPSCF/neighbor.surface                    
            CSAM=cell.surface-(CSAY+CSAG+CSAB+CSAR+CSAO+CSAW+CSAL+CSAK) #alternative method to calculate common surface area with medium                

            if (cell.type==1 or cell.type==2 or cell.type==3 or cell.type==6):
                DTREPO=(1/(1+np.exp(-((PTSO+PTSL+PTSR+PTSK)-BETA))))-(1/KAPPA)*cell.dict["PTS"][0]         #GFP       GRBY 0123
                cell.dict["PTS"][0]+=DTREPO                                                          
                DTREPT=(1/(1+np.exp(-((cell.dict["PTS"][0])-BETA))))-(1/KAPPA)*cell.dict["PTS"][2]    #BFP 
                cell.dict["PTS"][2]+=DTREPT                                                                                                      
            if (cell.type==5 or cell.type==4 or cell.type==7 or cell.type==8):
                DTREPH=(1/(1+np.exp(-((PTSB+PTSW)-BETA))))-(1/KAPPA)*cell.dict["PTS"][3]           #YFP  
                cell.dict["PTS"][3]+=DTREPH
                DTREPF=(1/(1+np.exp(-((cell.dict["PTS"][3])-BETA))))-(1/KAPPA)*cell.dict["PTS"][1]           #RFP mcherry      
                cell.dict["PTS"][1]+=DTREPF                
                                                             
            if (cell.type==1 or cell.type==2 or cell.type==3 or cell.type==6): #FACS PLOTS AND CADHERIN LINAKGE TO STATE
                if cell.dict["PTS"][0]>=THRESHOLD and cell.dict["PTS"][2]<THRESHOLD:   #GFP+ BFP- state
                    cell.type=2
                if cell.dict["PTS"][2]>=THRESHOLD and cell.dict["PTS"][0]<THRESHOLD:   #GFP- BFP+ state
                    cell.type=3                 
                if cell.dict["PTS"][0]>=THRESHOLD and cell.dict["PTS"][2]>=THRESHOLD:   #GFP+ BFP+ state
                    cell.type=6
                if cell.dict["PTS"][0]<THRESHOLD and cell.dict["PTS"][2]<THRESHOLD:   #GFP- BFP- state gray color
                    cell.type=1
            if (cell.type==5 or cell.type==4 or cell.type==7 or cell.type==8):
                if cell.dict["PTS"][3]>=THRESHOLD and cell.dict["PTS"][1]<THRESHOLD:   #YFP+ RFP-
                    cell.type=7
                if cell.dict["PTS"][1]>=THRESHOLD and cell.dict["PTS"][3]<THRESHOLD:   #YFP- RFP+ state
                    cell.type=4                 
                if cell.dict["PTS"][1]>=THRESHOLD and cell.dict["PTS"][3]>=THRESHOLD:   #YFP+ RFP+ state
                    cell.type=8
                if cell.dict["PTS"][1]<THRESHOLD and cell.dict["PTS"][3]<THRESHOLD:   #YFP- RFP- state orange color
                    cell.type=5
                    
            if cell.type==1: #gray cells                              
                cell.lambdaSurface=2.2           
                cell.lambdaVolume=2.2          
                NUMTY+=1                         
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+YtoY*CSAY+YtoG*CSAG+YtoB*CSAB+YtoR*CSAR+YtoO*CSAO+YtoW*CSAW+YtoL*CSAL+YtoK*CSAK)/cell.surface
                GFPPTS+=cell.dict["PTS"][0]
                BFPPTS+=cell.dict["PTS"][2]

            if cell.type==2: #green cells
                cell.lambdaSurface=2.2            
                cell.lambdaVolume=2.2                
                NUMTG+=1                          
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+GtoY*CSAY+GtoG*CSAG+GtoB*CSAB+GtoR*CSAR+GtoO*CSAO+GtoW*CSAW+GtoL*CSAL+GtoK*CSAK)/cell.surface
                GFPPTS+=cell.dict["PTS"][0]
                BFPPTS+=cell.dict["PTS"][2]
             
            if cell.type==3: #blue cells            
                CSABWYG+=(CSAY+CSAG)/cell.surface
                if (CSAY+CSAG)>0:
                    BWCYG+=1                
                CSABWBW+=(CSAB+CSAW)/cell.surface
                if (CSAB+CSAW)>0:
                    BWCBW+=1 
                CSABWOL+=(CSAO+CSAL)/cell.surface
                if (CSAO+CSAL)>0:
                    BWCOL+=1 
                CSABWRK+=(CSAR+CSAK)/cell.surface
                if (CSAR+CSAK)>0:
                    BWCRK+=1
                cell.lambdaSurface=1.0           
                cell.lambdaVolume=1.0              
                NUMTB+=1                         
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+BtoY*CSAY+BtoG*CSAG+BtoB*CSAB+BtoR*CSAR+BtoO*CSAO+BtoW*CSAW+BtoL*CSAL+BtoK*CSAK)/cell.surface 
                GFPPTS+=cell.dict["PTS"][0]
                BFPPTS+=cell.dict["PTS"][2]
                
            if cell.type==6: #syan
                CSABWYG+=(CSAY+CSAG)/cell.surface
                if (CSAY+CSAG)>0:
                    BWCYG+=1                
                CSABWBW+=(CSAB+CSAW)/cell.surface
                if (CSAB+CSAW)>0:
                    BWCBW+=1 
                CSABWOL+=(CSAO+CSAL)/cell.surface
                if (CSAO+CSAL)>0:
                    BWCOL+=1 
                CSABWRK+=(CSAR+CSAK)/cell.surface
                if (CSAR+CSAK)>0:
                    BWCRK+=1      
                cell.lambdaSurface=1.0            
                cell.lambdaVolume=1.0                 
                NUMTW+=1
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+WtoY*CSAY+WtoG*CSAG+WtoB*CSAB+WtoR*CSAR+WtoO*CSAO+WtoW*CSAW+WtoL*CSAL+WtoK*CSAK)/cell.surface #corrected cell motility, vetted    
                GFPPTS+=cell.dict["PTS"][0]
                BFPPTS+=cell.dict["PTS"][2]

            if cell.type==5: #orange cells                                                     
                cell.lambdaSurface=2.2            
                cell.lambdaVolume=2.2                   
                NUMTO+=1
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+OtoY*CSAY+OtoG*CSAG+OtoB*CSAB+OtoR*CSAR+OtoO*CSAO+OtoW*CSAW+OtoL*CSAL+OtoK*CSAK)/cell.surface #corrected cell motility, vetted           
                YFPPTS+=cell.dict["PTS"][3]
                RFPPTS+=cell.dict["PTS"][1] 

            if cell.type==7: #yellow cells                                                    
                cell.lambdaSurface=2.2           
                cell.lambdaVolume=2.2                 
                NUMTL+=1
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+LtoY*CSAY+LtoG*CSAG+LtoB*CSAB+LtoR*CSAR+LtoO*CSAO+LtoW*CSAW+LtoL*CSAL+LtoK*CSAK)/cell.surface #corrected cell motility, vetted           
                YFPPTS+=cell.dict["PTS"][3]
                RFPPTS+=cell.dict["PTS"][1] 

            if cell.type==4: #red cells                       
                CSARKYG+=(CSAY+CSAG)/cell.surface
                if (CSAY+CSAG)>0:
                    RKCYG+=1                
                CSARKBW+=(CSAB+CSAW)/cell.surface
                if (CSAB+CSAW)>0:
                    RKCBW+=1 
                CSARKOL+=(CSAO+CSAL)/cell.surface
                if (CSAO+CSAL)>0:
                    RKCOL+=1 
                CSARKRK+=(CSAR+CSAK)/cell.surface
                if (CSAR+CSAK)>0:
                    RKCRK+=1 
                cell.lambdaSurface=1.0            
                cell.lambdaVolume=1.0                  
                NUMTR+=1                          
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+RtoY*CSAY+RtoG*CSAG+RtoB*CSAB+RtoR*CSAR+RtoO*CSAO+RtoW*CSAW+RtoL*CSAL+RtoK*CSAK)/cell.surface #corrected cell motility, vetted
                YFPPTS+=cell.dict["PTS"][3]
                RFPPTS+=cell.dict["PTS"][1]

            if cell.type==8: #pink cells                         
                CSARKYG+=(CSAY+CSAG)/cell.surface
                if (CSAY+CSAG)>0:
                    RKCYG+=1                
                CSARKBW+=(CSAB+CSAW)/cell.surface
                if (CSAB+CSAW)>0:
                    RKCBW+=1 
                CSARKOL+=(CSAO+CSAL)/cell.surface
                if (CSAO+CSAL)>0:
                    RKCOL+=1 
                CSARKRK+=(CSAR+CSAK)/cell.surface
                if (CSAR+CSAK)>0:
                    RKCRK+=1 
                cell.lambdaSurface=1.0            
                cell.lambdaVolume=1.0                 
                NUMTK+=1                         
                cell.fluctAmpl=BASAL+SCF*(CtoM*CSAM+KtoY*CSAY+KtoG*CSAG+KtoB*CSAB+KtoR*CSAR+KtoO*CSAO+KtoW*CSAW+KtoL*CSAL+KtoK*CSAK)/cell.surface #corrected cell motility, vetted
                YFPPTS+=cell.dict["PTS"][3]
                RFPPTS+=cell.dict["PTS"][1] 
    
            SYPSI+=cell.fluctAmpl

        if mcs%RESOL==0:
            
            self.pW2.addDataPoint("Hamiltonian", mcs, SYPSI)
            
            self.pW3.addDataPoint("Y", mcs, NUMTY)
            self.pW3.addDataPoint("G", mcs, NUMTG) 
            self.pW3.addDataPoint("B", mcs, NUMTB)   
            self.pW3.addDataPoint("R", mcs, NUMTR)
            self.pW3.addDataPoint("O", mcs, NUMTO)
            self.pW3.addDataPoint("W", mcs, NUMTW)
            self.pW3.addDataPoint("L", mcs, NUMTL)
            self.pW3.addDataPoint("K", mcs, NUMTK)

            self.pW4.addDataPoint("GFP", mcs, GFPPTS/(NUMTY+NUMTG+NUMTB+NUMTW))
            self.pW4.addDataPoint("BFP", mcs, BFPPTS/(NUMTY+NUMTG+NUMTB+NUMTW))
            self.pW4.addDataPoint("YFP", mcs, YFPPTS/(NUMTO+NUMTL+NUMTR+NUMTK))
            self.pW4.addDataPoint("RFP", mcs, RFPPTS/(NUMTO+NUMTL+NUMTR+NUMTK))
            
            if BWCYG==0:
                self.pW1.addDataPoint("BWtoYG", mcs, 0)
            if BWCYG>0:
                self.pW1.addDataPoint("BWtoYG", mcs, CSABWYG/BWCYG)
            if BWCBW==0:
                self.pW1.addDataPoint("BWtoBW", mcs, 0)
            if BWCBW>0:
                self.pW1.addDataPoint("BWtoBW", mcs, CSABWBW/BWCBW)
            if BWCOL==0:
                self.pW1.addDataPoint("BWtoOL", mcs, 0)
            if BWCOL>0:
                self.pW1.addDataPoint("BWtoOL", mcs, CSABWOL/BWCOL)
            if BWCRK==0:
                self.pW1.addDataPoint("BWtoRK", mcs, 0)
            if BWCRK>0:
                self.pW1.addDataPoint("BWtoRK", mcs, CSABWRK/BWCRK)
                
            if RKCYG==0:
                self.pW1.addDataPoint("RKtoYG", mcs, 0)
            if RKCYG>0:
                self.pW1.addDataPoint("RKtoYG", mcs, CSARKYG/RKCYG)
            if RKCBW==0:
                self.pW1.addDataPoint("RKtoBW", mcs, 0)
            if RKCBW>0:
                self.pW1.addDataPoint("RKtoBW", mcs, CSARKBW/RKCBW)
            if RKCOL==0:
                self.pW1.addDataPoint("RKtoOL", mcs, 0)
            if RKCOL>0:
                self.pW1.addDataPoint("RKtoOL", mcs, CSARKOL/RKCOL)
            if RKCRK==0:
                self.pW1.addDataPoint("RKtoRK", mcs, 0)
            if RKCRK>0:
                self.pW1.addDataPoint("RKtoRK", mcs, CSARKRK/RKCRK)
             
            if mcs==ENDMCS:
                fileName = "Calibrate" + str(mcs) + ".txt"
                self.pW1.savePlotAsData(fileName)                
                fileName = "PSI" + str(mcs) + ".txt"
                self.pW2.savePlotAsData(fileName)
                fileName = "FOU" + str(mcs) + ".txt"
                self.pW3.savePlotAsData(fileName)
                fileName = "SIG" + str(mcs) + ".txt"
                self.pW4.savePlotAsData(fileName)                 
                self.stopSimulation()                 
    def finish(self):
        pass

class ExtraFields(SteppableBasePy):
  def __init__(self,_simulator,_frequency=1):
    SteppableBasePy.__init__(self,_simulator,_frequency)
    
    self.scalarFieldGFP=self.createScalarFieldCellLevelPy("GFP")
    self.scalarFieldmCherry=self.createScalarFieldCellLevelPy("mCherry")
    self.scalarFieldBFP=self.createScalarFieldCellLevelPy("BFP")
    self.scalarFieldYFP=self.createScalarFieldCellLevelPy("YFP")

  def step(self,mcs):
    self.scalarFieldGFP.clear()
    self.scalarFieldmCherry.clear()
    self.scalarFieldBFP.clear()
    self.scalarFieldYFP.clear()
    
    for cell in self.cellList:       
        self.scalarFieldGFP[cell]=cell.dict["PTS"][0]/KAPPA
        self.scalarFieldmCherry[cell]=cell.dict["PTS"][1]/KAPPA
        self.scalarFieldBFP[cell]=cell.dict["PTS"][2]/KAPPA
        self.scalarFieldYFP[cell]=cell.dict["PTS"][3]/KAPPA
        
from PySteppablesExamples import MitosisSteppableBase

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
        self.setParentChildPositionFlag(0) #randomize child cell position, see developer manual
    def step(self,mcs):        
        cells_to_divide=[]          #gen cells to divide list
        for cell in self.cellList:
            cell.dict["RDM"]+=RNG.uniform(MTFORCEMIN,MTFORCEMAX) #make cells grow in target radius by this much
            cell.targetSurface=4*math.pi*cell.dict["RDM"]**2 #spherical surface area
            cell.targetVolume=(4/3)*math.pi*cell.dict["RDM"]**3 #spherical volume
            if cell.volume>2*(4/3)*math.pi*RADAVG**3: #divide at two times the mean radius initialized with               
                cells_to_divide.append(cell)           #add these cells to divide list
                
        for cell in cells_to_divide:
            self.divideCellRandomOrientation(cell)  #divide the cells

    def updateAttributes(self):
        self.parentCell.dict["RDM"]=RNG.gauss(RADAVG,RADDEV) #reassign new target radius
        self.parentCell.targetVolume=(4/3)*math.pi*self.parentCell.dict["RDM"]**3 #new target volume
        self.parentCell.targetSurface=4*math.pi*self.parentCell.dict["RDM"]**2 #new target surface area
        self.cloneParent2Child()  #copy characterstics to child cell, indlucig signaling
        