# -*- coding: utf-8 -*-
"""
Created on Thu May  2 19:45:14 2013

@author: rob
"""
from pylab import *
from scipy.optimize import leastsq

import unittest

class TestSFRTFunctions(unittest.TestCase):

    def setUp(self):
        pass
    def testTotal_Input_Impedance(self):
        rl=ones(len(w))*50.0
        L1,L2=[5.460742789470963, 1.8772183303814558]
        C1,C2=[0.53078597959243212, 1.4523709082435476]
        R=1.00021453518
        s=1j*w
        ztrue=R+L2*s+1.0/(C2*s+1.0/(L1*s+1/(1/rl+C1*s)))
        s0=s[0]
        abcdtrue=matrix([[1,-R],[0,1]])*matrix([[1,-L2*s0],[0,1]])*matrix([[1,0],[-C2*s0,1]])*matrix([[1,-L1*s0],[0,1]])*matrix([[1,0],[-C1*s0,1]])*matrix([[1,0],[-1/50.0,1]])
        abcdz=-abcdtrue[0,0]/abcdtrue[1,0]
        ztot=Total_Input_Impedance(w,-1,Ls,Cs,rl,R)            
        print "ABCD matrix calc: ",abcdz
        print "Function calc: ",ztot[0]
        print "Analytic calc: ",ztrue[0]
        self.assertTrue(ztot[0]==ztrue[0])
        self.assertTrue(abcdz==ztrue[0])
        
def par(Z1,Z2):
    Zpar=1.0/(1.0/Z1+1.0/Z2)
    return(Zpar)
    
    
def short_monopole_antenna_Z(omega,Rl=50,f0=100e6):
    omega=omega*f0*2*pi    
    L1,L2,L3,L4=0.225e-6,0.4e-6,0.01e-6,0.1e-6
    C3,C4=11.4e-12,19.8e-12
    s=omega*1j
    Zin=zeros(len(omega),dtype='complex64')
    Zin=L4*s+1/(C4*s)+par(L3*s+1/(C3*s),L2*s+par(L1*s,Rl))
    return(Zin)

def example4_Z(omega):
    s=1j*omega    
    L,C,R=0.607,4,1    
    z=L*s+1/(C*s+1)
    return(z)

    
def reflection_from_Z(w,zfunc,Z0):
    Zin=zfunc(w)    
    y0=1.0/Z0    
    s=(1.0-Zin)/(Zin+1.0)
    return(s)
    

    
#def reflection(w):
#    magfile,phasefile="mag.txt","phase.txt"
#    omega=2*pi*loadtxt(phasefile,delimiter=',')[:,0]
#    omega=omega/(0.5*amax(omega))
#    zl1=loadtxt(magfile,delimiter=',')[:,1]*exp(1j*loadtxt(phasefile,delimiter=',')[:,1]*pi/180.0)    
#    zl=interp(w,omega,real(zl1))+1j*interp(w,omega,imag(zl1))
#    sl=(zl-1)/(zl+1)
#    return(sl)
#    
#def reflection3(w):
#    p=1j*w
#    YS=4*p
#    YL=1+YS
#    ZL=0.6078*p+1.0/YL
#    SL=(ZL-1)/(ZL+1)
#    return(SL)



#def reflection2(w):
#    #This is specific to the short dipole antenna network
#    #in chapter 10 of Yarman
#    p=1j*w
#    f0=100e6
#    w0=2*pi*f0
#    R0=50
#    L1=0.225e-6*w0/R0
#    Y1=1.0/(L1*p)+1.0
#    L2=0.4e-6*w0/R0
#    Z2=(L2*p)+1/Y1
#    Y2=1/Z2
#    L3=0.01e-6*w0/R0
#    C3=11.4e-12*w0*R0
#    Z3a=L3*p+1/(C3*p)
#    Y3=1/Z3a+Y2
#    Z3=1/Y3
#    L4=(0.1e-6*w0/R0)
#    C4=(19.8e-12)*w0*R0
#    Z4a=L4*p+1/(C4*p)
#    Z4=Z4a+Z3
#    ZL=Z4
#    SL=(ZL-1)/(ZL+1)
#    return(SL)

def sopt(x,T0,ntr,k,nopt,Wa,reflectance_function,passband):
    '''Inputs:
        T0: Desired gain level
        ntr: Control flag for equalizer design (with or without transformer)
        x: vector includes the unknown coefficients of the polynomial h(p)
        for ntr=0,k=0 len(x) is n
        for ntr=1, len(x)=n+1'''
    
    wopt=linspace(passband[0],passband[1],nopt)
    L11=reflectance_function(wopt)
    fun=zeros_like(wopt)
    if ntr==1:
        h=poly1d(x)
    else:
        h=poly1d(append(x,0))
    g=Hurwitzpoly_g(h,k)
    T=gain(wopt,h,g,k,L11)
    fun=real(T-T0)
    return(fun)

def Hurwitzpoly_g(h,k,verbose=False):
    '''Inputs: 
        k is the order of f(p)=p**k
        h is a poly1d object 
    '''
    #step 1: generate the even polynomial H(-p**2)=h(p)h(-p)    
    hminus=poly1d([(-1)**(i)*coeff for i,coeff in enumerate(h.coeffs[::-1])][::-1])
    H=h*hminus
    if verbose==True:
        print "H=",H
    #Step 2: Generate G(-p**2)=H(-p**2)+(-1)**k    
    vec=zeros(len(H.coeffs))
    vec[-(2*k+1)]=(-1)**k    
    G=H+poly1d(vec)
    if verbose==True:
        print "G=",G
   
    #Generate the vector X and compute the root of G(-p**2)
    G2=poly1d((G.coeffs[::-2])[::-1])
    X=poly1d([(-1)**(i)*coeff for i,coeff in enumerate(G2.coeffs[::-1])][::-1])
    #Compute the RHP and LHP roots of G(-p**2)    
    Xr=X.roots
    pr=sqrt(-Xr)
    if verbose==True:
        print "pr=",pr
   
    prm=-pr
    #Step 4:Generate the monic polynomial g'(p)
    C=poly1d(prm,True)
    #Step 5: Generate Scattering Hurwitz polynomial g(p)
    Cof=sqrt(abs(G.coeffs[0]))
    if verbose==True:
        print "Cof=",Cof
   
    g=poly1d(C*Cof)
    return(g)
    
def gain(W,h,g,k,L11):
    #L11 is the load scattering matrix element L11.
    #L11 = (Zl-1)/(Zl+1)
    s=1j*W
    gval=polyval(g,s)
    hval=polyval(h,s)
    fval=s**k
    Weight=fval*conj(fval)*(1.0-L11*conj(L11))
    #denominator of the gain function
    D=hval*conj(hval)*(1.0+L11*conj(L11))+fval*conj(fval)-2*real(L11*hval*conj(gval))
    T=Weight/D
    assert (imag(T)==0).all()
    
    return(T)

from scipy import poly1d,zeros

def cauer_decomp(n,d):
    '''Input is the numerator and denominator polynomials in an expression
    for z.  Format for polynomials is Python-style with highest order coefficient
    first, lowest order last.
    
    Returns: 
        circuittype: designating whether it starts with a series inductor (+1)
        or shunt capacitor (-1)
        a: The coefficients of the decomposition and also the values of the 
        circuit components        
        '''    
    ct=n.order-d.order
    #if circuittype==1 then the circuit starts with a series inductor.
    #if circuittype==2 then the circuit starts with a shunt capacitor
    Ls,Cs=[],[]    
    if ct==1:
        #Extract pole at infinity from z
        polyorder=n.order
        for i in range(n.order):      
            div=poly1d(n)/poly1d(d)
            print div
            if mod(i,2)==0:            
                Ls.append(div[0].c[0])
                print "Adding L ",div[0].c[0]
            else:
                Cs.append(div[0].c[0])
                print "Adding C ",div[0].c[0]
            n=d            
            d=div[1]
            if i==(polyorder-2):
                if len(Ls)!=len(Cs):
                    print "Circuit type LC"
                    circuittype='LC'
                    R=d.coeffs[0]/n.coeffs[1]
                else:
                    print "Circuit type LCL"
                    circuittype='LCL'
                    R=n.coeffs[1]/d.coeffs[0]
        
    elif ct==-1:
        #Extract pole at infinity from y
        polyorder=d.order
        for i in range(d.order):      
            div=poly1d(d)/poly1d(n)
            print div
            if mod(i,2)==0:            
                Cs.append(div[0].c[0])
                print "Adding C ",div[0].c[0]
            else:
                Ls.append(div[0].c[0])
                print "Adding L ",div[0].c[0]
            d=n            
            n=div[1]
            if i==(polyorder-2):
                if len(Ls)!=len(Cs):
                    print "Circuit type CL"
                    circuittype='CL'
                    R=d.coeffs[1]/n.coeffs[0]
                else:
                    print "Circuit type CLC"
                    circuittype='CLC'
                    R=n.coeffs[0]/d.coeffs[1]
      
    else:
        raise Exception("Order of numerator and denominator must differ by 1")            
#    assert (array(Ls)>0).all(), "Inductor values are not all positive"
        if not (array(Ls)>0).all():
            print "Inductor values are not all positive"

#    assert (array(Cs)>0).all(), "Capacitance values are not all positive"
        if not (array(Cs)>0).all():
            print "Capacitance values are not all positive"
    return(circuittype,Ls,Cs,R)
        
def circuit_tikz_code(circuittype,Ls,Cs):
    code=r'''
\documentclass{article}
\usepackage[active,tightpage]{preview}
\usepackage[siunitx]{circuitikz}
\usepackage{tikz}
\usepackage{verbatim}
\begin{document}
\begin{preview}
\[
\begin{circuitikz} \draw
(0,2) node[anchor=east]{$In_1$}to[short, o-*] (1,2)
(0,0) node[anchor=east]{$In_2$}to[short, o-*] (1,0)
    '''
    #Determines if the first element in the circuit is a shunt cap or 
    #series inductor
    first_c_node={'CL':1,'LC':3}
    for i,L in enumerate(Ls):
        #Ground stage
        code+=r'(%d,0) -- (%d,0)' % (2*i+1,2*i+3)
        code+='\n'
        #Series inductor 
        if L>=1:
            unit,Lnew=r'\henry',L*1
        elif 1e-3 <= L <1:
            unit,Lnew=r'\milli\henry',L*1e3                 
        elif 1e-6 <= L <1e-3:
            unit,Lnew=r'\micro\henry',L*1e6
        elif 1e-9 <= L <1e-6:
            unit,Lnew=r'\nano\henry',L*1e9
        code+=r'(%d,2) to [L, l=$%0.1f%s$] (%d,2)' % (2*i+1,Lnew,unit,2*i+3) 
        code+='\n'
    for i,C in enumerate(Cs):
        #Shunt capacitor inductor 
        if C>=1:
             unit,Cnew=r'\farad',C*1                 
        elif 1e-3 <= C <1:
             unit,Cnew=r'\milli\farad',C*1e3
        elif 1e-6 <= C <1e-3:
             unit,Cnew=r'\micro\farad',C*1e6
        elif 1e-9 < C <1e-6:
             unit,Cnew=r'\nano\farad',C*1e9
        currentloc=2*i+first_c_node[circuittype]
        code+=r'(%d,0) to [C, l=$%0.1f%s$] (%d,2)' % (currentloc,Cnew,unit,currentloc)
        code+='\n'

        lastnodex=len(Ls)*2+1
    
    code+='''    
(%d,2) node[anchor=west]{$Out_1$}to[short, o-*] (%d,2)
(%d,0) node[anchor=west]{$Out_2$}to[short, o-*] (%d,0)
(%d,0) to[generic,l=$Z_l$] (%d,2)
;
\end{circuitikz}
\]
\end{preview}
\end{document}
    ''' % (lastnodex+1,lastnodex,lastnodex+1,lastnodex,lastnodex+1,lastnodex+1)
    
    f=open("tmp.tex",'w')
    f.write(code)
    f.close()
    import subprocess
    p = subprocess.Popen('pdflatex tmp.tex', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        #print line,
        pass
    p = subprocess.Popen('acroread tmp.pdf', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return(code)
    
def Total_Input_Impedance(omega,circuittype,Ls,Cs,Zl):
    '''Calculates the load impedance of the matching circuit plus transducer'''
    #It's confusing because cauer_decomp returns the network looking in from the
    #load side and we calculate the impedance looking in from the source side.
    Ztot=zeros([len(omega),2,2])    
    if circuittype=='LC':
        print "Matching network starts with shunt capacitor"        
        abcd=CL(omega,Ls[::-1],Cs[::-1])
    elif circuittype=='LCL':
        print "Matching network starts with series inductor"                        
        abcd=LC(omega,Ls[::-1],Cs[::-1])
    elif circuittype=='CL':
        print "Network starts with series inductor"
        abcd=LC(omega,Ls[::-1],Cs[::-1])
    elif circuittype=='CLC':
        print "Network starts with shunt capacitor"
        abcd=CL(omega,Ls[::-1],Cs[::-1])
    for i in range(len(omega)): 
        abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,0.0],[1.0/Zl[i],1.0]])
    Ztot=abcd[:,0,0]/abcd[:,1,0]
    return(Ztot)

def LC(omega,Ls,Cs):
    '''Returns the impedance of an series L, shunt C network that starts with a
    series L
    parameters: 
        Ls - series inductances starting at the source
        Cs - shunt capacitances starting at the source
    returns:
        An ABCD matrix with the impedance of the network as seen from the source'''
    abcd=zeros([len(omega),2,2],dtype='complex64')
    for i,om in enumerate(omega):    
        abcd[i,:,:]=array([[1,0],[0,1]])
        if len(Ls)>len(Cs):
            #one more L than C, so network starts and ends with an L
            for L,C in zip(Ls[:-1],Cs):       
                abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,1j*L*om],[0.0,1.0]])*matrix([[1.0,0.0],[1j*C*om,1.0]])
            abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,1j*Ls[-1]*om],[0.0,1.0]])
        elif len(Ls)==len(Cs):
            for L,C in zip(Ls,Cs):       
                abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,1j*L*om],[0.0,1.0]])*matrix([[1.0,0.0],[1j*C*om,1.0]])
    return(abcd)

def CL(omega,Ls,Cs):    
    '''Returns the impedance of an series L, shunt C network that starts with a
    shunt C
    parameters: 
        Ls - series inductances starting at the source
        Cs - shunt capacitances starting at the source
    returns:
        An ABCD matrix with the impedance of the network as seen from the source
    '''
    abcd=zeros([len(omega),2,2],dtype='complex64')
    for i,om in enumerate(omega):    
        abcd[i,:,:]=array([[1.0,0.0],[0.0,1.0]])
        if len(Ls)>len(Cs):
            #one more C than L, so network starts and ends with an C
            for L,C in zip(Ls[:],Cs[:-1]):        
                abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,0.0],[1j*C*om,1.0]])*matrix([[1.0,1j*L*om],[0.0,1.0]])
            abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,0.0],[1j*C*om,1.0]])
        elif len(Ls)==len(Cs):
            for L,C in zip(Ls,Cs):        
                abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,0.0],[1j*C*om,1.0]])*matrix([[1.0,1j*L*om],[0.0,1.0]])
    return(abcd)
    
def Z22(omega,circuittype,Ls,Cs,R):
    '''Calculates the source impedance of the matching circuit plus source resistance'''
    #It's confusing because cauer_decomp returns the network looking in from the
    #load side and we calculate the impedance looking in from the source side.
    Ztot=zeros([len(omega),2,2])    
    if circuittype=='LC':
        print "Matching network starts with shunt capacitor"        
        abcd=LC(omega,Ls,Cs)
    elif circuittype=='LCL':
        print "Matching network starts with series inductor"                        
        abcd=LC(omega,Ls,Cs)
    elif circuittype=='CL':
        print "Network starts with series inductor"
        abcd=CL(omega,Ls,Cs)
    elif circuittype=='CLC':
        print "Network starts with shunt capacitor"
        abcd=CL(omega,Ls,Cs)
    for i in range(len(omega)): 
        abcd[i,:,:]=matrix(abcd[i,:,:])*matrix([[1.0,0.0],[1.0/R,1.0]])
    Ztot=abcd[:,0,0]/abcd[:,1,0]
    return(Ztot)

def ReportPhysicalCircuitParams(nomcenter,Z0,Ls, Cs, R):
    Ls_physical=array(Ls)/nomcenter
    Cs_physical=array(Cs)/nomcenter
    R_physical=R*Z0
    
    print "-"*50
    print "Physical circuit paramaters in matching network"
    print "Ls: ",Ls_physical
    print "Cs: ",Cs_physical
    print "R: ",R_physical
    print "-"*50
    
    return(Ls_physical,Cs_physical,R_physical)
    
    

if __name__=='__main__':
    ##Sanity Checks
    #htest=poly1d([1,2,1,0])
    #gtest=Hurwitzpoly_g(htest,0,verbose=True)
    #print "gtest=" , gtest.coeffs
    #gtest should be g(p) = p^3 + 2.6494p^2 + 2.5098p + 1.
    #
    #htest=poly1d([5,3,1,0.5,1,0])
    #gtest=Hurwitzpoly_g(htest,0,verbose=True)
    #print "gtest=" , gtest.coeffs
    #gtest should be g(p) = [5,12.15,14.86,10.71,4.73,1]
    #
    ##k=3
    #htest=poly1d([1,-2,1,-0,1,-1,-1.1,3,1,-1,1])
    #gtest=Hurwitzpoly_g(htest,3,verbose=True)
    #print "gtest=" , gtest.coeffs
    ##gtest should be g(p) = #gtest should be g(p) = [5,12.15,14.86,10.71,4.73,1]
    ##Define the function to be optimized
    
    #***************************
    f,zr,zi=loadtxt("./Sample data/green_probe.txt",unpack=True)
    Zorig=(zr+1j*zi)/50
    w=2*pi*f
    nomcenter=45e6*2*pi
    #Normalize the frequency
    w=w/nomcenter
    from scipy.interpolate import interp1d
    zfunc=interp1d(w,Zorig,kind='cubic')
    R0=1.0
    nopt=100
    ntr=0
    k=0
    T0=1.0
    passband=[0.9,1.5]
    
    #****************************
    
    
    #***************************
    #Zorig,R0=short_monopole_antenna_Z(w),50.0
    #zfunc=short_monopole_antenna
    #****************************
    
    #***************************
    #Example from chapter 9, page 220
#    w=linspace(0,1.6,1000)
#    Zorig,R0=example4_Z(w),1.0
#    zfunc=example4_Z
#    nopt=20
#    ntr=1
#    k=0
#    T0=0.8
#    passband=[0.2,1.0]
#    hin=poly1d([-1,1,1,1,-1])
    
    #h=[ -1.17844181e-03,-4.89334921e-01,-3.83191522e-01,-4.28596764e+00,-1.29137330e+00,-4.34124124e+00,4.04257482e-02]
    #Gives a negative T!
    #****************************
    
    
    close('all')    
    print "Using R0=",R0
    reflection=lambda w: reflection_from_Z(w,zfunc,R0)
    L11=reflection(w)
    #Setup 
    #Calculate T for no matching network
    T=gain(w,poly1d(1),poly1d(1),k,L11)
    
    order=5
        
    try_all_bit_patterns=True
    if try_all_bit_patterns:
        from bitstring import BitArray
        figure(13)
        plot(f/1e6,gain(w,poly1d(1),poly1d(1),k,L11),label="No match")
        xlabel("Frequency (MHz)")
        ylabel("Gain")
        sols=[]
        sos=zeros(2**order)
    
        for i in range(0,2**order):
            print "-"*50
            ba=sign(array(list(BitArray(uint=i,length=order)),dtype='float64')-0.5)
            hin=poly1d(sign(array(list(ba),dtype='float64')-0.5))
            print hin
            #Least squares solution
            print "Performing least-squares optimization..."
            sol,cov_x,infodict,mesg,ier=leastsq(sopt,hin.coeffs,maxfev=10000,args=(T0,ntr,k,nopt,w,reflection,passband),full_output=True)
            #print mesg
            #print "Error Code: %d" % ier
            #print "Solution is:", sol
            #print "Least squares sum=",sum(sopt(sol,T0,ntr,k,nopt,w,reflection,passband)**2)
            hbook=poly1d([-1.8884, -1.9863, -1.9372, -0.6279, -0.0468])
            gbook=Hurwitzpoly_g(hbook,k)
            hfit=poly1d(sol)
            gfit=Hurwitzpoly_g(hfit,k)
            sos[i]=sum(sopt(sol,T0,ntr,k,nopt,w,reflection,passband)**2)
            print "Least squares sum from book=",sum(sopt(hbook.coeffs,T0,ntr,k,nopt,w,reflection,passband)**2)
            circuittype,Ls,Cs,R=cauer_decomp((gfit-hfit),(gfit+hfit))
            if (array(Ls)>0).all() and (array(Cs)>0).all() and R>0:
                print "Valid solution found!"
                print "Ls: ", Ls
                print "Cs: ", Cs
                print "R: ",R
                figure(13)
                #plot(w,gain(w,hbook,gbook,k,L11),label="Book solution")
                plot(f/1e6,gain(w,hfit,gfit,k,L11),label="%d" %i)
                sols.append(sol)    
            
        sol=sols[argmin(sos)]
        print "Best solution achieves sum of squares of ", sols[argmin(sos)]
        grid()
        legend()
        hlines(T0,amin(w),amax(w))
        figure(83)
        scatter(range(0,2**order),sos)
        hlines(amin(sos),0,2**order)
    
    else:
        hin=poly1d([-1.94335442,1.11208621,-8.01659691,8.15526687,-9.11459217,16.75179636,-2.18424997,9.86945929])
        #hin=poly1d(rand(order))
        sol,cov_x,infodict,mesg,ier=leastsq(sopt,hin.coeffs,maxfev=10000,args=(T0,ntr,k,nopt,w,reflection,passband),full_output=True)
    
    
    hfit=poly1d(sol)
    #hbook=array([-1.0607, -2.3496, -2.7543, -1.8602, -0.5515])
    gfit=Hurwitzpoly_g(hfit,k)
    
    figure(11)
    title("Best solution")
    plot(f/1e6,gain(w,poly1d(1),poly1d(1),k,L11))
    plot(f/1e6,gain(w,hfit,gfit,k,L11),label="Best fit")
    xlabel("Frequency (MHz)")
    ylabel("Gain")
    legend()
    grid()
    print "Roots of gfit-hfit"
    print (gfit-hfit).roots    
    print "Roots of gfit-hfit"
    print (gfit+hfit).roots    
     
    print "Performing Cauer decomposition..."
    circuittype,Ls,Cs,R=cauer_decomp((gfit-hfit),(gfit+hfit))
    print "Ls=",Ls
    print "Cs=",Cs
    print "R=",R
    
    print "Calculating matching LC network..."
    
    figure(131)
    s=1j*w
    scatter(w,abs((gfit(s)-hfit(s))/(gfit(s)+hfit(s))),label="Desired matching network impedance",color='red')
    
    abcd=zeros([len(w),2,2],dtype='complex64')
    ztest=zeros(len(w),dtype='complex64')
    z2=zeros(len(w),dtype='complex64')  
    for i in range(len(w)):    
        if circuittype=='LC':
            abcd[i,:,:]=matrix([[1.0,Ls[0]*s[i]],[0,1.0]])*matrix([[1.0,0],[Cs[0]*s[i],1.0]])*matrix([[1.0,Ls[1]*s[i]],[0,1.0]])*matrix([[1.0,0],[Cs[1]*s[i],1.0]])*matrix([[1.0,0.0],[1.0/R,1.0]])                
        if circuittype=='CL':
            abcd[i,:,:]=matrix([[1.0,0],[Cs[0]*s[i],1.0]])*matrix([[1.0,Ls[0]*s[i]],[0.0,1.0]])*matrix([[1.0,0.0],[Cs[1]*s[i],1.0]])*matrix([[1.0,Ls[1]*s[i]],[0.0,1.0]])*matrix([[1.0,0.0],[1.0/R,1.0]])
        ztest[i]=abcd[i,0,0]/abcd[i,1,0]
    E22net=(1.0-ztest)/(1.0+ztest)    
    plot(f/1e6,abs(ztest),label="Manual ABCD")
  
    
    
    ztot=Total_Input_Impedance(w,circuittype,Ls,Cs,Zorig)
   
    print "Plotting impedance"
    figure(2)
    title("Impedances")
    subplot(211)
    title("Re{Z}")
    plot(f/1e6,real(Zorig)*50,label="Unmatched")
    plot(f/1e6,real(ztot)*50,label="Matched")
    xlabel("Frequency (MHz)")
    ylabel(r"Z ($\Omega$)")
    grid('on',which='both')
    legend()
    subplot(212)
    title("Im{Z}")
    plot(f/1e6,imag(Zorig)*50,label="Unmatched")
    plot(f/1e6,imag(ztot)*50,label="Matched")
    xlabel("Frequency (MHz)")
    ylabel(r"Z ($\Omega$)")
    grid('on',which='both')
    legend()
    
    
    
    zs=Z22(w,circuittype,Ls,Cs,R)
    E22=(1-zs)/(1+zs)    
    Tman3=4*real(zs)*real(Zorig)/abs(zs+Zorig)**2  
    
    figure(20)
    Tmatch=((1-abs(E22)**2)*(1-abs(L11)**2))/(abs(1-E22*L11)**2)
    title("Transmission")
    plot(f/1e6,T,label="Unmatched")
    plot(f/1e6,Tmatch,label="Matched",lw=4)
    
    plot(f/1e6,Tman3,label="From Impedance 3",lw=1)    
    grid('on',which='both')
    ylabel("Transmission")
    xlabel("Frequency (MHz)")
    
    legend()
    show()
    #Order of Cs and Ls must be reversed to display the circuit looking in from 
    #a source on the left
    #print "Generating Tikz code..."
    #circuit_tikz_code(ct,Ls[::-1],Cs[::-1])
    
    #TODO - check and make sure that gain formula works.  Plug in the fit value os
    #E22.  Try to figure out that big resonance in the total impedance..
    
    ReportPhysicalCircuitParams(nomcenter,50.0,Ls,Cs,R)
    
    print "Optimization complete."
    
    #Four circuit types CLC, CL,LC, LCL
