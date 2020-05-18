#requires numpy
#import numpy as np
#functions used for plotting of appendix fit figures in https://ui.adsabs.harvard.edu/abs/2020arXiv200107785L/abstract

def nbr_avg_fit_fesc_z10(Ms): 
    
    """
    For a given mass in solar masses Ms, return the fit of the average CoDa ii
    escape fraction at z=10.1
    """
    
    assert np.any([Ms<1e12,Ms>1e8]), "Found masses outside of [1e8,1e12] solar mass range!"
    
    slope=(np.log10(0.77)-np.log10(0.1))/(np.log10(1)-np.log10(1000))
    return(0.77*(Ms/1e8)**(slope))
    
    
    
def nbr_avg_fit_fesc(Ms,Zs):                                           

    
    """
    For a given mass in solar masses Ms, return the fit of the number weighted average CoDa ii
    escape fraction at z=Zs. With 6<=z<=14.9
    """
    
    fesc_max0,fesc_knee0,beta,zeta,gamma=[0.98,  0.8, -0.5, -2,  -3.5]
    
    assert Zs>=6, "Redshift too low! Minimum of z=6."
    assert Zs<=14.9, "Redshift too high! Maximum of z=14.9"
    
    assert np.any([Ms<1e12,Ms>1e8]), "Found masses outside of [1e8,1e12] solar mass range!"
             
    z10=nbr_avg_fit_fesc_z10(Ms)
    
    if Zs<10.1:
    
        Mknee=2.5e9/(((1.+Zs)/(1.+6))**zeta)                                                          
        fesc_max=fesc_max0*((1.+Zs)/(1.+6))**beta                                                          
        fesc_knee=fesc_knee0*((1.+Zs)/(1.+6))**gamma                                                          

        fesc_knee=np.max([fesc_knee,0.23])    
        delta_knee=(-np.log10(fesc_max)+np.log10(fesc_knee))/(-np.log10(1e8)+np.log10(Mknee))
        delta_max=(-np.log10(0.1)+np.log10(fesc_knee))/(-np.log10(1e11)+np.log10(Mknee))
    
        down=fesc_knee*(Ms/Mknee)**(delta_knee)                                                                                      
        up=fesc_knee*(Ms/Mknee)**(delta_max)                                                                
        
        out=np.where(Ms>=Mknee,up,down)        
        return(np.where(Ms<=1e11,out,z10))  
    
    else :
        
        return(z10**(10.1/Zs)**-1.64)
        
def SFR_avg_fit_fesc(Ms,Zs):
    
    """
    For a given mass in solar masses Ms, return the fit of the SFR weighted average CoDa ii
    escape fraction at z=Zs. With 6<=z<=14.9
    """
    
    assert Zs>=6, "Redshift too low! Minimum of z=6."
    assert Zs<=14.9, "Redshift too high! Maximum of z=14.9"
    
    assert np.any([Ms<1e12,Ms>1e8]), "Found masses outside of [1e8,1e12] solar mass range!"
    
    fesc_max0,fesc_knee0,delta,beta,gamma=[1.0,0.65,-0.75,-0.5,-2.5] 
    
    Mknee=3e9/(((1.+Zs)/(1.+6)))                                                           
    fesc_max=fesc_max0*((1.+Zs)/(1.+6))**beta                                                           
    fesc_knee=fesc_knee0*((1.+Zs)/(1.+6))**gamma                                                           
                                                                                   
    fesc_max=fesc_knee*(Ms/Mknee)**((-np.log10(fesc_max)+np.log10(fesc_knee))/(-np.log10(1e8/Mknee)+\
np.log10(1)))                                                                                        
    fesc_knee=fesc_knee*(Ms/Mknee)**delta                                                                 
                                                                                                     
    return(np.where(Ms<Mknee,fesc_max,fesc_knee))  
