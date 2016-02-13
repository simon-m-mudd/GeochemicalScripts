## Weathering_zone_thickness.py
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## Makes plots of weathering zone thickness as a function of various
## parameters following Lebedeva et al 2010 ESPL
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
## SMM 13/02/2016
##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

import matplotlib.pyplot as pp
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import rcParams
import matplotlib.lines as mpllines
from matplotlib.ticker import FormatStrFormatter

# This is equation 24 in Lebedeva et al 2010
# ums is velocity squared
def calculate_beta(k, D, phi, ums):
  beta = np.sqrt(1+4*k*D*phi/(ums*ums))
  return beta
  
# This is equation 4 in Lebedeva et al 2010  
def calculate_k(k_ab, phi_ab, s_ab, Psi_ab):
  k = 2*k_ab*phi_ab*s_ab*Psi_ab
  
  print k
  return k

# This calculates the specific surface area for a given diameter
def calculate_geom_ssa(pDiam):
  ssa = 6/pDiam
  return ssa
 
# calculates the velocity in metres per second given velocity in metres per year 
def calculate_ums(umyr):
  ums = umyr/(365*24*60*60)
  return ums

# This is equation 23 in Lebedeva et al 2010 
def calculate_weath_thick(u,ssa, k_ab, phi, phi_ab, Psi_ab, D,Ce,Cr,Cl):


  k = calculate_k(k_ab, phi_ab, ssa, Psi_ab)
  beta = calculate_beta(k,D,phi,u)    
  leading_term = 2*D*phi/(u*(beta-1))
  log_term = (2*beta*(Ce-Cr))/((beta+1)*(Ce-Cl))
  
  delta = leading_term*np.log(log_term)
  return delta


# Plots the sensitivity of the weathering thickness to flow rate
def Plot_flow_rate_sensitivity():

  # These set the font size
  label_size = 20
  #title_size = 30
  axis_size = 28

  # set some variables for the fonts
  #rcParams['font.family'] = 'sans-serif'
  #rcParams['font.sans-serif'] = ['arial']
  
  rcParams['font.size'] = label_size  

  logu = np.arange(-1,1,0.01)
  print logu
  
  u = np.power(10,logu)
  #print u

  Dp = [0.000064, 0.000128,0.000256]
  Dpart = np.asarray(Dp)
  #print Dpart
  ssa = calculate_geom_ssa(Dpart)
  print "ssa"  
  print ssa
  
  
  k_ab = 3.87*np.power(10.0,-10)
  phi = 0.2
  phi_ab = 0.4
  Psi_ab = 1.2*np.power(10.0,-2)
  D = 0.8*np.power(10.0,-9)
  #print "Yo, checking"
  #print D
  umyr = 0.2
  this_Dp = 0.001;  
  
  Ce = 0.2
  Cr = 0
  Cl = 0.9*Ce
  
  #ssa = calculate_geom_ssa(Dpart)
  ums = calculate_ums(u)
  
  #ssa = 3.5*np.power(10.0,4)
  
    
  
  delta1 = calculate_weath_thick(ums,ssa[0], k_ab, phi, phi_ab, Psi_ab, D,Ce,Cr,Cl)
  delta2 = calculate_weath_thick(ums,ssa[1], k_ab, phi, phi_ab, Psi_ab, D,Ce,Cr,Cl)
  delta3 = calculate_weath_thick(ums,ssa[2], k_ab, phi, phi_ab, Psi_ab, D,Ce,Cr,Cl)
  #print delta1
  #print delta2
  #print delta3
  
  
    
  # now plot the results    
  fig = pp.figure(1, facecolor='white',figsize=(10,7.5))
  ax1 = fig.add_subplot(1,1,1)  
  
  pp.plot(u,delta1,'k',u,delta2, u,delta3,linewidth=3)

  
  pp.rcParams['xtick.direction'] = 'out'
  pp.rcParams['ytick.direction'] = 'out'
  ax1.set_xscale('log')
  #ax1.set_yscale('log')
  
  
  ax1.spines['top'].set_linewidth(2.5)
  ax1.spines['left'].set_linewidth(2.5)
  ax1.spines['right'].set_linewidth(2.5)
  ax1.spines['bottom'].set_linewidth(2.5) 
  ax1.tick_params(axis='both', width=2.5)    
  ax1.xaxis.set_major_formatter(FormatStrFormatter('%.01f'))

  for line in ax1.get_xticklines():
    line.set_marker(mpllines.TICKDOWN)

  for line in ax1.get_yticklines():
    line.set_marker(mpllines.TICKLEFT)

  pp.ylabel('Weathering zone thickness ($m$)',fontsize = axis_size)
  pp.xlabel('Flow velocity, $m$ $yr^{-1}$',fontsize = axis_size) 

  pp.tight_layout()
      
  #pp.savefig("Particle_size_thick.svg", format='svg')  
  pp.show()
  
  #pp.cla()
  #pp.clf()
  #pp.close(fig)

  
if __name__ == "__main__":
    Plot_flow_rate_sensitivity()  