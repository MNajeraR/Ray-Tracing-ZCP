#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 13:45:51 2021

@author: Morgan
"""

import numpy as np
import matplotlib.pyplot as plt

###########################################################################
def surface_Derivative(x,y,z,i):
    delta=0.000001
    
    Dx=(surface(x+delta,y,z,i)-surface(x-delta,y,z,i))/(2.0*delta)
    Dy=(surface(x,y+delta,z,i)-surface(x,y-delta,z,i))/(2.0*delta)
    Dz=(surface(x,y,z+delta,i)-surface(x,y,z-delta,i))/(2.0*delta)   
    
    Dr=np.sqrt((Dx*Dx)+(Dy*Dy)+(Dz*Dz))
    return Dx/Dr,Dy/Dr,Dz/Dr         
    
    
###########################################################################
def Snell_refraction_vector(Ray_Vect,Surf_Normal,n_0,n_1):
    
        
    Ray_Vect=np.asarray(Ray_Vect)
    Surf_Normal=np.asarray(Surf_Normal)       
    
    Nsurf_Cros_s1=np.cross(Surf_Normal,Ray_Vect)
        
    NN=n_0/n_1
    
    R=(NN*NN)*np.dot(Nsurf_Cros_s1,Nsurf_Cros_s1)
            
    if R>1:
        print("---- reflexion interna= ", R)
        
    s2=NN*(np.cross(Surf_Normal,np.cross(-Surf_Normal,Ray_Vect)))-Surf_Normal*np.sqrt(1.0-((NN*NN)*np.dot(Nsurf_Cros_s1,Nsurf_Cros_s1)))
   
    return s2

###########################################################################

def surface(x_1,y_1,z_1,i):
    
    R_c=Rc[i]
    k=kk[i]
    Alpha_1=aalpha_1[i]
    Alpha_2=aalpha_2[i]
    Alpha_3=aalpha_3[i]
    Alpha_4=aalpha_4[i]
    Alpha_5=aalpha_5[i]

  
    c=1.0/R_c
    s2=(x_1*x_1)+(y_1*y_1)
    s4=s2*s2
    s6=s4*s2
    s8=s6*s2
    s10=s8*s2
    
    f1=c*s2/(1.0+np.sqrt(1.0-((k+1.0)*c*c*s2)))
    f2=(Alpha_1*s2)+(Alpha_2*s4)+(Alpha_3*s6)+(Alpha_4*s8)+(Alpha_5*s10)
    return f1+f2-z_1
###########################################################################

def F_RS(z_1,rayorigin,i):
    
    
    x_0=rayorigin[0]
    y_0=rayorigin[1]
    z_0=rayorigin[2]
    L=rayorigin[3]
    M=rayorigin[4]
    N=rayorigin[5]
    
    
    x_1=(L/N)*(z_1-z_0)+x_0
    y_1=(M/N)*(z_1-z_0)+y_0
    
    f_RS=surface(x_1,y_1,z_1,i)    
    
    return f_RS

################################################################################

def cosdir(x,y,z):
  magni = np.sqrt(x**2+y**2+z**2)
  L = x/magni
  M = y/magni
  N = z/magni

  return L, M, N

###########################################################################
def Der_F_RS(z_1,rayorigin,i):
    h=0.00000000001
    der=(F_RS(z_1+h, rayorigin,i)-F_RS(z_1-h, rayorigin,i))/(2*h)
    return der
###########################################################################

def Newthon_Raphson(rayorigin,i):
    z_1=0
    cnt=0
    while True:
        z_2=z_1-(F_RS(z_1,rayorigin,i)/Der_F_RS(z_1,rayorigin,i))
        
        if np.abs(z_2-z_1)<0.00001:
            
            break
        else:
            z_1=z_2
        if cnt == 30:
            break    
        cnt=cnt+1    

        
    [x_0, y_0, z_0,L,M,N]=rayorigin
    x_1=(L/N)*(z_1-z_0)+x_0
    y_1=(M/N)*(z_1-z_0)+y_0
        
    return x_1, y_1, z_1

###########################################################################


def Trax(Ray_Param,i):
    ## punto de intercección en la superficie
    x_1, y_1, z_1 = Newthon_Raphson(Ray_Param,i)
    
    ## vector normal en el punto de intercección en la superficie
    NORM=surface_Derivative(x_1,y_1,z_1,i)
    [x_0, y_0, z_0,L,M,N]=Ray_Param
    Ray_Vect=[L,M,N]
    Surf_Normal=NORM
    
    
    
    n_0=nn0[i]
    n_1=nn1[i]
        
    lmn=Snell_refraction_vector(Ray_Vect,Surf_Normal,n_0,n_1)
    
    return([x_1, y_1, z_1, lmn[0], lmn[1], lmn[2]])

################################################################################

def trasform(TRF):
    [TiltX, TiltY, TiltZ, dx, dy, dz]=TRF
    Tx=np.deg2rad(TiltX)
    Ty=np.deg2rad(TiltY)
    Tz=np.deg2rad(TiltZ)    
    
        
    # https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    Rx_1A = np.matrix([[1.0, 0.0, 0.0 , 0.0], [0.0, np.cos(-Tx), np.sin(-Tx),0.0],[0.0,-np.sin(-Tx), np.cos(-Tx),0.0],[0.0, 0.0, 0.0 ,1.0]])
    Ry_1A = np.matrix([[np.cos(-Ty),0.0, -np.sin(-Ty),0.0], [0, 1.0, 0.0, 0.0],[np.sin(-Ty),0.0, np.cos(-Ty),0.0],[0.0, 0.0, 0.0 ,1.0]])
    Rz_1A = np.matrix([[np.cos(-Tz),-np.sin(-Tz), 0.0, 0.0],[np.sin(-Tz),np.cos(-Tz), 0.0, 0.0],[0.0, 0.0, 1.0, 0.0],[0, 0, 0 ,1.0]])    
    Dxyz_1A = np.matrix([[1.0, 0.0, 0.0, -dx],[0.0, 1.0, 0.0, -dy],[0.0, 0.0, 1.0, -dz],[0.0, 0.0, 0.0 ,1.0]])
    
    Start_trans1=Rz_1A*Ry_1A*Rx_1A*Dxyz_1A
    Start_trans2=np.linalg.inv(Start_trans1)
    
        
    return Start_trans1, Start_trans2

###############################################################################

#def Rayos_int(x_0, y_0, z_0, L, M, N, R_c, k, Alpha_1, Alpha_2, Alpha_3, Alpha_4, Alpha_5, n0, n1):
def Rayos_int(x_0, y_0, z_0, L, M, N, i):    

    Ray_Param = [x_0, y_0, z_0, L, M, N]
    #Surf_Param=[R_c,k, Alpha_1, Alpha_2, Alpha_3, Alpha_4, Alpha_5, n0, n1]
    xyz_1=Trax(Ray_Param,i)
    
    return([xyz_1[0], xyz_1[1], xyz_1[2], xyz_1[3], xyz_1[4], xyz_1[5]])

###############################################################################

def traceray(inx,iny,inz,inL,inM,inN, MT):
    

    coorx_0[0] = inx
    coory_0[0] = iny
    coorz_0[0] = inz
    cosL[0] = inL
    cosM[0] = inM
    cosN[0] = inN
    
    for i in range (5):
        A1_Trans = MT[i][0]
        A2_Trans = MT[i][1]
        
        x_0 = coorx_0[i] + dxps[i]
        y_0 = coory_0[i] + dyps[i]
        z_0 = coorz_0[i] - d[i] 
        
        L = cosL[i]
        M = cosM[i]
        N = cosN[i]
                
        Ray_Param = [x_0, y_0, z_0, L, M, N]
        xyz_1=Trax(Ray_Param,i)
        
        # #################################################

        x_1 = xyz_1[0]
        y_1 = xyz_1[1]
        z_1 = xyz_1[2]
        

        coor_pm0 = np.dot(A1_Trans, ([x_0, y_0, z_0, 1.0]))
        coor_pm1 = np.dot(A1_Trans,([x_1, y_1, z_1, 1.0]))
               
        scosD_pm = cos_prim(coor_pm0, coor_pm1) 
         
        x_0 = coor_pm0[0,0]
        y_0 = coor_pm0[0,1]
        z_0 = coor_pm0[0,2]
        


        # print(" .--.-.-.-.-.-.-.-.-.-.-..")
        # print(i)
        # print(L, M, N)
        
        s=np.sign(d[i])
        
        if i==0:
            s=1
            
        L = scosD_pm[0]*s
        M = scosD_pm[1]*s
        N = scosD_pm[2]*s
        #print(L, M, N)
        # print(s)
        
        
        Ray_Param = [x_0, y_0, z_0, L, M, N]
        xyz_1=Trax(Ray_Param,i)
        
        x_1 = xyz_1[0]
        y_1 = xyz_1[1]
        z_1 = xyz_1[2]
     
        L = xyz_1[3]
        M = xyz_1[4]
        N = xyz_1[5]
        
        z_2=z_0
        
        x_2=(L/N)*(z_2-z_1)+x_1
        y_2=(M/N)*(z_2-z_1)+y_1
        
        coor_pm0 = np.dot(A2_Trans, ([x_1, y_1, z_1, 1.0]))
        coor_pm1 = np.dot(A2_Trans, ([x_2, y_2, z_2, 1.0]))
        scosD_pm = cos_prim(coor_pm0,coor_pm1) 
        
        x_1 = coor_pm0[0,0]
        y_1 = coor_pm0[0,1]
        z_1 = coor_pm0[0,2]
        
        L = scosD_pm[0]
        M = scosD_pm[1]
        N = scosD_pm[2]
        
        
        
        coorx_1[i] = x_1
        coory_1[i] = y_1
        coorz_1[i] = z_1
        cosL[i+1] = L
        cosM[i+1] = M
        cosN[i+1] = N
        # ######################################################
        
        
        # coorx_1[i] = xyz_1[0]
        # coory_1[i] = xyz_1[1]
        # coorz_1[i] = xyz_1[2]
        # cosL[i+1] = xyz_1[3]
        # cosM[i+1] = xyz_1[4]
        # cosN[i+1] = xyz_1[5]

        
        coorx_0[i+1] = coorx_1[i]
        coory_0[i+1] = coory_1[i]
        coorz_0[i+1] = coorz_1[i]

    return coorx_1,coory_1, coorz_1, cosL, cosM, cosN

##############################################################################

def cos_prim(sistema_0,sistema_1):

    x_0 = sistema_0[0,0]
    y_0 = sistema_0[0,1]
    z_0 = sistema_0[0,2]
    
    x_1 = sistema_1[0,0]
    y_1 = sistema_1[0,1]
    z_1 = sistema_1[0,2]

    magni = np.sqrt((x_1-x_0)**2+(y_1-y_0)**2+(z_1-z_0)**2)
    
    L_pm = (x_1-x_0)/magni
    M_pm = (y_1-y_0)/magni
    N_pm = (z_1-z_0)/magni
    
    return L_pm, M_pm, N_pm

###############################################################################

def Trazado_rayos(ray_param,Matrix):

  coorinc_0 = ray_param[0]
  coorinc_1 = ray_param[1]
  coorinc_2 = ray_param[2]
  cosedir_0 = ray_param[3]
  cosedir_1 = ray_param[4]
  cosedir_2 = ray_param[5]
  NR2 = len(coorinc_0)
  sistcoorx_5 = np.ones(NR2)
  sistcoory_5 = np.ones(NR2)
  sistcoorz_5 = np.ones(NR2)
  sistcosL_5 = np.ones(NR2)
  sistcosM_5 = np.ones(NR2)
  sistcosN_5 = np.ones(NR2) 
  sistcoorx_4 = np.ones(NR2)
  sistcoory_4 = np.ones(NR2)
  sistcoorz_4 = np.ones(NR2)
  sistcosL_4 = np.ones(NR2)
  sistcosM_4 = np.ones(NR2)
  sistcosN_4 = np.ones(NR2) 
  sistcoorx_3 = np.ones(NR2)
  sistcoory_3 = np.ones(NR2)
  sistcoorz_3 = np.ones(NR2)
  sistcosL_3 = np.ones(NR2)
  sistcosM_3 = np.ones(NR2)
  sistcosN_3 = np.ones(NR2) 
  sistcoorx_2 = np.ones(NR2)
  sistcoory_2 = np.ones(NR2)
  sistcoorz_2 = np.ones(NR2)
  sistcosL_2 = np.ones(NR2)
  sistcosM_2 = np.ones(NR2)
  sistcosN_2 = np.ones(NR2) 
  sistcoorx_1 = np.ones(NR2)
  sistcoory_1 = np.ones(NR2)
  sistcoorz_1 = np.ones(NR2)
  sistcosL_1 = np.ones(NR2)
  sistcosM_1 = np.ones(NR2)
  sistcosN_1 = np.ones(NR2)
  planoimagen = np.ones(NR2)

  for i in range(NR2):

      inx = coorinc_0[i]
      iny = coorinc_1[i]
      inz = coorinc_2[i]
      inL = cosedir_0[i]
      inM = cosedir_1[i]
      inN = cosedir_2[i]

      cx_1, cy_1, cz_1, cL, cM, cN = traceray(inx,iny,inz,inL, inM, inN, Matrix)

      sistcoorx_5[i] = cx_1[4]
      sistcoory_5[i] = cy_1[4]
      sistcoorz_5[i] = cz_1[4]
      sistcosL_5[i] = cL[5]
      sistcosM_5[i] = cM[5]
      sistcosN_5[i] = cN[5]
      sistcoorx_4[i] = cx_1[3]
      sistcoory_4[i] = cy_1[3]
      sistcoorz_4[i] = cz_1[3]
      sistcosL_4[i] = cL[4]
      sistcosM_4[i] = cM[4]
      sistcosN_4[i] = cN[4]
      sistcoorx_3[i] = cx_1[2]
      sistcoory_3[i] = cy_1[2]
      sistcoorz_3[i] = cz_1[2]
      sistcosL_3[i] = cL[3]
      sistcosM_3[i] = cM[3]
      sistcosN_3[i] = cN[3]
      sistcoorx_2[i] = cx_1[1]
      sistcoory_2[i] = cy_1[1]
      sistcoorz_2[i] = cz_1[1]
      sistcosL_2[i] = cL[2]
      sistcosM_2[i] = cM[2]
      sistcosN_2[i] = cN[2]
      sistcoorx_1[i] = cx_1[0]
      sistcoory_1[i] = cy_1[0]
      sistcoorz_1[i] = cz_1[0]
      sistcosL_1[i] = cL[1]
      sistcosM_1[i] = cM[1]
      sistcosN_1[i] = cN[1]

  # #########################################################################
    
  planoobjeto = (coorinc_0, coorinc_1, coorinc_2, cosedir_0, cosedir_1, cosedir_2)
  espejoprimario = (sistcoorx_1, sistcoory_1, sistcoorz_1, sistcosL_1, sistcosM_1, sistcosN_1)
  espejosecundario = (sistcoorx_2, sistcoory_2, sistcoorz_2, sistcosL_2, sistcosM_2, sistcosN_2)
  lentecorrectora_cara1 = (sistcoorx_3, sistcoory_3, sistcoorz_3, sistcosL_3, sistcosM_3, sistcosN_3)
  lentecorrectora = (sistcoorx_4, sistcoory_4, sistcoorz_4, sistcosL_4, sistcosM_4, sistcosN_4)
  planoimagen = (sistcoorx_5, sistcoory_5, sistcoorz_5, sistcosL_5, sistcosM_5, sistcosN_5)
  surf_opt = (espejoprimario, espejosecundario, lentecorrectora_cara1)
  plane_opt = (planoobjeto, lentecorrectora, planoimagen)

  return surf_opt, plane_opt


##############################################################################

def rotat_trasl(tip_1, tilt_1, set_rayos):
  tipx = tip_1[0]
  tipy = tip_1[1]
  tiltx = tilt_1[0]
  tilty = tilt_1[1]
  MT=[]
  ma_rt_s0 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s1 = [tipx , tipy, 0.0, tiltx, tilty, 0.0]
  ma_rt_s2 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s3 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s4 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s5 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s6 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s7 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  M = trasform(ma_rt_s0)
  MT.append(M)
  M = trasform(ma_rt_s1)
  MT.append(M)
  M = trasform(ma_rt_s2)
  MT.append(M)
  M = trasform(ma_rt_s3)
  MT.append(M)
  M = trasform(ma_rt_s4)
  MT.append(M)
  surf_rottras, plane_rottras = Trazado_rayos(set_rayos,MT)
  return surf_rottras, plane_rottras


###############################################################################

def Height_yp(Set_lsurface):
    last_surface = Set_lsurface
    valx_LS = last_surface[0]
    valy_LS = last_surface[1]
    valz_LS = last_surface[2]
    valL_LS = last_surface[3]
    valM_LS = last_surface[4]
    valN_LS = last_surface[5]


    yp0 = valy_LS[0]
    yma = valy_LS[1]
    ymb = valy_LS[2]

    zp0 = valz_LS[0]  
    zma = valz_LS[1]
    zmb = valz_LS[2]

    Mp0 = valM_LS[0] 
    Mma = valM_LS[1]
    Mmb = valM_LS[2]

    Np0 = valN_LS[0] 
    Nma = valN_LS[1]
    Nmb = valN_LS[2]

    ym = ((Mma*Mmb)/(Nma*Mmb-Nmb*Mma))*(zmb-zma-ymb*(Nmb/Mmb)+yma*(Nma/Mma))
    zm = (ym-ymb)*(Nmb/Mmb)+zmb

    yp = yp0+(Mp0/Np0)*(zm-zp0)

    dt_yp = np.abs(yp-ym)

    return dt_yp

###############################################################################

def Der_highyp(tilt_1, set_rayos):
    tilt_1 = np.array(tilt_1)
    tip_1 = (0.0, 0.0)
    h = 0.000000001
    tilt1 = tilt_1 + h
    tilt2 = tilt_1 - h
    surf_rph, plane_rph = rotat_trasl(tip_1,tilt1,set_rayos)
    surf_rnh, plane_rnh= rotat_trasl(tip_1,tilt2,set_rayos)
    l_surfaceph = plane_rph[2]
    l_surfacenh = plane_rnh[2]
    del_ymph = Height_yp(l_surfaceph)
    del_ymnh= Height_yp(l_surfacenh)
    der_delym = (del_ymph-del_ymnh)/(2*h)

    return der_delym

###############################################################################

def Newton_Rapshonyp(tilt_1, set_rayos):
    tip_1 = (0.0,0.0)
    tilt_1 = np.array(tilt_1)
    tip_1 = np.array(tip_1)
    numcnt = 50000
    tip_01 = 0.0
    cnt = 0
    while True:
      surf_1, plane_1 = rotat_trasl(tip_1,tilt_1,set_rayos)
      l_surface = plane_1 [2]
      del_ym = Height_yp(l_surface)
      Der_ym = Der_highyp(tilt_1, set_rayos)
      tip_1[0] = tip_01 - del_ym/Der_ym
      print(del_ym)
      print(del_ym/Der_ym)
      print(tip_01 - del_ym/Der_ym)
      if del_ym < 1.5E-3:
          break
      else: 
          tip_01 = tip_1[0]
      if cnt == numcnt:
        break
      cnt = cnt + 1

    Tx = np.deg2rad(tip_1[0])
    pointzcoma = np.abs(tilt_1[1]/np.sin(Tx))

    return tip_1, pointzcoma

###########################################################################

def setof_rayos(option):
  option = option
  if option == 0:

    '''
    Parameters for surfaces
    '''
    NR_1 = 100
    LL_1 = 1.071721E+003

    sistx1 = np.linspace(-LL_1, LL_1, num=NR_1)
    sisty1 = np.random.randint(-LL_1,LL_1,NR_1)*0.0
    sistz1 = np.ones(NR_1) * -100
    sistL1 = np.zeros(NR_1)
    sistM1 = np.zeros(NR_1)
    sistN1 = np.ones(NR_1)

    conrayos_surf = []
    conrayos_surf.append(sistx1)
    conrayos_surf.append(sisty1)
    conrayos_surf.append(sistz1)
    conrayos_surf.append(sistL1)
    conrayos_surf.append(sistM1)
    conrayos_surf.append(sistN1)

    '''
    Parameters for rays
    '''
    meshcoorx = ( 0.0, 0.0, 0.0)
    meshcoory = ( 0.0, 1.071721E+003, -1.071721E+003)
    meshcoorx = np.array(meshcoorx)
    meshcoory = np.array(meshcoory)
    ###############################################################################

    NR_1 = len(meshcoorx)
    sistz = np.ones(NR_1) * -100
    sistL = np.zeros(NR_1)
    sistM = np.zeros(NR_1)
    sistN = np.ones(NR_1) 

    conrayos_ray = []
    conrayos_ray.append(meshcoorx)
    conrayos_ray.append(meshcoory)
    conrayos_ray.append(sistz)
    conrayos_ray.append(sistL)
    conrayos_ray.append(sistM)
    conrayos_ray.append(sistN)
  
  return(conrayos_surf, conrayos_ray)

###############################################################################


def setof_rayos1(option):
  option = option
  if option == 0:
    '''
    Parameters for surfaces
    '''
    NR_1 = 100
    LL_1 = 1.071721E+003

    sistx1 = np.linspace(-LL_1, LL_1, num=NR_1)
    sisty1 = np.random.randint(-LL_1,LL_1,NR_1)*0.0
    sistz1 = np.ones(NR_1) * -100
    sistL1 = np.zeros(NR_1)
    sistM1 = np.zeros(NR_1)
    sistN1 = np.ones(NR_1)

    conrayos_surf = []
    conrayos_surf.append(sistx1)
    conrayos_surf.append(sisty1)
    conrayos_surf.append(sistz1)
    conrayos_surf.append(sistL1)
    conrayos_surf.append(sistM1)
    conrayos_surf.append(sistN1)

    '''
    Parameters for rays
    '''

    num_0 = 20
    LL = 2 * 1.071721E+003
    rad = []
    NR = 10
    sistx = []
    sisty = []
    NR = []
    num_i = 6.0

    div_0 = LL/num_0


    for i in range(11):
      rad.append(LL-i*div_0)
      NR.append(num_i * i)

    NR = np.array(NR)
    rad = np.array(rad)
    rad = -(rad - LL)
    #############################################################################
    for i in range(11):
      num_segmentos = int(NR[i])
      angulo = np.linspace(0, 2*np.pi, num_segmentos+1)
      sistx.append(rad[i] * np.cos(angulo))
      sisty.append(rad[i] * np.sin(angulo))

    #############################################################################
    meshcoorx = []
    meshcoory = []

    for i in range(11):
      numsuma = int(NR[i]+1)
      for j in range(numsuma):
        meshcoorx.append(sistx[i][j])
        meshcoory.append(sisty[i][j])

    meshcoorx = np.array(meshcoorx) 
    meshcoory = np.array(meshcoory)
    ###############################################################################

    NR_1 = len(meshcoorx)
    sistz = np.ones(NR_1) * -100
    sistL = np.zeros(NR_1) 
    sistM = np.zeros(NR_1) 
    sistN = np.ones(NR_1) 

    conrayos_ray = []
    conrayos_ray.append(meshcoorx)
    conrayos_ray.append(meshcoory)
    conrayos_ray.append(sistz)
    conrayos_ray.append(sistL)
    conrayos_ray.append(sistM)
    conrayos_ray.append(sistN)

  else:
    '''
    Parameters for surfaces
    '''
    NR_1 = 100
    LL_1 = 1.071721E+003

    sistx1 = np.linspace(-LL_1, LL_1, num=NR_1)
    sisty1 = np.random.randint(-LL_1,LL_1,NR_1)*0.0
    sistz1 = np.ones(NR_1) * -100
    sistL1 = np.zeros(NR_1) 
    sistM1 = np.zeros(NR_1)
    sistN1 = np.ones(NR_1)

    conrayos_surf = []
    conrayos_surf.append(sistx1)
    conrayos_surf.append(sisty1)
    conrayos_surf.append(sistz1)
    conrayos_surf.append(sistL1)
    conrayos_surf.append(sistM1)
    conrayos_surf.append(sistN1)

    '''
   Parameters for rays
    '''

    NR = 45
    LL = 2 * 1.071721E+003

    R = []
    teta = []

    sistx = np.linspace(LL, -LL, num=NR)
    sisty = np.linspace(LL, -LL, num=NR)
    sistxx, sistyy = np.meshgrid(sistx, sisty)
    R = np.sqrt(sistxx**2 + sistyy**2)
    R[(R>=1.071721E+003)] = np.NaN
    cos = sistxx/R
    sin = sistyy/R
    teta = np.where( sin >= 0. , np.arccos(cos) , -np.arccos(cos) )
    teta = np.where(R == 0. , 0., teta)
    systxx = R*np.cos(teta)                                                             
    systyy = R*np.sin(teta) 

    meshcoorx = []
    meshcoory = []

    for i in range(NR):
      for j in range(NR):

        meshcoorx.append(systxx[i][j])
        meshcoory.append(systyy[i][j])


    meshcoorx = np.array(meshcoorx)
    meshcoory = np.array(meshcoory)
    meshcoorx = meshcoorx[~np.isnan(meshcoorx)]
    meshcoory = meshcoory[~np.isnan(meshcoory)]
    NR_1 = len(meshcoorx)
    sistz = np.ones(NR_1) * -100
    sistL = np.zeros(NR_1) 
    sistM = np.zeros(NR_1)
    sistN = np.ones(NR_1) 

    conrayos_ray = []
    conrayos_ray.append(meshcoorx)
    conrayos_ray.append(meshcoory)
    conrayos_ray.append(sistz)
    conrayos_ray.append(sistL)
    conrayos_ray.append(sistM)
    conrayos_ray.append(sistN)
  
  return(conrayos_surf, conrayos_ray)

###############################################################################
###############################################################################
###############################################################################

cosL = np.zeros(10) * 0.0
cosM = np.zeros(10) * 0.0
cosN = np.zeros(10) * 0.0
coorx_0 = np.ones(10) * 1.0
coory_0 = np.ones(10) * 1.0
coorz_0 = np.ones(10) * 1.0
coorx_1 = np.ones(30) * 1.0
coory_1 = np.ones(30) * 1.0
coorz_1 = np.ones(30) * 1.0
Rc = np.zeros(10)
kk = np.zeros(10)

Rc[0] = -7670.112
Rc[1] = -7670.112
Rc[2] = 9999999999999999.
Rc[3] = 9999999999999999.
Rc[4] = 9999999999999999

d = np.zeros(10) * 0.0
d[0] = 0.0 
d[1] = - 2272.41
d[2] = -(- 2272.41) + 56.538
d[3] = 7.0   
d[4] = 300.0

dxps = np.zeros(10) * 0.0
dxps[4] =  0.0

dyps = np.zeros(10) * 0.0
dyps[4] =  0.0 

kk[0] = -1.597
kk[1] = -37.027
kk[2] = 0.0
kk[3] = 0.0
kk[4] = 0.0

aalpha_1 = np.zeros(10) * 0.0
aalpha_2 = np.zeros(10) * 0.0
aalpha_3 = np.zeros(10) * 0.0
aalpha_4 = np.zeros(10) * 0.0
aalpha_5 = np.zeros(10) * 0.0 

aalpha_3[1] = 5.943E-019
aalpha_1[2] = 3.955E-005
aalpha_2[2] = -1.021E-009

nn0 = np.ones(10) * 1.0 
nn1 = np.ones(10) * -1.0 


lmd=0.4

K_1 = 6.69422575E-001
L_1 = 4.48011239E-003
K_2 = 4.34583937E-001
L_2 = 1.32847049E-002
K_3 = 8.71694723E-001
L_3 = 9.53414824E+001

nn1[2] = np.sqrt(1+(K_1*lmd**2)/(lmd**2-L_1)+(K_2*lmd**2)/(lmd**2-L_2)+(K_3*lmd**2)/(lmd**2-L_3))
nn0[3] = np.sqrt(1+(K_1*lmd**2)/(lmd**2-L_1)+(K_2*lmd**2)/(lmd**2-L_2)+(K_3*lmd**2)/(lmd**2-L_3))
nn1[3] = 1.0

##############################################################################

tilt_m = (0.0, 1.0)
tip_m = (0.0, 0.0)
opc = 0

surf_rset , ray_rset = setof_rayos(opc)
surf_surf, surf_plane = rotat_trasl(tip_m, tilt_m, ray_rset)
ray_surf, ray_plane = rotat_trasl(tip_m,tilt_m, ray_rset)

tip_cm, pointncoma = Newton_Rapshonyp(tilt_m, ray_rset)

surf_cset , ray_cset = setof_rayos1(opc)
surf_msurf, surf_mplane = rotat_trasl(tip_m, tilt_m, ray_cset)
surf_cmsurf, surf_cmplane = rotat_trasl(tip_cm, tilt_m, surf_cset)
ray_cmsurf, ray_cmplane = rotat_trasl(tip_cm,tilt_m, ray_cset)

print(tip_cm, pointncoma)

plt.plot(surf_mplane[2][0], surf_mplane[2][1],'.', color='k', linestyle='none')
plt.axis('equal')
plt.xlabel('X')
plt.title('Spot diagram')

plt.show()

plt.plot(ray_cmplane [2][0], ray_cmplane[2][1],'.', color='k', linestyle='none')
plt.title('Corrected Spot diagram ')
plt.axis('equal')

plt.show()

plt.plot(ray_cmsurf[0][0], ray_cmsurf[0][1],'.', color='k', linestyle='none')
plt.axis('equal')
plt.ylabel('Y')
plt.title('Primary mirror')

plt.show()