import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

 
'''
Developer: Najera Roa, Morgan Rhai.
email: mnajera@astro.unam.mx
orcid: http://orcid.org/0000-0003-3283-0407 
'''


'''
This example uses the angle calculated using the SoS_ZCP_Non_classical_telescope_example.py
to visualize nine different fields where the optic system has been compensated. We also
show two additional plots, one displaying the optical system perfectly aligned and 
the other with the secondary mirror shifted. Finally, one last plot is show of the compensated out-of-focus 
optical system. This example require to define a new function described in Section c.7) and we 
rewrite function n.1).


'''


'''

Telescope's parameters

\begin{deluxetable*}{ccccccccc}
\label{tab:tabletelescope}
\tablecaption{Telescope $\phi$~2.1~m  $f/3$ (Theoretical example)}
\tablewidth{0pt}
\tablehead{
\colhead{Element} & \colhead{Radius} & \colhead{Thickness} & \colhead{Glass} & \colhead{Semi-Diameter} & \colhead{$k$} & \colhead{$\alpha_1$} & \colhead{$\alpha_2$}& \colhead{$\alpha_3$}\\
\colhead{ } & \colhead{(mm)} & \colhead{(mm)} & \colhead{} & \colhead{(mm)} & \colhead{} & \colhead{} & \colhead{}& \colhead{}
}
\startdata
$M_1$ & $-7670.112$ &$-227241$  &Mirror &1005.0 &$-1.597$       & & & \\
$M_2$ & $-7670.112$ &227241   &Mirror &465.0    &$-37.027$  & & &5.940e$-19$  \\
$M_1$ Vertex    &   &56.538   & Air   &     &           & & & \\
Corrector Front &   &7.0    &Silica Schott  & 191.0 &   &3.955e$-5$ &$-$1.021e$-9$ &\\
Corrector Back  &   &300.0    &Air  & 191.0 &   &    &  & \\
Image Plane &   &     &                &    78.4    &           & & & \\
\enddata
\end{deluxetable*}
'''

'''
Part I.
General exact ray tracing
 
'''
 
 
###############################################################################
 
'''
 
a) Parameters of the Optical Elements
 
In this Section we define some parameters for the optical system, such as:
 
1.- Number of surfaces: this is fundamental in accordance with the general exact
ray tracing procedure. This procedure uses vectors to represent the rays of light,
hence, we know the coordinates of the origin point (P_{-1}(x_{-1}, y_{-1}, z_{-1})
of each ray and its direction. Subsequently, we found the intersection point on 
the next surface defined by the sagitta equation (see Section c.2), which is centered 
on a new coordinate system positioned at a distance d; however, it is necessary to 
to first find the intersection with the coordinate system (P_{0}(x_{0}, y_{0}, z_{0}). 
If we already know the intersection points with the i-th surface, the procedure is 
the same to get the intersection points of the i-th + 1 surface.
 
2.- Radius of curvature of the different surfaces that we are going to use in the
sagitta equation.
 
3.- Conic constant for the surfaces: this constant, as well as the radius of curvature, 
is going to be used in the sagitta equation.
 
4.- Semi-diameters of the primary and the secondary mirrors to define the set of 
the rays that we will intersect with the primary mirror. Furthermore, given its
position, the semi-diameter of the secondary mirror defines a surface that is blocking 
some rays of our original set of rays.
 
5.- Distance between surfaces: as we mentioned, every surface is positioned at a
distance d from each other.
 
6.- Asphericity polynomial terms: a non-classical telescope is composed of conical 
surfaces (as the sagitta equation) plus asphericity polynomial terms; these terms 
allow one or more aspherical surface or correcting lenses to be contemplated in 
the optical system. For classical telescopes, the asphericity polynomial terms 
are equal to zero.
 
7.- Wavelength: This parameter allows us to use differents wavelengths according 
to our requirements.
 
8.- Refractive index: In the case of classical telescopes the refractive index 
is specified as 1 and -1, because we are working only with mirrors. However,
we can change the refractive index if we have correcting aspherical lenses. To 
define the refractive index, you have to set the dispersion coefficients.
 
'''
 
 
 
'''
Parameters of the optical elements and the rays
'''
 
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

 
'''
1a.- Number of surfaces of the optical system 
'''
 
number_of_surfaces = 5
 
'''
2a.- Radius of curvature of the different surfaces
'''
 
Rc[0] = -7670.112
Rc[1] = -7670.112
Rc[2] = 9999999999999999.
Rc[3] = 9999999999999999.
Rc[4] = 9999999999999999
 
"""
3a.- Conic constant for the surfaces
"""
 
kk[0] = -1.597
kk[1] = -37.027
kk[2] = 0.0
kk[3] = 0.0
kk[4] = 0.0
 
"""
4a.- Semi-diameters of the primary mirror and the secondary mirror
"""
 
semidiameter_pm = 1.071721E+003
semidiameter_sm = 4.64365E+002
 
'''
5a.- Distance between surfaces
'''
d = np.zeros(10) * 0.0
d[0] = 0.0 
d[1] = - 2272.41
d[2] = -(- 2272.41) + 56.538
d[3] = 7.0   
d[4] = 300.0
"""
6a.- Asphericity polynomial terms
"""
 
aalpha_1 = np.zeros(10) * 0.0
aalpha_2 = np.zeros(10) * 0.0
aalpha_3 = np.zeros(10) * 0.0
aalpha_4 = np.zeros(10) * 0.0
aalpha_5 = np.zeros(10) * 0.0 
 
aalpha_3[1] = 5.943E-019
aalpha_1[2] = 3.955E-005
aalpha_2[2] = -1.021E-009
 
'''
7a.- Wavelength
'''
lmd=0.4
 
"""
8a.- Refractive index
"""
nn0 = np.ones(10) * 1.0 
nn1 = np.ones(10) * -1.0 
 
"""
Coefficients for an specific glass: Silica SCHOOTT
"""
K_1 = 6.69422575E-001
L_1 = 4.48011239E-003
K_2 = 4.34583937E-001
L_2 = 1.32847049E-002
K_3 = 8.71694723E-001
L_3 = 9.53414824E+001
 
"""
Dispersion formula:
In this special case, we use the Sellmeier formula for an specific wavelength 
"""
 
nn1[2] = np.sqrt(1+(K_1*lmd**2)/(lmd**2-L_1)+(K_2*lmd**2)/(lmd**2-L_2)+(K_3*lmd**2)/(lmd**2-L_3))
nn0[3] = np.sqrt(1+(K_1*lmd**2)/(lmd**2-L_1)+(K_2*lmd**2)/(lmd**2-L_2)+(K_3*lmd**2)/(lmd**2-L_3))
nn1[3] = 1.0
 
##############################################################################
##############################################################################
print('This process may run slower depending on your hardware, please wait...')
 
"""
 
b) Ray parameters
 
Here we are going to focus on defining two sets of rays, where each one will contain 
two subsets. One subset will define the surfaces and the second subset the grid-mesh.
We have two kinds of mesh, the first is an hexapolar mesh and the other is a cartesian 
circular mesh.
 
We are going to briefly introduce some important variables that are used in this Section. 
First, the variable called option, allows us to choose what kind of subset we want, if 
option is equal to 0, then we are going to work with the hexapolar mesh. On the other 
hand, if option is equal to 1 it will use the circular mesh. NR and NR_1 set the number 
of elements that we want to use in our sets. LL and LL_1 define the range of the points, 
for this reason, we chose a semidiameter_pm (semidiameter of the primary mirror) as the 
maximum value. For the hexapolar mesh we have num_i, which indicates the number of multiple 
points of that parameter. 

Moreover, rad and R define the radius of the surface which, in the same way as before, 
the values of this two variables are going to be set between semidiameter_pm and 
semidiameter_sm. In order to calculate the coordinates of the origin point of each 
ray and its direction, the variables meshcoorx, meshcoory, systz, systL, systM 
and systN are used to define the value of (P_{-1}(x_{-1}, y_{-1}, z_{-1}). Finally, 
the two subsets will be saved in assemblagesurf_surf and assemblageray_ray.
 
"""
 
###############################################################################
 
def setof_rayos1(option):
  option = option
  if option == 0:
    '''
    Parameters for surfaces
    '''
    NR_1 = 100
    LL_1 = semidiameter_pm 
 
    systx1 = np.linspace(-LL_1, LL_1, num=NR_1)
    systy1 = np.random.randint(-LL_1,LL_1,NR_1)*0.0
    systz1 = np.ones(NR_1) * -100
    systL1 = np.zeros(NR_1)
    systM1 = np.zeros(NR_1)
    systN1 = np.ones(NR_1)
 
    assemblagesurf_surf = [ ]
    assemblagesurf_surf.append(systx1)
    assemblagesurf_surf.append(systy1)
    assemblagesurf_surf.append(systz1)
    assemblagesurf_surf.append(systL1)
    assemblagesurf_surf.append(systM1)
    assemblagesurf_surf.append(systN1)
 
    '''
    Parameters for rays
    '''
 
    num_0 = 10
    LL = semidiameter_pm
    rad = []
    NR = 10
    systx = []
    systy = []
    NR = []
    num_i = 6.0
 
    div_0 = LL/num_0
 
 
    for i in range(11):
      rad.append(LL-i*div_0)
      NR.append(num_i * i)
 
    NR = np.array(NR)
    rad = np.array(rad)
    rad = -(rad-semidiameter_pm)
    rad[(rad<=semidiameter_sm)] = np.NaN
    rad = rad[~np.isnan(rad)]
    numrad = len(rad)
    numrad1 = numrad-1
    for i in range(numrad1):
      NR=np.delete(NR,abs(i-i))
 
 
    #############################################################################
    for i in range(numrad):
      num_segmentos = int(NR[i])
      angle = np.linspace(0, 2*np.pi, num_segmentos+1)
      systx.append(rad[i] * np.cos(angle))
      systy.append(rad[i] * np.sin(angle))
 
   #############################################################################
    meshcoorx = []
    meshcoory = []
 
    for i in range(numrad):
      numsum = int(NR[i]+1)
      for j in range(numsum):
        meshcoorx.append(systx[i][j])
        meshcoory.append(systy[i][j])
 
    meshcoorx = np.array(meshcoorx) 
    meshcoory = np.array(meshcoory)
###############################################################################
 
    NR_1 = len(meshcoorx)
    systz = np.ones(NR_1) * -100
    systL = np.zeros(NR_1) 
    systM = np.zeros(NR_1) 
    systN = np.ones(NR_1) 
 
    assemblageray_ray = []
    assemblageray_ray.append(meshcoorx)
    assemblageray_ray.append(meshcoory)
    assemblageray_ray.append(systz)
    assemblageray_ray.append(systL)
    assemblageray_ray.append(systM)
    assemblageray_ray.append(systN)
 
  else:
    '''
    Parameters for surfaces
    '''
    NR_1 = 100
    LL_1 = semidiameter_pm
 
    systx1 = np.linspace(-LL_1, LL_1, num=NR_1)
    systy1 = np.random.randint(-LL_1,LL_1,NR_1)*0.0
    systz1 = np.ones(NR_1) * -100
    systL1 = np.zeros(NR_1) 
    systM1 = np.zeros(NR_1)
    systN1 = np.ones(NR_1)
 
    assemblagesurf_surf = [ ]
    assemblagesurf_surf.append(systx1)
    assemblagesurf_surf.append(systy1)
    assemblagesurf_surf.append(systz1)
    assemblagesurf_surf.append(systL1)
    assemblagesurf_surf.append(systM1)
    assemblagesurf_surf.append(systN1)
 
    '''
   Parameters for rays
    '''
 
    NR = 45
    LL = 2 * semidiameter_pm
 
    R = []
    teta = []
 
    systx = np.linspace(LL, -LL, num=NR)
    systy = np.linspace(LL, -LL, num=NR)
    systxx, systyy = np.meshgrid(systx, systy)
    R = np.sqrt(systxx**2 + systyy**2)
    R[(R<=semidiameter_sm)] = np.NaN
    R[(R>=semidiameter_pm )] = np.NaN
    cos = systxx/R
    sin = systyy/R
    theta = np.where( sin >= 0. , np.arccos(cos) , -np.arccos(cos) )
    theta = np.where( R == 0. , 0., theta)
    systxx = R*np.cos(theta)                                                             
    systyy = R*np.sin(theta) 
 
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
    systz = np.ones(NR_1) * -100
    systL = np.zeros(NR_1) 
    systM = np.zeros(NR_1)
    systN = np.ones(NR_1) 
 
    assemblageray_ray = []
    assemblageray_ray.append(meshcoorx)
    assemblageray_ray.append(meshcoory)
    assemblageray_ray.append(systz)
    assemblageray_ray.append(systL)
    assemblageray_ray.append(systM)
    assemblageray_ray.append(systN)
  
  return (assemblagesurf_surf, assemblageray_ray)
 
###############################################################################
###############################################################################
 
"""
 
c) Find the intersection of the ray in the next coordinate system
 
As we mentioned in Section a), the main goal in general exact ray tracing is to find
the intersection of the ray in every coordinate system, because it determines the
optical path of the ray, passing through a set of surfaces. As a first step, we
define a function that will help us to find the intersection point. Hereafter this 
Section is divided into Subsections, because these subsections contain the functions 
that we will use. We present a brief explanation of each function below.  
 
"""
 
###############################################################################
 
'''
 
c.1) Coordinates of the P_{0} function
 
Before we find the intersection point P_{1}, which is the intersection point on 
a surface, we need to find the origin point P_{0}. To find this point, we use 
equation x_1 = (L/N)*(z_{1} - z_{0}) + x_{0} and y_{1} = (M/N)*(z_{1} - z_{0}) + y_{0}.
Then, we iteratively solve the ray-surface function.
 
'''
 
def F_RS(z_1,rayorigin,i):
    
    
    x_0 = rayorigin[0]
    y_0 = rayorigin[1]
    z_0 = rayorigin[2]
    L = rayorigin[3]
    M = rayorigin[4]
    N = rayorigin[5]
    
    
    x_1 = (L/N)*(z_1-z_0)+x_0
    y_1 = (M/N)*(z_1-z_0)+y_0
    
    f_RS = surface(x_1,y_1,z_1,i)    
    
    return f_RS
       
    
###########################################################################
 
'''
 
c.2) Ray-Surface function
 
After finding the coordinates x_1 and y_1, we iteratively solve the ray-surface
function. This function contains the sagitta function which defines a surface. 
 
'''
 
def surface(x_1,y_1,z_1,i):
    
    R_c = Rc[i]
    k = kk[i]
    Alpha_1 = aalpha_1[i]
    Alpha_2 = aalpha_2[i]
    Alpha_3 = aalpha_3[i]
    Alpha_4 = aalpha_4[i]
    Alpha_5 = aalpha_5[i]
 
  
    c = 1.0/R_c
    s2 = (x_1*x_1)+(y_1*y_1)
    s4 = s2*s2
    s6 = s4*s2
    s8 = s6*s2
    s10 = s8*s2
    
    f1 = c*s2/(1.0+np.sqrt(1.0-((k+1.0)*c*c*s2)))
    f2 = (Alpha_1*s2)+(Alpha_2*s4)+(Alpha_3*s6)+(Alpha_4*s8)+(Alpha_5*s10)
 
    return f1+f2-z_1
  
 
###########################################################################
 
'''
 
c.3) Derivative of the ray-surface function
 
We compute a numerical derivative of the ray-surface function, instead of perfor-
ming a simple derivative, in order to calculate the intersection point.
 
'''
 
def Der_F_RS(z_1,rayorigin,i):
    h = 0.00000000001
    der = (F_RS(z_1+h, rayorigin,i)-F_RS(z_1-h, rayorigin,i))/(2*h)
    return der
 
###########################################################################
 
'''
 
c.4) Newton_Raphson function
 
To find the point P_{1}, we use the Newton-Rhapson method within a loop. This method
is useful for finding the roots of the function that defines the intersection with
the surface (c.3).
 
'''
 
 
def Newton_Raphson(rayorigin,i):
    z_1 = 0
    cnt = 0
    while True:
        z_2 = z_1-(F_RS(z_1,rayorigin,i)/Der_F_RS(z_1,rayorigin,i))
        
        if np.abs(z_2-z_1)<0.00001:
            
            break
        else:
            z_1=z_2
        if cnt == 30:
            break    
        cnt=cnt+1    
 
        
    [x_0, y_0, z_0,L,M,N]=rayorigin
    x_1 = (L/N)*(z_1-z_0)+x_0
    y_1 = (M/N)*(z_1-z_0)+y_0
        
    return x_1, y_1, z_1
 
 
################################################################################
 
'''
 
c.5) Normal Vector function
 
Once the position of the intersection point is known, we calculate the normal 
vector using the numerical derivative of the ray-surface function.
 
'''
 
def surface_Derivative(x,y,z,i):
    delta=0.000001
    
    Dx = (surface(x+delta,y,z,i)-surface(x-delta,y,z,i))/(2.0*delta)
    Dy = (surface(x,y+delta,z,i)-surface(x,y-delta,z,i))/(2.0*delta)
    Dz = (surface(x,y,z+delta,i)-surface(x,y,z-delta,i))/(2.0*delta)   
    
    Dr = np.sqrt((Dx*Dx)+(Dy*Dy)+(Dz*Dz))
    return Dx/Dr,Dy/Dr,Dz/Dr
 
 
###########################################################################
 
'''
 
c.6) Snell's function
 
When the ray hits a surface, we already know the direction cosines and the normal 
vector to that surface, but we don't know the direction of the ray that is refracted
or reflected. To know that we have to use Snell's law in its vector form.
 
 
'''
 
def Snell_refraction_vector(Ray_Vect,Surf_Normal,n_0,n_1):
    
        
    Ray_Vect = np.asarray(Ray_Vect)
    Surf_Normal = np.asarray(Surf_Normal)       
    
    Nsurf_Cros_s1 = np.cross(Surf_Normal,Ray_Vect)
        
    NN = n_0/n_1
    
    R = (NN*NN)*np.dot(Nsurf_Cros_s1,Nsurf_Cros_s1)
            
    if R>1:
        print("---- internal reflection= ", R)
        
    s2 = NN*(np.cross(Surf_Normal,np.cross(-Surf_Normal,Ray_Vect)))-Surf_Normal*np.sqrt(1.0-((NN*NN)*np.dot(Nsurf_Cros_s1,Nsurf_Cros_s1)))
   
    return s2
 
 
###########################################################################
 
'''
 
c.7) Skew ray tracing function
 
This function uses all the previous functions to calculate the point of intersection
and the direction of the ray when it interacts with the surface. This function is
used as many times as there are surfaces.
 
'''
 
def Trax(Ray_Param,i):
    ## Point of intersection on the surface 
    x_1, y_1, z_1 = Newton_Raphson(Ray_Param,i)
    
    ## Normal vector at the point of intersection on the surface 
    NORM = surface_Derivative(x_1,y_1,z_1,i)
    [x_0, y_0, z_0,L,M,N] = Ray_Param
    Ray_Vect = [L,M,N]
    Surf_Normal = NORM
    
    
    
    n_0 = nn0[i]
    n_1 = nn1[i]
        
    lmn = Snell_refraction_vector(Ray_Vect,Surf_Normal,n_0,n_1)
    
    return ([x_1, y_1, z_1, lmn[0], lmn[1], lmn[2]])
 
###########################################################################
 
'''
 
c.6) Director Cosines Function
 
If we don't know the director cosines of some point, we can calculate them with this function.
 
'''
 
def director_cosines(x,y,z):
  magnitude = np.sqrt(x**2+y**2+z**2)
  L = x/magnitude
  M = y/magnitude
  N = z/magnitude
 
  return L, M, N

###############################################################################
 
'''
 
c.7) Arrangement of fields function.
 
This function helps to know the position and the direction of determined arrangement 
of fields.
 
'''

def arrangement_ray(option, angle, amount, points):
    option = option
    LL = angle * (np.pi/180)
    Arc_l = LL * 0.1071721
    cnt = amount
    num_i = points
    if option == 0:
        point_01 = np.linspace(-Arc_l, Arc_l, num=cnt)
        point_02 = np.linspace(-Arc_l, Arc_l, num=cnt)
        num_point = len(point_01)


        point_001,point_002 = np.meshgrid(point_01, point_02)
        meshcoorx = []
        meshcoory = []
################################################################################
        for i in range(num_point):
          for j in range(num_point):

            meshcoorx.append(point_001[i][j])
            meshcoory.append(point_002[i][j])

        meshcoorx = np.array(meshcoorx)
        meshcoory = np.array(meshcoory)
        lencos = len(meshcoorx)
        L_1 = np.ones(lencos)
        M_1 = np.ones(lencos)
        N_1 = np.ones(lencos)
        for i in range(lencos):
            L_1[i], M_1[i], N_1[i] = director_cosines( meshcoorx[i], meshcoory[i], 1)

        cosdir_pack = (L_1,M_1,N_1) 
        
    elif option == 1:
        num_0 = 10
        rad = []
        NR = 10
        systx = []
        systy = []
        NR = []


        div_0 = LL/num_0



        for i in range(amount):
          rad.append(LL-i*div_0)
          NR.append(num_i * i)

        NR = np.array(NR)
        rad = np.array(rad)
        rad = -(rad-LL)
#############################################################################
        for i in range(amount):
          num_segmentos = int(NR[i])
          angle_0 = np.linspace(0, 2*np.pi, num_segmentos+1)
          systx.append(rad[i] * np.cos(angle_0))
          systy.append(rad[i] * np.sin(angle_0))
#############################################################################
        meshcoorx = []
        meshcoory = []

        for i in range(amount):
          numsume = int(NR[i]+1)
          for j in range(numsume):
            meshcoorx.append(systx[i][j])
            meshcoory.append(systy[i][j])

        meshcoorx = np.array(meshcoorx) 
        meshcoory = np.array(meshcoory)
        lencos = len(meshcoorx)
        L_1 = np.ones(lencos)
        M_1 = np.ones(lencos)
        N_1 = np.ones(lencos)
        for i in range(lencos):
            L_1[i], M_1[i], N_1[i] = director_cosines( meshcoorx[i], meshcoory[i], 1)

        cosdir_pack = (L_1,M_1,N_1) 

    return (meshcoorx,meshcoory, cosdir_pack)
 
 
################################################################################
 
'''
d) Change the ray to the transformed coordinate system
 
We want to perform a displacement and a rotation on the secondary mirror. For this 
we calculate the values of a segment of a ray that travels between two surfaces. 
The segment is defined by two extreme points, P_0 and P_1, both have a coordinate 
system (x,y,z). The shift and rotation on the secondary mirror make the original 
coordinate system to become (x',y',z'). Then, we calculate the intersection point 
with the surface in the new system that has undergone the transformation.
 
The information of the point of intersection in the new system can be obtained with 
different functions that we will explain below.
 
'''
 
################################################################################
################################################################################
 
'''
 
Transformation function.
 
The information in the transformed space is known through the translation 
matrix (T_ {xyz}) and the rotation matrices (R_x, R_y and R_z). The first matrix 
contains the three displacements of the three spatial directions and the other 
matrices contain the angles at which the rotation around the axes was applied.
 
Consequently, the matrix with the resulting transformation is the product of these
four matrices. Also, this matrix will involve an inverse transformation to return
to the original coordinate system.
 
 
'''
 
def Transform(TRF):
    [TiltX, TiltY, TiltZ, dx, dy, dz]=TRF
    Tx=np.deg2rad(TiltX)
    Ty=np.deg2rad(TiltY)
    Tz=np.deg2rad(TiltZ)    
    
        
    #https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    Rx_1A = np.matrix([[1.0, 0.0, 0.0 , 0.0], [0.0, np.cos(-Tx), np.sin(-Tx),0.0],[0.0,-np.sin(-Tx), np.cos(-Tx),0.0],[0.0, 0.0, 0.0 ,1.0]])
    Ry_1A = np.matrix([[np.cos(-Ty),0.0, -np.sin(-Ty),0.0], [0, 1.0, 0.0, 0.0],[np.sin(-Ty),0.0, np.cos(-Ty),0.0],[0.0, 0.0, 0.0 ,1.0]])
    Rz_1A = np.matrix([[np.cos(-Tz),-np.sin(-Tz), 0.0, 0.0],[np.sin(-Tz),np.cos(-Tz), 0.0, 0.0],[0.0, 0.0, 1.0, 0.0],[0, 0, 0 ,1.0]])    
    Dxyz_1A = np.matrix([[1.0, 0.0, 0.0, -dx],[0.0, 1.0, 0.0, -dy],[0.0, 0.0, 1.0, -dz],[0.0, 0.0, 0.0 ,1.0]])
    
    Start_trans1=Rz_1A*Ry_1A*Rx_1A*Dxyz_1A
    Start_trans2=np.linalg.inv(Start_trans1)
    
        
    return Start_trans1, Start_trans2
 
 
###############################################################################
 
'''
Operations in the transformed space system
 
'''
 
def Traceray(inx,iny,inz,inL,inM,inN, MT, dfocus):
    
 
    coorx_0[0] = inx
    coory_0[0] = iny
    coorz_0[0] = inz
    cosL[0] = inL
    cosM[0] = inM
    cosN[0] = inN
    arraydfocus = np.zeros(number_of_surfaces)
    arraydfocus[4] = dfocus
    
    for i in range (number_of_surfaces):
 
        
        '''
        A1_Trans is the transformation matrix to the transformed space and A2_Trans
        contains the inverse transformation.
        '''
 
        A1_Trans = MT[i][0]
        A2_Trans = MT[i][1]
        
        x_0 = coorx_0[i] 
        y_0 = coory_0[i] 
        z_0 = coorz_0[i] - d[i] - arraydfocus[i]
        L = cosL[i]
        M = cosM[i]
        N = cosN[i]
 
        
        '''
        e), f) and g) We already know the points of origin and their direction, thus, 
        we can calculate the second point, in order to find the extreme points P_0 and P_1.
        '''
        
 
        Ray_Param = [x_0, y_0, z_0, L, M, N]
        xyz_1=Trax(Ray_Param,i)
        
        # #################################################
 
        x_1 = xyz_1[0]
        y_1 = xyz_1[1]
        z_1 = xyz_1[2]
        
       
        '''
        We transform the components of P_0 and P_1 using the transformation matrix
        as follows:
        '''
 
 
        prime_coor0 = np.dot(A1_Trans, ([x_0, y_0, z_0, 1.0]))
        prime_coor1 = np.dot(A1_Trans,([x_1, y_1, z_1, 1.0]))
 
        '''
        We have to calculate the director cosines in the transformed space because
        the director cosines are not transformed automatically.
        '''
 
        prime_scosinD = Cosine_transform(prime_coor0, prime_coor1) 
 
        '''                                
        Now, we know the points of origin in the transformed space and their direction:   
        '''
 
        x_0 = prime_coor0[0,0]
        y_0 = prime_coor0[0,1]
        z_0 = prime_coor0[0,2]
        
        
        s=np.sign(d[i])
        
        if i==0:
            s=1
            
        L = prime_scosinD[0]*s
        M = prime_scosinD[1]*s
        N = prime_scosinD[2]*s
      
        
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
 
        
        '''
        h) At this point we know the point of intersection with the system that underwent 
        a displacement and a rotation. We have to calculate two points because we 
        need to transform the director cosines to the original system. 
 
        j) In this way, prime_coor0 and prime_coor1 are the intersection points in the 
        original system. 
 
        k) Prime_scosinD are the director cosines of the ray in the original system.
        '''
        
        
        prime_coor0 = np.dot(A2_Trans, ([x_1, y_1, z_1, 1.0]))
        prime_coor1 = np.dot(A2_Trans, ([x_2, y_2, z_2, 1.0]))
        prime_scosinD = Cosine_transform(prime_coor0,prime_coor1)  
        
        x_1 = prime_coor0[0,0]
        y_1 = prime_coor0[0,1]
        z_1 = prime_coor0[0,2]
        
        L = prime_scosinD[0]
        M = prime_scosinD[1]
        N = prime_scosinD[2]
        
        
        coorx_1[i] = x_1
        coory_1[i] = y_1
        coorz_1[i] = z_1
        cosL[i+1] = L
        cosM[i+1] = M
        cosN[i+1] = N
        
        coorx_0[i+1] = coorx_1[i]
        coory_0[i+1] = coory_1[i]
        coorz_0[i+1] = coorz_1[i]
 
    return coorx_1,coory_1, coorz_1, cosL, cosM, cosN
 
##############################################################################
 
'''
 
Cosine transformation function.
 
As we mention, we have to calculate the director cosines of the ray using two points,
one in each systems. This function allows us to do it.
 
'''
 
def Cosine_transform(system_cos0,system_cos1):
 
    x_0 = system_cos0[0,0]
    y_0 = system_cos0[0,1]
    z_0 = system_cos0[0,2]
    
    x_1 = system_cos1[0,0]
    y_1 = system_cos1[0,1]
    z_1 = system_cos1[0,2]
 
    magnitude_transform = np.sqrt((x_1-x_0)**2+(y_1-y_0)**2+(z_1-z_0)**2)
    
    pm_L = (x_1-x_0)/magnitude_transform
    pm_M = (y_1-y_0)/magnitude_transform
    pm_N = (z_1-z_0)/magnitude_transform
    
    return pm_L, pm_M, pm_N
 

      
##############################################################################
##############################################################################
 
'''
 
n) Ray information through the system
 
All the information up to this point will be stored in two variables, each of which 
will contain ray and surface information.
 
We defined two functions (n.1 and n.2), the first function calls all the previous 
functions that we define and it carries out the exact ray tracing. Moreover, this 
function will store all the information that it generates. The second function 
specifies the displacement and rotation of all the surfaces that will undergo a 
change and calculates the two transformation matrices. Finally, we use the (n.1) 
function to compute the parameters of the rays and the surface.
 
'''
 
##############################################################################
 
'''
n.1 Arrangement function
 
We introduce the arrangement_ray(0, 1.5, 3, 4.0) function which, as we mentioned
in Section c.7), helps us to calculate the direction and the position of an arrangement 
of rays.

'''
 
def Rays_infomation(ray_param,Matrix,dfocus):

  star_pos1, starpos2, cos_pack= arrangement_ray(0, 1.5, 3, 4.0)
  cos_packL = cos_pack[0]
  cos_packM = cos_pack[1]
  cos_packN = cos_pack[2] 
  amount_cos = len(cos_packL)
  systcoorx_5pack = []
  systcoory_5pack = [] 
  systcoorz_5pack = [] 
  systcosL_5pack = [] 
  systcosM_5pack = [] 
  systcosN_5pack = [] 
  systcoorx_4pack = [] 
  systcoory_4pack = []
  systcoorz_4pack = []
  systcosL_4pack = []
  systcosM_4pack = []
  systcosN_4pack = []
  systcoorx_3pack = []
  systcoory_3pack = []
  systcoorz_3pack = []
  systcosL_3pack = []
  systcosM_3pack = []
  systcosN_3pack = []
  systcoorx_2pack = []
  systcoory_2pack = []
  systcoorz_2pack = []
  systcosL_2pack = []
  systcosM_2pack = []
  systcosN_2pack = []
  systcoorx_1pack = []
  systcoory_1pack = []
  systcoorz_1pack = []
  systcosL_1pack = []
  systcosM_1pack = []
  systcosN_1pack = []
  initial_coor0 = ray_param[0]
  initial_coor1 = ray_param[1]
  NR2 = len(initial_coor0)

  for i in range(amount_cos):
      initial_coor2 = np.ones(NR2) * -100
      director_cos0 = np.ones(NR2) * cos_packL[i]
      director_cos1 = np.ones(NR2) * cos_packM[i]
      director_cos2 = np.ones(NR2) * cos_packN[i]
      systcoorx_5 = np.ones(NR2)
      systcoory_5 = np.ones(NR2)
      systcoorz_5 = np.ones(NR2)
      systcosL_5 = np.ones(NR2)
      systcosM_5 = np.ones(NR2)
      systcosN_5 = np.ones(NR2)
      systcoorx_4 = np.ones(NR2)
      systcoory_4 = np.ones(NR2)
      systcoorz_4 = np.ones(NR2)
      systcosL_4 = np.ones(NR2)
      systcosM_4 = np.ones(NR2)
      systcosN_4 = np.ones(NR2) 
      systcoorx_3 = np.ones(NR2)
      systcoory_3 = np.ones(NR2)
      systcoorz_3 = np.ones(NR2)
      systcosL_3 = np.ones(NR2)
      systcosM_3 = np.ones(NR2)
      systcosN_3 = np.ones(NR2) 
      systcoorx_2 = np.ones(NR2)
      systcoory_2 = np.ones(NR2)
      systcoorz_2 = np.ones(NR2)
      systcosL_2 = np.ones(NR2)
      systcosM_2 = np.ones(NR2)
      systcosN_2 = np.ones(NR2) 
      systcoorx_1 = np.ones(NR2)
      systcoory_1 = np.ones(NR2)
      systcoorz_1 = np.ones(NR2)
      systcosL_1 = np.ones(NR2)
      systcosM_1 = np.ones(NR2)
      systcosN_1 = np.ones(NR2)

      ################################################################################
 
      for i in range(NR2):
    
          inx = initial_coor0[i]
          iny = initial_coor1[i]
          inz = initial_coor2[i]
          inL = director_cos0[i]
          inM = director_cos1[i]
          inN = director_cos2[i]
    
          cx_1, cy_1, cz_1, cL, cM, cN = Traceray(inx,iny,inz,inL, inM, inN, Matrix,dfocus)
    
          systcoorx_5[i] = cx_1[4]
          systcoory_5[i] = cy_1[4]
          systcoorz_5[i] = cz_1[4]
          systcosL_5[i] = cL[5]
          systcosM_5[i] = cM[5]
          systcosN_5[i] = cN[5]
          systcoorx_4[i] = cx_1[3]
          systcoory_4[i] = cy_1[3]
          systcoorz_4[i] = cz_1[3]
          systcosL_4[i] = cL[4]
          systcosM_4[i] = cM[4]
          systcosN_4[i] = cN[4]
          systcoorx_3[i] = cx_1[2]
          systcoory_3[i] = cy_1[2]
          systcoorz_3[i] = cz_1[2]
          systcosL_3[i] = cL[3]
          systcosM_3[i] = cM[3]
          systcosN_3[i] = cN[3]
          systcoorx_2[i] = cx_1[1]
          systcoory_2[i] = cy_1[1]
          systcoorz_2[i] = cz_1[1]
          systcosL_2[i] = cL[2]
          systcosM_2[i] = cM[2]
          systcosN_2[i] = cN[2]
          systcoorx_1[i] = cx_1[0]
          systcoory_1[i] = cy_1[0]
          systcoorz_1[i] = cz_1[0]
          systcosL_1[i] = cL[1]
          systcosM_1[i] = cM[1]
          systcosN_1[i] = cN[1]

      systcoorx_5pack.append(systcoorx_5)
      systcoory_5pack.append(systcoory_5) 
      systcoorz_5pack.append(systcoorz_5) 
      systcosL_5pack.append(systcosL_1) 
      systcosM_5pack.append(systcosM_1) 
      systcosN_5pack.append(systcosN_1) 
      systcoorx_4pack.append(systcoorx_4)
      systcoory_4pack.append(systcoory_4) 
      systcoorz_4pack.append(systcoorz_4) 
      systcosL_4pack.append(systcosL_1)
      systcosM_4pack.append(systcosM_1)
      systcosN_4pack.append(systcosN_1)
      systcoorx_3pack.append(systcoorx_3)
      systcoory_3pack.append(systcoory_3)
      systcoorz_3pack.append(systcoorz_3)
      systcosL_3pack.append(systcosL_1)
      systcosM_3pack.append(systcosM_1)
      systcosN_3pack.append(systcosN_1)
      systcoorx_2pack.append(systcoorx_2)
      systcoory_2pack.append(systcoory_2)
      systcoorz_2pack.append(systcoorz_2)
      systcosL_2pack.append(systcosL_1)
      systcosM_2pack.append(systcosM_1)
      systcosN_2pack.append(systcosN_1)
      systcoorx_1pack.append(systcoorx_1)
      systcoory_1pack.append(systcoory_1)
      systcoorz_1pack.append(systcoorz_1)
      systcosL_1pack.append(systcosL_1)
      systcosM_1pack.append(systcosM_1)
      systcosN_1pack.append(systcosN_1)
 
  # #########################################################################
    
  Object_plane = (initial_coor0, initial_coor1, initial_coor2, director_cos0, director_cos1, director_cos2)
  Primary_mirror = (systcoorx_1pack, systcoory_1pack, systcoorz_1pack, systcosL_1pack, systcosM_1pack, systcosN_1pack)
  Secondary_mirror = (systcoorx_2pack, systcoory_2pack, systcoorz_2pack, systcosL_2pack, systcosM_2pack, systcosN_2pack)
  Correctingasphlens = (systcoorx_3pack, systcoory_3pack, systcoorz_3pack, systcosL_3pack, systcosM_3pack, systcosN_3pack)
  Corrector_plane = (systcoorx_4pack, systcoory_4pack, systcoorz_4pack, systcosL_4pack, systcosM_4pack, systcosN_4pack)
  Image_plane = (systcoorx_5pack, systcoory_5pack, systcoorz_5pack, systcosL_5pack, systcosM_5, systcosN_5pack)
  surf_opt = (Primary_mirror, Secondary_mirror, Correctingasphlens)
  plane_opt = (Object_plane, Corrector_plane, Image_plane)
 
  return surf_opt, plane_opt
 
 
###############################################################################
 
'''
 
n.2 function
 
'''
 
def TiltandShift_transtormation(tilt_smirror, shift_smirror, set_rayos, dfocus):
  tiltthetax = tilt_smirror[0]
  tiltthetay = tilt_smirror[1]
  shiftx = shift_smirror[0]
  shifty = shift_smirror[1]
  MT=[]
  ma_rt_s0 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s1 = [tiltthetax , tiltthetay, 0.0, shiftx, shifty, 0.0]
  ma_rt_s2 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s3 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s4 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s5 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s6 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  ma_rt_s7 = [0.0 , 0.0, 0.0, 0.0, 0.0, 0.0]
  M = Transform(ma_rt_s0)
  MT.append(M)
  M = Transform(ma_rt_s1)
  MT.append(M)
  M = Transform(ma_rt_s2)
  MT.append(M)
  M = Transform(ma_rt_s3)
  MT.append(M)
  M = Transform(ma_rt_s4)
  MT.append(M)
  surf_transformate, plane_transformate = Rays_infomation(set_rayos,MT, dfocus)
  return surf_transformate, plane_transformate
 
 
###############################################################################
###############################################################################


surf_set , ray_set = setof_rayos1(0)

nonshift_smirror = (0.0, 0.0)
nontilt_smirror = (0.0, 0.0)

shift_smirror = (0.0, 1.0)
tilt_smirror = (0.0, 0.0)

shift_smirror = (0.0, 1.0)
tiltcompensated_smirror = (0.10153185221278913, 0.0)

displacement_imageplane = 0.0
ray_surf, ray_plane = TiltandShift_transtormation(nontilt_smirror,nonshift_smirror, ray_set, displacement_imageplane)

displacement_imageplane = 0.0
ray_shiftedsurf, ray_shiftedplane = TiltandShift_transtormation(tilt_smirror,shift_smirror, ray_set,displacement_imageplane)

displacement_imageplane = 0.0
ray_correctedsurf, ray_correctedplane = TiltandShift_transtormation(tiltcompensated_smirror ,shift_smirror, ray_set, displacement_imageplane)

displacement_imageplane = 0.4
ray_correcteddonutsurf, ray_correcteddonutplane = TiltandShift_transtormation(tiltcompensated_smirror ,shift_smirror, ray_set, displacement_imageplane)

###############################################################################

'''
Part IV.
Results
 
In this last Section, nine fields of the focal plane are shown, divided into four cases: 
the first shows the optic system perfectly aligned, the second shows the optic system 
shifted, the third plot displays a compensated out-of-focus optical system, and the 
last plot shows the compensated system. These four plots set out the nine fields. 

'''
print('')
print('')
###############################################################################

'''
IV.a) Optic system perfectly aligned
'''

'''
coordinate x
'''

Plane_imagenninefieldx = ray_plane[2][0]
Plane_imagenninefieldy = ray_plane[2][1]
coordinatex_star1 = Plane_imagenninefieldx [6]
coordinatex_star2 = Plane_imagenninefieldx [7]
coordinatex_star3 = Plane_imagenninefieldx [8]
coordinatex_star4 = Plane_imagenninefieldx [3]
coordinatex_star5 = Plane_imagenninefieldx [4]
coordinatex_star6 = Plane_imagenninefieldx [5]
coordinatex_star7 = Plane_imagenninefieldx [0]
coordinatex_star8 = Plane_imagenninefieldx [1]
coordinatex_star9 = Plane_imagenninefieldx [2]

'''
coordinate y
'''

coordinatey_star1 = Plane_imagenninefieldy [6]
coordinatey_star2 = Plane_imagenninefieldy [7]
coordinatey_star3 = Plane_imagenninefieldy [8]
coordinatey_star4 = Plane_imagenninefieldy [3]
coordinatey_star5 = Plane_imagenninefieldy [4]
coordinatey_star6 = Plane_imagenninefieldy [5]
coordinatey_star7 = Plane_imagenninefieldy [0]
coordinatey_star8 = Plane_imagenninefieldy [1]
coordinatey_star9 = Plane_imagenninefieldy [2]

###############################################################################
'''
Plot of the optic system perfectly aligned
'''

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(nrows=3,ncols=3,figsize=(8,8),sharex='col' ,sharey='row')
ax1.plot(coordinatex_star1, coordinatey_star1,'.', color='k', linestyle='none')
ax2.plot(coordinatex_star2, coordinatey_star2,'.', color='k', linestyle='none')
ax3.plot(coordinatex_star3, coordinatey_star3,'.', color='k', linestyle='none')
ax4.plot(coordinatex_star4, coordinatey_star4,'.', color='k', linestyle='none')
ax5.plot(coordinatex_star5, coordinatey_star5,'.', color='k', linestyle='none')
ax6.plot(coordinatex_star6, coordinatey_star6,'.', color='k', linestyle='none')
ax7.plot(coordinatex_star7, coordinatey_star7,'.', color='k', linestyle='none')
ax8.plot(coordinatex_star8, coordinatey_star8,'.', color='k', linestyle='none')
ax9.plot(coordinatex_star9, coordinatey_star9,'.', color='k', linestyle='none')

ax1.set_yticks([17.96])
ax1.set_yticklabels(["1.5"])
ax1.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax2.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax3.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax4.set_yticks([0.0])
ax4.set_yticklabels(["0.0"])
ax4.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax5.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax6.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax7.set_yticks([-17.96])
ax7.set_yticklabels(["-1.5"])
majors7 = [ -17.96]
ax7.xaxis.set_major_locator(ticker.FixedLocator(majors7))
ax7.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors8 = [ 0.0]
ax8.xaxis.set_major_locator(ticker.FixedLocator(majors8))
ax8.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors9 = [ 17.96]
ax9.xaxis.set_major_locator(ticker.FixedLocator(majors9))
ax9.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')

ax8.set_xlabel(r'Distance from paraxial focal plane [mm]', fontsize=22, fontstyle="normal", family='serif')
ax4.set_ylabel(r'Field angle [degrees]', fontsize=22, fontstyle="normal", family='serif')
fig.tight_layout()
fig.subplots_adjust(hspace=0.00, wspace=0.00, top=0.93)
plt.figure(1) 
plt.show()

print('')
print('')

###############################################################################
###############################################################################

'''
IV.b) Optic system shifted
'''

'''
coordinate x
'''

Plane_imagenninefieldshiftedx = ray_shiftedplane[2][0]
Plane_imagenninefieldshiftedy = ray_shiftedplane[2][1]
coordinateshiftedx_star1 = Plane_imagenninefieldshiftedx [6]
coordinateshiftedx_star2 = Plane_imagenninefieldshiftedx [7]
coordinateshiftedx_star3 = Plane_imagenninefieldshiftedx [8]
coordinateshiftedx_star4 = Plane_imagenninefieldshiftedx [3]
coordinateshiftedx_star5 = Plane_imagenninefieldshiftedx [4]
coordinateshiftedx_star6 = Plane_imagenninefieldshiftedx [5]
coordinateshiftedx_star7 = Plane_imagenninefieldshiftedx [0]
coordinateshiftedx_star8 = Plane_imagenninefieldshiftedx [1]
coordinateshiftedx_star9 = Plane_imagenninefieldshiftedx [2]

'''
coordinate y
'''

coordinateshiftedy_star1 = Plane_imagenninefieldshiftedy [6]
coordinateshiftedy_star2 = Plane_imagenninefieldshiftedy [7]
coordinateshiftedy_star3 = Plane_imagenninefieldshiftedy [8]
coordinateshiftedy_star4 = Plane_imagenninefieldshiftedy [3]
coordinateshiftedy_star5 = Plane_imagenninefieldshiftedy [4]
coordinateshiftedy_star6 = Plane_imagenninefieldshiftedy [5]
coordinateshiftedy_star7 = Plane_imagenninefieldshiftedy [0]
coordinateshiftedy_star8 = Plane_imagenninefieldshiftedy [1]
coordinateshiftedy_star9 = Plane_imagenninefieldshiftedy [2]

################################################################################
'''
Plot of the optic system shifted
'''

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(nrows=3,ncols=3,figsize=(8,8),sharex='col' ,sharey='row')
ax1.plot(coordinateshiftedx_star1, coordinateshiftedy_star1,'.', color='k', linestyle='none')
ax2.plot(coordinateshiftedx_star2, coordinateshiftedy_star2,'.', color='k', linestyle='none')
ax3.plot(coordinateshiftedx_star3, coordinateshiftedy_star3,'.', color='k', linestyle='none')
ax4.plot(coordinateshiftedx_star4, coordinateshiftedy_star4,'.', color='k', linestyle='none')
ax5.plot(coordinateshiftedx_star5, coordinateshiftedy_star5,'.', color='k', linestyle='none')
ax6.plot(coordinateshiftedx_star6, coordinateshiftedy_star6,'.', color='k', linestyle='none')
ax7.plot(coordinateshiftedx_star7, coordinateshiftedy_star7,'.', color='k', linestyle='none')
ax8.plot(coordinateshiftedx_star8, coordinateshiftedy_star8,'.', color='k', linestyle='none')
ax9.plot(coordinateshiftedx_star9, coordinateshiftedy_star9,'.', color='k', linestyle='none')

ax1.set_yticks([17.36])
ax1.set_yticklabels(["1.5"])
ax1.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax2.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax3.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax4.set_yticks([-0.6])
ax4.set_yticklabels(["0.0"])
ax4.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax5.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax6.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax7.set_yticks([-18.560])
ax7.set_yticklabels(["-1.5"])
majors7 = [ -17.96]
ax7.xaxis.set_major_locator(ticker.FixedLocator(majors7))
ax7.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors8 = [ 0.0]
ax8.xaxis.set_major_locator(ticker.FixedLocator(majors8))
ax8.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors9 = [ 17.96]
ax9.xaxis.set_major_locator(ticker.FixedLocator(majors9))
ax9.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')

ax8.set_xlabel(r'Distance from paraxial focal plane [mm]', fontsize=22, fontstyle="normal", family='serif')
ax4.set_ylabel(r'Field angle [degrees]', fontsize=22, fontstyle="normal", family='serif')
fig.tight_layout()
fig.subplots_adjust(hspace=0.00, wspace=0.00, top=0.93)
plt.figure(2) 
plt.show()

print('')
print('')

###############################################################################
###############################################################################
'''
IV.c) Compensated out-of-focus optical system
'''

'''
coordinate x
'''

Plane_imagedonutx = ray_correcteddonutplane[2][0]
coordinatecompensateddonutx_star1 = Plane_imagedonutx [6]
coordinatecompensateddonutx_star2 = Plane_imagedonutx [7]
coordinatecompensateddonutx_star3 = Plane_imagedonutx [8]
coordinatecompensateddonutx_star4 = Plane_imagedonutx [3]
coordinatecompensateddonutx_star5 = Plane_imagedonutx [4]
coordinatecompensateddonutx_star6 = Plane_imagedonutx [5]
coordinatecompensateddonutx_star7 = Plane_imagedonutx [0]
coordinatecompensateddonutx_star8 = Plane_imagedonutx [1]
coordinatecompensateddonutx_star9 = Plane_imagedonutx [2]

'''
coordinate y
'''

Plane_imagedonuty = ray_correcteddonutplane[2][1]
coordinatecompensateddonuty_star1 = Plane_imagedonuty [6]
coordinatecompensateddonuty_star2 = Plane_imagedonuty [7]
coordinatecompensateddonuty_star3 = Plane_imagedonuty [8]
coordinatecompensateddonuty_star4 = Plane_imagedonuty [3]
coordinatecompensateddonuty_star5 = Plane_imagedonuty [4]
coordinatecompensateddonuty_star6 = Plane_imagedonuty [5]
coordinatecompensateddonuty_star7 = Plane_imagedonuty [0]
coordinatecompensateddonuty_star8 = Plane_imagedonuty [1]
coordinatecompensateddonuty_star9 = Plane_imagedonuty [2]

###############################################################################
###############################################################################

'''
Plot of the compensated out-of-focus optical system
'''

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(nrows=3,ncols=3,figsize=(8,8),sharex='col' ,sharey='row')
ax1.plot(coordinatecompensateddonutx_star1, coordinatecompensateddonuty_star1,'.', color='k', linestyle='none')
ax2.plot(coordinatecompensateddonutx_star2, coordinatecompensateddonuty_star2,'.', color='k', linestyle='none')
ax3.plot(coordinatecompensateddonutx_star3, coordinatecompensateddonuty_star3,'.', color='k', linestyle='none')
ax4.plot(coordinatecompensateddonutx_star4, coordinatecompensateddonuty_star4,'.', color='k', linestyle='none')
ax5.plot(coordinatecompensateddonutx_star5, coordinatecompensateddonuty_star5,'.', color='k', linestyle='none')
ax6.plot(coordinatecompensateddonutx_star6, coordinatecompensateddonuty_star6,'.', color='k', linestyle='none')
ax7.plot(coordinatecompensateddonutx_star7, coordinatecompensateddonuty_star7,'.', color='k', linestyle='none')
ax8.plot(coordinatecompensateddonutx_star8, coordinatecompensateddonuty_star8,'.', color='k', linestyle='none')
ax9.plot(coordinatecompensateddonutx_star9, coordinatecompensateddonuty_star9,'.', color='k', linestyle='none')

ax1.set_yticks([26.53])
ax1.set_yticklabels(["1.5"])
ax1.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax2.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax3.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax4.set_yticks([8.57])
ax4.set_yticklabels(["0.0"])
ax4.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax5.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax6.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax7.set_yticks([-9.40])
ax7.set_yticklabels(["-1.5"])
majors7 = [ -17.96]
ax7.xaxis.set_major_locator(ticker.FixedLocator(majors7))
ax7.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors8 = [ 0.0]
ax8.xaxis.set_major_locator(ticker.FixedLocator(majors8))
ax8.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors9 = [ 17.96]
ax9.xaxis.set_major_locator(ticker.FixedLocator(majors9))
ax9.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')

ax8.set_xlabel(r'Distance from paraxial focal plane [mm]', fontsize=22, fontstyle="normal", family='serif')
ax4.set_ylabel(r'Field angle [degrees]', fontsize=22, fontstyle="normal", family='serif')
fig.tight_layout()
fig.subplots_adjust(hspace=0.00, wspace=0.00, top=0.93)
plt.figure(3)
plt.show()


###############################################################################
###############################################################################

'''
IV.d) Compensated optic system 
'''

'''
coordinate x
'''

Plane_imageninefieldsx = ray_correctedplane[2][0]
Plane_imageninefieldsy = ray_correctedplane[2][1]
coordinatecompensatedx_star1 = Plane_imageninefieldsx [6]
coordinatecompensatedx_star2 = Plane_imageninefieldsx [7]
coordinatecompensatedx_star3 = Plane_imageninefieldsx [8]
coordinatecompensatedx_star4 = Plane_imageninefieldsx [3]
coordinatecompensatedx_star5 = Plane_imageninefieldsx [4]
coordinatecompensatedx_star6 = Plane_imageninefieldsx [5]
coordinatecompensatedx_star7 = Plane_imageninefieldsx [0]
coordinatecompensatedx_star8 = Plane_imageninefieldsx [1]
coordinatecompensatedx_star9 = Plane_imageninefieldsx [2]

'''
coordinate y
'''

coordinatecompensatedy_star1 = Plane_imageninefieldsy  [6]
coordinatecompensatedy_star2 = Plane_imageninefieldsy  [7]
coordinatecompensatedy_star3 = Plane_imageninefieldsy  [8]
coordinatecompensatedy_star4 = Plane_imageninefieldsy  [3]
coordinatecompensatedy_star5 = Plane_imageninefieldsy  [4]
coordinatecompensatedy_star6 = Plane_imageninefieldsy  [5]
coordinatecompensatedy_star7 = Plane_imageninefieldsy  [0]
coordinatecompensatedy_star8 = Plane_imageninefieldsy  [1]
coordinatecompensatedy_star9 = Plane_imageninefieldsy  [2]

###############################################################################
'''
Plot of the compensated optical system
'''

fig, ((ax1,ax2,ax3),(ax4,ax5,ax6),(ax7,ax8,ax9)) = plt.subplots(nrows=3,ncols=3,figsize=(8,8),sharex='col' ,sharey='row')
ax1.plot(coordinatecompensatedx_star1, coordinatecompensatedy_star1,'.', color='k', linestyle='none')
ax2.plot(coordinatecompensatedx_star2, coordinatecompensatedy_star2,'.', color='k', linestyle='none')
ax3.plot(coordinatecompensatedx_star3, coordinatecompensatedy_star3,'.', color='k', linestyle='none')
ax4.plot(coordinatecompensatedx_star4, coordinatecompensatedy_star4,'.', color='k', linestyle='none')
ax5.plot(coordinatecompensatedx_star5, coordinatecompensatedy_star5,'.', color='k', linestyle='none')
ax6.plot(coordinatecompensatedx_star6, coordinatecompensatedy_star6,'.', color='k', linestyle='none')
ax7.plot(coordinatecompensatedx_star7, coordinatecompensatedy_star7,'.', color='k', linestyle='none')
ax8.plot(coordinatecompensatedx_star8, coordinatecompensatedy_star8,'.', color='k', linestyle='none')
ax9.plot(coordinatecompensatedx_star9, coordinatecompensatedy_star9,'.', color='k', linestyle='none')

ax1.set_yticks([26.525])
ax1.set_yticklabels(["1.5"])
ax1.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax2.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax3.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle = '--')
ax4.set_yticks([8.565])
ax4.set_yticklabels(["0.0"])
ax4.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax5.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax6.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
ax7.set_yticks([-9.395])
ax7.set_yticklabels(["-1.5"])
majors7 = [ -17.965]
ax7.xaxis.set_major_locator(ticker.FixedLocator(majors7))
ax7.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors8 = [ 0.0]
ax8.xaxis.set_major_locator(ticker.FixedLocator(majors8))
ax8.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')
majors9 = [ 17.965]
ax9.xaxis.set_major_locator(ticker.FixedLocator(majors9))
ax9.tick_params(axis='both', width=2, length=5, labelsize=14, gridOn= True, grid_linestyle ='--')

ax8.set_xlabel(r'Distance from paraxial focal plane [mm]', fontsize=22, fontstyle="normal", family='serif')
ax4.set_ylabel(r'Field angle [degrees]', fontsize=22, fontstyle="normal", family='serif')
fig.tight_layout()
fig.subplots_adjust(hspace=0.00, wspace=0.00, top=0.93)
plt.figure(4) 
plt.show()

