from math import *
import matplotlib.pyplot as plt


##_T represents the time to calculate until
_T = 8*24*60*60 #[s] 
#dt represents time steps  (delta time)
dt= 3600 #15[s] 
##thus the number of time steps is 5760

##AU is 1 astronomical unit
AU = 1.495978707*pow(10, 11) #[m] distance from COM(Earth) to COM(Sun)
univ_grav= float(6.67408*pow(10,-11))  #G = 6.67*10^-11[N*m^2/kg^2]

#variables m_sun and m_earth can represent any objects in the universe

m_sun=float(1.9891*pow(10,30))   #[kg] 
m_earth=float(5.97219*pow(10,24))#[kg]


###Helper Functions for i step
def unit_vect( x ):
        '''Given array x (w/ 2 elements), return the value of x as a unit vector [rx, ry]'''
        return [x[0]/dist(x), x[1]/dist(x)]

##calculating ri -> distance sun to earth
def dist( x ):
        '''Return the distance c, using the Pythagorean Theorem
        pos is parameter containing pos[0] = x and pos[1] as y'''
        return sqrt(pow(x[0], 2) + pow(x[1], 2))

    
##initial conditions
r=[AU, 0]     #x=1AU, y=0

v=[0, 32700]  #vx=0,   vy=29800 (launch vel)


##calculating acceleration given position
def calc_accel(r):
    '''return the acceleration as a vector (list)
    acceleration is force / mass, uses force function to find angular acceleration'''
    ax= -(univ_grav*m_sun*r[0])/(pow(dist(r), 3))
    ay= -(univ_grav*m_sun*r[1])/(pow(dist(r), 3))
    return [ax, ay]

#calculating next_pos (in time step) using position, velocity, and time
def next_pos(r, v,  dt):
    '''calculate the next position using kinematics'''
    acc=calc_accel(r) #bring in previous acceleration (ai)
    xf=r[0] + v[0]*dt + 0.5*acc[0]*pow(dt, 2)
    yf=r[1] + v[1]*dt + 0.5*acc[1]*pow(dt, 2)
    return [xf, yf]

#calculating next_vel (in time step) using position, velocity, and time
def next_vel(r, v, dt):
    '''calculates the next velocity using kinematic equations'''
    acc=calc_accel(r)
    vxf = v[0] + acc[0]*dt
    vyf = v[1] + acc[1]*dt
    return [vxf, vyf]


##store the x,y coordinates
pos_coordinates = open("351coordinates.txt", "w")
#graph shows distance of about 1.5 AU's around the sun
ax = plt.axes(xlim=(-2*AU, 2*AU), ylim =(-2*AU, 2*AU))
ax.set_facecolor('black') #space <- black
plt.title('Earths Orbit Around the Sun - Eulers Method')
plt.xlabel('x-pos')
plt.ylabel('y-pos')


#plt.grid()  #-> to show the gridding
#plt.axhline().set_color('indigo')  #horizontal
#plt.axvline().set_color('indigo') #vertical


#store x and y values 
x = []
y = []

#start at time 0 
ti=0   
#obtain all x,y values while (ti<!T)
while(ti<=_T):
    #write (x,y) coordinates to file for matplotlib starting here
    r2=next_pos(r, v, ti) #calculate next pos
    v2=next_vel(r, v, ti) #calculate next vel
    r= r2 #changing to next pos
    v=v2  #changing to next velo
    x.append(r[0])
    y.append(r[1])
    ti = ti + dt #adding time step


plt.plot(x,y,"cyan") #plot all x,y in greenish
plt.scatter( 0 , 0 , s = 20 ).set_facecolor('orange')

plt.show() #show