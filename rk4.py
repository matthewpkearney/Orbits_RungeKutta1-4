from math import *
import matplotlib.pyplot as plt

""" 
RK2 
New variables (k) are used to iterate more precisely through time steps
"""

##_T represents the time to calculate until
_T = 400*24*60*60 #86400[s] One week
#dt represents time steps  (delta time)
dt= 86400 #60[s]
##thus the number of time steps is 5760

##AU is 1 astronomical unit
AU = 1.495978707*pow(10, 11) #[m] distance from COM(Earth) to COM(Sun)
univ_grav= float(6.67*pow(10,-11))  #G = 6.67*10^-11[N*m^2/kg^2]

#variables m_sun and m_earth can represent any objects in the universe

m_sun=float(1.9891*pow(10,30))   #[kg] 
m_earth=float(5.97219*pow(10,24))#[kg]

##initial conditions, r=(x,y) and v=(vx,vy)
r=[AU, 0]     #x=1AU, y=0
v=[0, 29800]  #vx=0,   vy=29800 (launch vel)


""" 
    Procedure for Runge Kutta 2 goes as:
    inital conditions -> k1 -> k2 -> updated value (x,y) & (vx, vy)
"""

""" Step One: Calculate the first set of k's from initial condiitons
(ri, vi) -> k (represent as a tuple to get kx/ky values)

k1x = vix*t 
k1y = viy*t

k1vx = (-G*M*x*t/((x^2 + y^2)^(3/2)))
k1vy = (-G*M*y*t/((x^2 + y^2)^(3/2))) """




#calculating next_pos (in time step) using position, velocity, and time
def k1r(v, dt):
    '''calculate k1 position'''
    xf1= v[0]*dt
    yf1= v[1]*dt 
    return [xf1, yf1]

#calculating next_vel (in time step) using position, velocity, and time
def k1v(r, dt):
    '''calculate k1 velocity'''
    vx1 = (-univ_grav*m_sun*r[0]*dt)/(pow ( (pow (r[0], 2) + pow (r[1], 2)), 1.5) )
    vy1 = (-univ_grav*m_sun*r[1]*dt)/(pow ( (pow (r[0], 2) + pow (r[1], 2)), 1.5) )
    return [vx1, vy1]



"""Step Two: Calculate the next set of k's from initial conditions
k1 -> k2 


k2x = (vix + k1x)*t
k2y = (xiy + k1y)*t

k2vx = (-G*M*(x + k1x)*t) / ((x+k1x)^2 + (y + k1y)^2)^3/2)
k2vy = (-G*M*(y + k1y)*t) / ((x+k1x)^2 + (y + k1y)^2)^3/2)
"""


def k2r(r, v, dt):
    """calculate position in Runge Kutta 2nd order, need velocity and t"""
    k1vel = k1v(r, dt)
    k2x = (v[0] + k1vel[0]/2)*dt
    k2y = (v[1] + k1vel[1]/2)*dt
    return [k2x, k2y]

def k2v(r, dt):
    """error here, should be r[1] + k1s[1]..."""
    k1s = k1r(r, dt)
    k2vx = -((univ_grav*m_sun*(r[0] + k1s[0]/2)*dt)/(pow (( pow  (r[0]+ k1s[0]/2 , 2)  + pow (r[1]+k1s[1]/2 , 2 )) , 1.5  )))
    k2vy = -((univ_grav*m_sun*(r[1] + k1s[1]/2)*dt)/(pow (( pow  (r[0]+ k1s[0]/2 , 2)  +  pow (r[1]+k1s[1]/2 , 2 )) , 1.5  )))
    return [k2vx, k2vy]



'''Step Three: 

k3x=( vx + k1vx)*dt
k3y=( vy + k1vy)*dt


k3vx=-G*M*x*dt/((x+k2x)^2+(y+k2y)^2)^3/2
k3vy=-G*M*y*dt/((x+k2x)^2+(y+k2y)^2)^3/2
'''

def k3r(r, v, dt):
    '''return pos in k3'''
    k2vel = k2v(r, dt)
    k3x = (v[0] + k2vel[0]*0.5)*dt
    k3y = (v[1] + k2vel[1]*0.5)*dt
    return [k3x, k3y]

def k3v(r, v, dt):
    '''return velocity in k3'''
    k2s = k2r(r, v, dt)
    k3vx = -(univ_grav*m_sun*(r[0] + k2s[0])*0.5) / (pow (    (pow  ( r[0]+ k2s[0]*0.5 , 2)  +  pow (r[1]+k2s[1]*0.5 , 2  )) , 1.5  ))*dt
    k3vy = -(univ_grav*m_sun*(r[1] + k2s[1])*0.5) / (pow (    (pow  ( r[0]+ k2s[0]*0.5 , 2)  +  pow (r[1]+k2s[1]*0.5 , 2  )) , 1.5  ))*dt
    return [k3vx, k3vy]



'''Step Four: 

k4x = (vx + k3vx)*dt
k4y = (vy + k3vy)*dt

k4vx = -GM(x+k3x)*dt/((x+k3x)^2 + (y + k3y)^2)^3/2
k4vy = -GM(y+k3y)*dt/((x+k3x)^2 + (y + k3y)^2)^3/2'''

def k4r(r, v, dt):
    '''Calculate pos in 4th order'''
    k3vel = k3v(r, v, dt)
    k3x = (v[0] + k3vel[0])*dt
    k3y = (v[1] + k3vel[1])*dt
    return [k3x, k3y]

def k4v(r, v, dt):
    '''calculate velocity in k4'''
    k3s = k3r(r, v, dt)
    k4vx = -(univ_grav*m_sun*(r[0] + k3s[0])*dt)/(pow ( ( pow  ( r[0]+ k3s[0], 2)  +  pow (r[1]+k3s[1] , 2 )) , 1.5  ))
    k4vy = -(univ_grav*m_sun*(r[1] + k3s[1])*dt)/(pow ( ( pow  ( r[0]+ k3s[0], 2)  +  pow (r[1]+k3s[1] , 2 )) , 1.5  ))
    return [k4vx, k4vy]



"""Step Three: Calculate new variables
x= x + (1/2)*(k1x + k2x)
y = y + (1/2)*(k1y + k2y)
vx = vx + (1/2)*(k1vx + k2vx)
vy = vy + (1/2)*k1vx + k2vx
"""

def update_pos(r, v, dt):
    """calculate updated position"""
    k1s= k1r(v, dt)
    k2s = k2r(r, v,  dt)
    k3s = k3r(r, v, dt)
    k4s = k4r(r, v, dt)
    # print("k1s : " + str(k1s))
    # print("k2s : " + str(k2s))
    r[0] = r[0] + (1/6)*(k1s[0] + 2*k2s[0] + 2*k3s[0] + k4s[0])
    r[1] = r[1] + (1/6)*(k1s[1] + 2*k2s[1] + 2*k3s[1] + k4s[1])
    return [r[0], r[1]]

def update_vel(r, v, dt):
    """calculate updated velocity"""
    k1vel = k1v(r, dt)
    k2vel = k2v(r, dt)
    k3vel = k3v(r, v, dt)
    k4vel = k4v(r, v, dt)
    v[0] = v[0] + (1/6)*(k1vel[0] + 2*k2vel[0] + 2*k3vel[0] + k4vel[0])
    v[1] = v[1] + (1/6)*(k1vel[1] + 2*k2vel[1] + 2*k3vel[1] + k4vel[1])
    return [v[0], v[1]]

##store the x,y coordinates
pos_coordinates = open("351coordinates.txt", "w")
#graph shows distance of about 1.5 AU's around the sun
ax = plt.axes(xlim=(-30*AU, 10*AU), ylim =(-5*AU, 10*AU))
ax.set_facecolor('black') #space <- black
plt.title('Earths Orbit Around the Sun - RK4')
plt.xlabel('x-pos')
plt.ylabel('y-pos')

#store x and y values 
x = []
y = []

"""Error here: not iterating through time correctly, 
should be using delta time as a multiplier (small t)
but am using the current t"""

#start at time 0 
ti=0   
#obtain all x,y values while (ti<!T)
while(ti<=_T):
    #write (x,y) coordinates to file for matplotlib starting here
    x.append(r[0])
    y.append(r[1])
    r = update_pos(r, v, ti) #calculate next conditions
    v = update_vel(r, v, ti)
    ti = ti + dt #adding time step

plt.plot(x,y,"cyan") #plot all x,y in greenish
plt.scatter( 0 , 0 , s = 20 ).set_facecolor('orange') #sun

plt.show() #show




