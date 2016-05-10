# Diffusion
# from myIterativeMethods import explicit_euler, implicit_euler, ...

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import time
#import scipy.linalg
#from scipy.special import erf
from heapq import heappush, heappop, heapify

# --- Import my functions.
#from tdma import tdma

# --- Class.
class Particle():

    particle_count = 0

    def __init__(self):
        self.x, self.y   = 0.5, 0.5
        self.vx, self.vy = 0.5, 0.5
        self.radius      = 1.e-2
        self.mass        = 1.0
        self.colcount    = 0
        self.id          = Particle.particle_count
        Particle.particle_count += 1

    def __init__(self, x, y, vx, vy, r, m):
        self.x, self.y   = x, y
        self.vx, self.vy = vx, vy
        self.radius      = r
        self.mass        = m
        self.colcount    = 0
        self.id          = Particle.particle_count
        Particle.particle_count += 1


    # --- Gettere.
    def position(self):
        return np.array([self.x, self.y])

    def velocity(self):
        return np.array([self.vx, self.vy])

    def radius(self):
        return self.radius

    def mass(self):
        return self.mass

    def increment_colcount(self):
        self.colcount += 1

    # --- Settere.
    def set_velocity(self, v):
        self.vx = v[0]
        self.vy = v[1]

    # --- Functions.
    def update(self, dt):
        self.x = self.x + dt*self.vx
        self.y = self.y + dt*self.vy

# --- Class: Collision
class Collision():
    def __init__(self, time, involved, colcount, new_velocities):  
        self.time = time
        self.involved = involved
        self.colcount = colcount
        self.new_velocities = new_velocities

    def __cmp__(self, other):
        return -1.0*(self.time < other.time) + 1.0*(self.time >= other.time)

    #Collision(temp, np.array([i, j]), np.array([particles[i].colcount, particles[j].colcount]), np.array([v1, v2]))      
   
# --- Function: Calculate time until collision with a wall.
#     input : Particle
#     output: time until collision with wall, corresponding new velocity after collision.
def time_until_wall_collision(particle):

    # --- Calculate time until collision with vertical wall.
    if particle.vx > 0.0:
        timex = (1 - particle.radius - particle.x)/particle.vx
    elif particle.vx < 0.0:
        timex = (particle.radius - particle.x)/particle.vx
    else:
        timex = float('inf')

    # --- Calculate time until collision with horizontal wall.
    if particle.vy > 0.0:
        timey = (1 - particle.radius - particle.y)/particle.vy
    elif particle.vy < 0.0:
        timey = (particle.radius - particle.y)/particle.vy
    else:
        timey = float('inf')

    # --- Return the smallest time with the corresponding new velocity.
    if timex < timey:
        return timex, np.array([-xi*particle.vx, xi*particle.vy])
    else:
        return timey, np.array([xi*particle.vx, -xi*particle.vy])

# --- Calculate time until collision between particles.
#     input : two Particles
#     output: time until collision between particles, corresponding new velocity after collision.
def time_until_collision(particle1, particle2):
    dx = particle2.position() - particle1.position()
    dv = particle2.velocity() - particle1.velocity()
    R  = particle1.radius + particle2.radius
    d  = (np.dot(dv,dx))**2 - (np.dot(dv,dv)) * ((np.dot(dx,dx)) - R**2)
    if np.dot(dv,dx) >= 0:
        return float('inf'), -1, -1
    elif d <= 0:
        return float('inf'), -1, -1
    else:
        dt = -(np.dot(dv,dx) + np.sqrt(d))/(np.dot(dv,dv))
        dx2 = particle2.position() + dt*particle2.velocity() - (particle1.position()+dt*particle1.velocity())
        velocity1 = particle1.velocity() + ((1.0+xi)*particle2.mass/(particle1.mass+particle2.mass)*np.dot(dv,dx2)/R**2)*dx2
        velocity2 = particle2.velocity() - ((1.0+xi)*particle1.mass/(particle1.mass+particle2.mass)*np.dot(dv,dx2)/R**2)*dx2
        return dt, velocity1, velocity2

# --- Function: Calculate if a new particle overlaps one of 
#     the existing particles.
#     input : Particle, list of Particles.
#     output: boolean.
def is_overlap(particle, particles):
    for p in particles:
        distance = np.sqrt((particle.x - p.x)**2 + (particle.y - p.y)**2)
        if distance < particle.radius + p.radius:
            return True
    return False

# --- Function: Plot the distribution of speeds in a histogram.
#     input : list of speeds, filename of savefile.
#     output: none, but saves a figure of the histogram.
def plot_speed_distribution(speeds, filename):

    # --- Divide data into Nbins.
    Nbins = 101
    #minimum = -1.0; maximum = 3.0
    minimum = -1.0; maximum = 5.0
    hist, bins = np.histogram(speeds, bins=Nbins, range=(minimum, maximum), normed=True)
    bin_values = (bins[:-1] + bins[1:]) / 2

    # --- Plot histogram.
    histogram = plt.figure()
    #plt.scatter(bin_values, hist, s=50, facecolors='none', edgecolors='r')
    plt.bar(bin_values, hist, width=bins[1]-bins[0], color='blue')

    plt.xlabel(r'Speeds $v$ [m/s]', fontsize=20)
    plt.ylabel(r'$P(v)$', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.axis([minimum, maximum, 0.0, 1.1*max(hist)])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig(filename)

# --- Function: Calculate the average of the variable 
#     colcount in Particle.
#     input : list of Particles.
#     output: average colcount.
def calculate_avg_colcount(particles):
    avg = 0.0
    for p in particles:
        avg += p.colcount
    return avg/len(particles)

# --- Function: Evaluates if the colcount hasn't 
#     changed since the collision was calculated.
#     input : indexes of the involved particles, list of all particles, old colcount.
#     output: boolean.
def is_colcountOK(involved, particles, colcount):
    for i in range(len(involved)):
        if colcount[i] != particles[involved[i]].colcount: 
            return False
    return True

# --- Function: Get the corresponding vector 
#     of speeds from a list of Particles.
#     input : list of Particles
#     output: list of speeds.
def get_speed_vector(particles):
    Nparticles = len(particles)
    speeds = np.zeros(Nparticles)
    for i in range(Nparticles):
        speeds[i] = np.sqrt(particles[i].vx**2 + particles[i].vy**2)
    return speeds

def run(particles, start_time, stop_colcount, plotter=False):

    # --- Initialise square box.
    '''xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0'''

    # --- Set start time.
    t0 = start_time

    #print('Start1: ', particles[0].position(), particles[0].velocity(), particles[0].radius, particles[0].mass)
    #print('Start2: ', particles[1].position(), particles[1].velocity(), particles[1].radius, particles[1].mass)
    
    # --- Initialise priority queue.
    collisions = []
    for i in range (len(particles)):

        # Calculate if and when particle will collide with all other particles,
        # and store the collision times.
        for j in range(i+1, len(particles)):
            #print(i, j)
            dt, v1, v2 = time_until_collision(particles[i], particles[j])
            if dt < float('inf'):
                temp = t0 + dt
                heappush(collisions, Collision(temp, np.array([i, j]), np.array([particles[i].colcount, particles[j].colcount]), np.array([v1, v2])))

        # Calculate when particle will collide with a wall,
        # and store the collision times.
        dt, v1 = time_until_wall_collision(particles[i])
        if dt < float('inf'):
            temp = t0 + dt
            heappush(collisions, Collision(temp, np.array([i]), np.array([particles[i].colcount]), np.array([v1])))

    avg_colcount = calculate_avg_colcount(particles)

    # --- Calculate initial total energy.
    energy0 = 0.0
    for p in particles:
        energy0 += 0.5*p.mass*(p.vx**2 + p.vy**2)

    # --- Identify the earliest collision.
    #t, involved, colcount, new_velocities = heappop(collisions)
    collision = heappop(collisions)
    #print('First collision at t = ', collision.time)

    # --- Plot start positions of all particles.
    if plotter == True:
        plt.ion()
        fig = plt.figure(figsize = (10,10))
        for p in particles:
            if p.radius < 0.1:
                c = 'r'
            else:
                c = 'b'
            circle = plt.Circle(p.position(), p.radius, color=c)
            fig.gca().add_artist(circle)
        plt.draw()
        time.sleep(0.01)
    a = particles[0].colcount
    b = particles[1].colcount
    
    # --- Loop:
    while(collisions and avg_colcount < stop_colcount):
    #for loop in range(1):

        # --- Move all particles forward in time (straight lines, constant 
        #     velocity) until the earliest collision.
        for i in range (len(particles)):
            particles[i].update(collision.time - t0)

        # --- For the particle(s) involved in the collision, calculate new 
        #     velocities.
        for i in range(len(collision.involved)):
            particles[collision.involved[i]].increment_colcount()
            particles[collision.involved[i]].set_velocity(collision.new_velocities[i])

        avg_colcount = calculate_avg_colcount(particles)
        #print(avg_colcount)

        # --- Calculate total energy.
        energy = 0.0
        for p in particles:
            energy += 0.5*p.mass*(p.vx**2 + p.vy**2)
        #print('Conservation of energy: ', energy/energy0)

        # --- Set new time.
        t0 = collision.time

        #print('1: ', particles[0].position(), particles[0].velocity())
        #print('2: ', particles[1].position(), particles[1].velocity())
        #print(particles[0].colcount, particles[1].colcount)

        # --- Update priority queue.
        for i in collision.involved:

            # Calculate if and when the particle(s) involved in the collision 
            # will collide with all other particles, and store the collision times.
            for j in range(0, len(particles)):
                if j != i:
                    dt, v1, v2 = time_until_collision(particles[i], particles[j])
                    if dt < float('inf'):
                        temp = t0 + dt
                        heappush(collisions, Collision(temp, np.array([i, j]), np.array([particles[i].colcount, particles[j].colcount]), np.array([v1, v2])))

            # Calculate when the particle(s) involved in the collision
            # will collide with a wall, and store the collision times.
            dt, v1 = time_until_wall_collision(particles[i])
            if dt < float('inf'):
                temp = t0 + dt
                heappush(collisions, Collision(temp, np.array([i]), np.array([particles[i].colcount]), np.array([v1])))

        # --- Identify the earliest collision.
        collision = heappop(collisions)
        while not (is_colcountOK(collision.involved, particles, collision.colcount)):
            #print('INVALID!')
            collision = heappop(collisions)
        #print('Collision at t = ', t)

        # --- Having resolved the collision, identify the new earliest collision,
        #     and check if it is still valid (if the particle(s) involved have 
        #     collided with something since the collision was calculated, it is 
        #     no longer valid). If the collision is invalid, discard it and move
        #     to the next earliest. Once a valid collision is found, repeat the 
        #     steps in the loop.

        # --- Plot current positions of particles.
        if plotter == True:
            plt.clf()
            for p in particles:
                if p.radius < 0.1:
                    c = 'r'
                else:
                    c = 'b'
                circle = plt.Circle(p.position(), p.radius, color=c)
                fig.gca().add_artist(circle)
            plt.draw()
            time.sleep(0.01)

    return particles, t0

def problem1():
    # in run: run loop only once
    stop_colcount = 0.5

    Nb = 1000
    angles = np.zeros(Nb)
    i=0
    bs = np.linspace(0.0, 0.101, Nb)
    for i in range(Nb):
        particles = [Particle(0.5, 0.5, 0.0, 0.0, 1.e-1, 1.e6), Particle(0.1, 0.5+bs[i], 1.0, 0.0, 1.e-3, 1.e0)]
        particles, time = run(particles, 0.0, stop_colcount, False)
        angle = np.arccos(np.dot(np.array([1.0,0.0]), particles[1].velocity())/1.0/np.sqrt(np.dot(particles[1].velocity(), particles[1].velocity())))
        angles[i] = angle*(180/np.pi)
    angles[-1] = 0.0

    plt.figure()
    plt.plot(np.divide(bs, 0.101), angles, 'bo', label='empty')
    plt.axis([0.0, 1.0, 0.0, 180.0])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlabel('Impact parameter $b$ [m/$R_{12}$]', fontsize=20); plt.ylabel('Scattering angle [deg]', fontsize=20)    
    plt.savefig('impact_parameter_02.png')
    plt.show()

def problem2():
    # --- Initialisation of particles.
    Nparticles = 500
    v0 = 1.0
    radius = 1.e-3
    mass = 1.0
    particles1 = []
    print('Initialising 1.')
    for i in range(Nparticles):
        x  = np.random.uniform(0.0+radius, 1.0-radius)
        y  = np.random.uniform(0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles1)):
            x  = np.random.uniform(0.0+radius, 1.0-radius)
            y  = np.random.uniform(0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles1.append(Particle(x, y, vx, vy, radius, mass))
    print('Initialising 2.')
    particles2 = []
    for i in range(Nparticles):
        x  = np.random.uniform(0.0+radius, 1.0-radius)
        y  = np.random.uniform(0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles2)):
            x  = np.random.uniform(0.0+radius, 1.0-radius)
            y  = np.random.uniform(0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles2.append(Particle(x, y, vx, vy, radius, mass))
    print('Initialising 3.')
    particles3 = []
    for i in range(Nparticles):
        x  = np.random.uniform(0.0+radius, 1.0-radius)
        y  = np.random.uniform(0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles3)):
            x  = np.random.uniform(0.0+radius, 1.0-radius)
            y  = np.random.uniform(0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles3.append(Particle(x, y, vx, vy, radius, mass))

    # --- Calculate speed of particles.
    speeds = get_speed_vector(particles1+particles2+particles3)

    # --- Plot initial speed distribution.
    plot_speed_distribution(speeds, 'speed_distribution_initial.png')

    # --- Run system.
    stop_colcount = 100
    print('Starting 1.')
    particles1, time = run(particles1, 0.0, stop_colcount)
    print('Starting 2.')
    particles2, time = run(particles2, 0.0, stop_colcount)
    print('Starting 3.')
    particles3, time = run(particles3, 0.0, stop_colcount)

    # --- Calculate speed of particles.
    speeds = get_speed_vector(particles1+particles2+particles3)

    # --- Save speeds to file just in case.
    print('Save speeds to file: ', 'my_speeds_Nparticles500.dat')
    fout = open('my_speeds_temp.dat', 'ab'); np.savetxt(fout, speeds); fout.close()
    
    # --- Load speeds from file.
    speeds = np.loadtxt('my_speeds_Nparticles500.dat')

    # --- Plot final speed distribution.
    plot_speed_distribution(speeds, 'speed_distribution_final.png')

def problem3():

    # --- Initialisation of particles.
    Nparticles = 2*500
    v0 = 1.0
    radius = 1.e-3
    mass = 1.0
    particles1 = []
    print('Initialising 1.')
    for i in range(Nparticles/2):
        x  = np.random.uniform(0.0+radius, 1.0-radius)
        y  = np.random.uniform(0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles1)):
            x  = np.random.uniform(0.0+radius, 1.0-radius)
            y  = np.random.uniform(0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles1.append(Particle(x, y, vx, vy, radius, mass))
    for i in range(Nparticles/2):
        x  = np.random.uniform(0.0+radius, 1.0-radius)
        y  = np.random.uniform(0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, 4*mass), particles1)):
            x  = np.random.uniform(0.0+radius, 1.0-radius)
            y  = np.random.uniform(0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles1.append(Particle(x, y, vx, vy, radius, 4*mass))

    # --- Calculate speed of particles.
    speeds_light = get_speed_vector(particles1[:Nparticles/2])
    speeds_heavy = get_speed_vector(particles1[Nparticles/2:])

    # --- Plot initial speed distribution.
    plot_speed_distribution(speeds_light, 'mixture_light_speed_distribution_initial.png')
    plot_speed_distribution(speeds_heavy, 'mixture_heavy_speed_distribution_initial.png')

    # --- Run system.
    stop_colcount = 100
    print('Starting 1.')
    particles1, time = run(particles1, 0.0, stop_colcount, False)

    # --- Calculate speed of particles.
    speeds_light = get_speed_vector(particles1[:Nparticles/2])
    speeds_heavy = get_speed_vector(particles1[Nparticles/2:])

    # --- Save speeds to file just in case.
    print('Save speeds to file: ', 'my_light_speeds.dat', 'and', 'my_heavy_speeds.dat')
    fout = open('my_light_speeds.dat', 'ab'); np.savetxt(fout, speeds_light); fout.close()
    fout = open('my_heavy_speeds.dat', 'ab'); np.savetxt(fout, speeds_heavy); fout.close()

    # --- Load final speeds from file.
    speeds_avg_light = np.loadtxt('my_light_speeds.dat')
    speeds_avg_heavy = np.loadtxt('my_heavy_speeds.dat')
    print('Number of speeds each: ', len(speeds_light), len(speeds_avg_heavy))

    # --- Calculate average speed.
    speeds_avg_light = np.sum(speeds_light)/(len(speeds_light))
    speeds_avg_heavy = np.sum(speeds_heavy)/(len(speeds_heavy))
    print('Average speed: ', speeds_avg_light, speeds_avg_heavy)

    # --- Calculate average kinetic energy.
    energy_avg_light = 0.5*mass*np.sum(np.power(speeds_light, 2))/(len(speeds_light))
    energy_avg_heavy = 0.5*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))
    print('Average kinetic energy: ', energy_avg_light, energy_avg_heavy)

    # --- Plot final speed distribution.
    plot_speed_distribution(speeds_light, 'mixture_light_speed_distribution_final.png')
    plot_speed_distribution(speeds_heavy, 'mixture_heavy_speed_distribution_final.png')

def problem4():

    # --- Initialisation of particles.
    print('Initialising 1.')
    Nparticles = 2*1000
    v0 = 1.0
    radius = 1.e-3
    mass = 1.0
    particles = []
    for i in range(Nparticles/2):
        x  = np.random.uniform(0.0+radius, 1.0-radius)
        y  = np.random.uniform(0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles)):
            x  = np.random.uniform(0.0+radius, 1.0-radius)
            y  = np.random.uniform(0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles.append(Particle(x, y, vx, vy, radius, mass))
    for i in range(Nparticles/2):
        x  = np.random.uniform(0.0+radius, 1.0-radius)
        y  = np.random.uniform(0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, 4*mass), particles)):
            x  = np.random.uniform(0.0+radius, 1.0-radius)
            y  = np.random.uniform(0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles.append(Particle(x, y, vx, vy, radius, 4*mass))

    # --- Calculate initial speed of particles.
    speeds_light = get_speed_vector(particles[:Nparticles/2])
    speeds_heavy = get_speed_vector(particles[Nparticles/2:])

    # --- Initialise time and kinetic energy arrays.
    Nsamples = 100
    times = np.zeros(Nsamples); energy_avg_light = np.zeros(Nsamples); energy_avg_heavy = np.zeros(Nsamples)
    times[0] = 0.0
    energy_avg_light[0] = 0.5*mass*np.sum(np.power(speeds_light, 2))/(len(speeds_light))
    energy_avg_heavy[0] = 0.5*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))

    # --- Run system.
    for i, stop_colcount in enumerate(np.linspace(1.0, 20, Nsamples-1)):#enumerate(np.logspace(0.0, 1.0, Nsamples-1)):
        print(i+1)
        particles, times[i+1] = run(particles, times[i], stop_colcount, False)

        # --- Calculate speed of particles.
        speeds_light = get_speed_vector(particles[:Nparticles/2])
        speeds_heavy = get_speed_vector(particles[Nparticles/2:])

        # --- Calculate average kinetic energy.
        energy_avg_light[i+1] = 0.5*mass*np.sum(np.power(speeds_light, 2))/(len(speeds_light))
        energy_avg_heavy[i+1] = 0.5*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))

    # --- Plot development of kinetic energy over time.
    fig = plt.figure()
    plt.scatter(times, energy_avg_light, s=50, facecolors='none', edgecolors='b', label='Light particles $m=m_0$')
    plt.scatter(times, energy_avg_heavy, s=50, facecolors='none', edgecolors='r', label='Heavy particles $m=4m_0$')
    #plt.plot(times, energy_avg_light, label='Light particles $m=m_0$')
    #plt.plot(times, energy_avg_heavy, label='Heavy particles $m=4m_0$')
    plt.xlabel(r'Time $t$', fontsize=20)
    plt.ylabel(r'Kinetic energy $E_k$', fontsize=20)
    plt.legend(loc='upper left')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.axis([min(times), max(times), 0.9*min(energy_avg_heavy), 1.1*max(energy_avg_light)])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('kinetic_energy.png')

def problem5():
    print('Hei')
    

# ---------------------------------------------- #
# -------------------- MAIN -------------------- #
# ---------------------------------------------- #
if __name__ == "__main__":

    # --- Initialise elastisity constants.
    xi = 1.0

    # --- PROBLEM 1: Calculate angle as a function of impact parameter.
    #problem1()

    # --- PROBLEM 2: Speed distribution with one particle type.
    #problem2()

    # --- PROBLEM 3: Speed distribution with two different particles.
    #problem3()

    # --- PROBLEM 4: Speed distribution with two different particles.
    #problem4()

    # --- PROBLEM 5: Crater formation.
    problem5()

    


   

    
 

    # --- Test of heapq.
    '''items = [(3, "Clear drains"), (4, "Feed cat"), (5, "Make tea"), (1, "Solve RC tasks"), (2, "Tax return")]
    heappush(items, (3, "test"))
    heapify(items)

    while items:
        print(heappop(items))'''