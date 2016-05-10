# Diffusion
# from myIterativeMethods import explicit_euler, implicit_euler, ...

import numpy as np
import copy
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

# --- Function: Calculates kinetic energy of all particles.
#     input : list of Particles.
#     output: list of kinetic energies.
def get_kinetic_energy_array(particles):
    Nparticles = len(particles)
    kinetic_energies = np.zeros(Nparticles)
    for i in range(Nparticles):
        kinetic_energies[i] = 0.5*particles[i].mass*(particles[i].vx**2 + particles[i].vy**2)
    return kinetic_energies

def get_average_kinetic_energy(particles):
    return sum(get_kinetic_energy_array(particles))/Nparticles

def get_total_kinetic_energy(particles):
    return sum(get_kinetic_energy_array(particles))

# --- Function: Calculated the size of the crater by comparing how
#     many of the particles that have moved more than one radius 
#     away from its initial starting position.
#     input : list of Particles
#     output: the crater size as a number between 0.0 and 1.0
#               - 0.0 if none of the particles got moved
#               - 1.0 if all of the particles got moved
def get_crater_size(initial_particles, final_particles):
    affected = 0.0
    for i in range(len(initial_particles)):
        ip = initial_particles[i]
        fp = final_particles[i]
        distance = np.sqrt((fp.x - ip.x)**2 + (fp.y - ip.y)**2)
        if distance > ip.radius:
            affected += 1.0
    return affected/len(initial_particles)

# --- Function: Plot the particles.
#     input : list of Particles
#     output: none, but save a figure of the particles.
def plot_particles(particles, filename):
    fig = plt.figure(figsize = (10,10))
    for p in particles:
        if p.radius < 1.e-2:
            c = 'b'
        else:
            c = 'r'
        circle = plt.Circle(p.position(), p.radius, color=c)
        fig.gca().add_artist(circle)
    
    # --- Saving figure.
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.tight_layout()
    plt.savefig(filename)

# --- Function: Run the event-driven simulation of particles.
#     input : 
#     output: 
def run(particles, start_time, stop_colcount, plotter=False):

    # --- Set start time.
    t0 = start_time

    # --- Initialise priority queue.
    collisions = []
    for i in range (len(particles)):

        # Calculate if and when particle will collide with all other particles,
        # and store the collision times.
        for j in range(i+1, len(particles)):
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
    collision = heappop(collisions)

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
    
    # --- Loop:
    while(collisions and avg_colcount < stop_colcount):
    #for loop in range(1):

        # --- Move all particles forward in time (straight lines, constant 
        #     velocity) until the earliest collision.
        for i in range (len(particles)):
            particles[i].update(collision.time - t0)

        # --- Set new time.
        t0 = collision.time

        # --- For the particle(s) involved in the collision, increment
        #     colcount and calculate new velocities.
        for i in range(len(collision.involved)):
            particles[collision.involved[i]].increment_colcount()
            particles[collision.involved[i]].set_velocity(collision.new_velocities[i])
        
        # --- Update average colcount.
        avg_colcount = calculate_avg_colcount(particles)

        # --- Calculate total energy.
        energy = 0.0
        for p in particles:
            energy += 0.5*p.mass*(p.vx**2 + p.vy**2)
        #print('Conservation of energy: ', energy/energy0)        

        # --- Update priority queue for the particle(s) involved in the collision.
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

        # --- Identify the new earliest collision, and check if it is still valid. 
        #     If the collision is invalid, discard it and move to the next earliest, and so on.
        collision = heappop(collisions)
        while not (is_colcountOK(collision.involved, particles, collision.colcount)):
            collision = heappop(collisions)

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
    # in run(): run loop only once
    stop_colcount = 0.5
    start_time = 0.0

    Nb = 1000
    angles = np.zeros(Nb)
    for i, b in enumerate(np.linspace(0.0, 0.101, Nb)):
        particles = [Particle(0.5, 0.5, 0.0, 0.0, 1.e-1, 1.e6), Particle(0.1, 0.5+b, 1.0, 0.0, 1.e-3, 1.e0)]
        particles, time = run(particles, start_time, stop_colcount)
        angle = np.arccos(np.dot(np.array([1.0,0.0]), particles[1].velocity())/1.0/np.sqrt(np.dot(particles[1].velocity(), particles[1].velocity())))
        angles[i] = angle*(180/np.pi)
    angles[-1] = 0.0

    # --- Plot scattering angle as a function of impact parameter.
    plt.figure()
    plt.plot(np.divide(bs, 0.101), angles, 'bo', label='empty')
    plt.axis([0.0, 1.0, 0.0, 180.0])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlabel('Impact parameter $b$ [m/$R_{12}$]', fontsize=20); plt.ylabel('Scattering angle [deg]', fontsize=20)

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('impact_parameter_02.png')

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
    energy_avg_heavy = 0.5*4*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))
    print('Average kinetic energy: ', energy_avg_light, energy_avg_heavy)

    # --- Plot final speed distribution.
    plot_speed_distribution(speeds_light, 'mixture_light_speed_distribution_final.png')
    plot_speed_distribution(speeds_heavy, 'mixture_heavy_speed_distribution_final.png')

def problem4():

    # --- Initialisation of particles.
    print('Initialising 1.')
    Nparticles = 2*500
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
    Nsamples = 150
    times = np.zeros(Nsamples); energy_avg_light = np.zeros(Nsamples); energy_avg_heavy = np.zeros(Nsamples)
    times[0] = 0.0
    energy_avg_light[0] = 0.5*mass*np.sum(np.power(speeds_light, 2))/(len(speeds_light))
    energy_avg_heavy[0] = 0.5*4*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))

    # --- Run system.
    for i, stop_colcount in enumerate(np.linspace(2/Nparticles, 20, Nsamples-1)):#enumerate(np.logspace(0.0, 1.0, Nsamples-1)):
        print(i+1)
        particles, times[i+1] = run(particles, times[i], stop_colcount, False)

        # --- Calculate speed of particles.
        speeds_light = get_speed_vector(particles[:Nparticles/2])
        speeds_heavy = get_speed_vector(particles[Nparticles/2:])

        # --- Calculate average kinetic energy.
        energy_avg_light[i+1] = 0.5*mass*np.sum(np.power(speeds_light, 2))/(len(speeds_light))
        energy_avg_heavy[i+1] = 0.5*4*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))

    # --- Plot development of kinetic energy over time.
    fig = plt.figure()
    plt.scatter(times, energy_avg_light, s=50, facecolors='none', edgecolors='b', label='Light particles')
    plt.scatter(times, energy_avg_heavy, s=50, facecolors='none', edgecolors='r', label='Heavy particles')
    plt.scatter(times, (energy_avg_light+energy_avg_heavy)/2, s=25, facecolors='none', edgecolors='g', label='All particles')
    #plt.plot(times, energy_avg_light, label='Light particles $m=m_0$')
    #plt.plot(times, energy_avg_heavy, label='Heavy particles $m=4m_0$')
    plt.xlabel(r'Time $t$', fontsize=20)
    plt.ylabel(r'Kinetic energy $E_k$', fontsize=20)
    plt.legend(loc='upper right')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.axis([min(times), max(times), 0.9*min(min(energy_avg_light),min(energy_avg_heavy)), 1.1*max(max(energy_avg_light), max(energy_avg_heavy))])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('kinetic_energy.png')

def problem5():

    # --- Initialisation of wall of particles in (0.0, 1.0)x(0.0, 0.5).
    print('Initialising particles.')
    Nparticles = 5000
    radius1 = 4.e-3 #np.sqrt(0.5*0.5/(Nparticles*np.pi)) #
    mass1 = 1.0
    initial_particles = []
    for i in range(Nparticles/2):
        x  = np.random.uniform(0.0+radius1, 1.0-radius1)
        y  = np.random.uniform(0.0+radius1, 0.5-radius1)
        while (is_overlap(Particle(x, y, 0.0, 0.0, radius1, mass1), initial_particles)):
            x  = np.random.uniform(0.0+radius1, 1.0-radius1)
            y  = np.random.uniform(0.0+radius1, 0.5-radius1)
        initial_particles.append(Particle(x, y, 0.0, 0.0, radius1, mass1))

    # --- Append the projectile with larger mass and larger radius.
    radius2 = 5*radius1
    mass2 = 25*mass1
    v0 = 5.0
    initial_particles.append(Particle(0.5, 0.75, 0, -v0, radius2, mass2))

    # --- Plot initial positions of all particles.
    plot_particles(initial_particles, 'crater_formation_initial.png')

    # --- Run the simulation several times to get a
    #     parameter scan of the mass of the projectile.
    Nfactors = 10
    crater_size = np.zeros(Nfactors)
    factors = np.linspace(1, 30, Nfactors)
    for i, f in enumerate(factors):
        print('Run for mass of projectile ', f, ' times bigger than the mass of the smaller particles.')
        particles = copy.deepcopy(initial_particles)
        particles[-1].mass = f*mass1

        # --- Calculate intial energy.
        initial_energy = get_total_kinetic_energy(particles)
        energy = initial_energy

        # --- Plot initial positions of all particles.
        plot_particles(initial_particles, 'crater_formation_initial' + str(i) + '.png')

        # --- Run simulation until only 10 % of the initial energy remains.
        stop_colcount = 3.0
        time = 0.0
        while (energy/initial_energy > 0.1):
            particles, time = run(particles, time, stop_colcount, False)

            # --- Calculate total energy.
            energy = get_total_kinetic_energy(particles)
            print('Fraction of energy left: ', energy/initial_energy)
            stop_colcount += 1.e-1

        # --- Plot final positions of all particles.
        plot_particles(particles, 'crater_formation_final' + str(i) + '.png')

        # --- Calculate size of crater.
        crater_size[i] = get_crater_size(initial_particles[:Nparticles], particles[:Nparticles])

        # --- Save crater size and corresponding mass factor to file (just in case).
        fout = open('crater_size_vs_mass.out','ab'); np.savetxt(fout, np.atleast_2d(np.array([f, crater_size[i]]))); fout.close()

    # --- Plot size of crater as a function of mass.
    fig = plt.figure()
    plt.scatter(factors, crater_size, s=50, facecolors='none', edgecolors='b')
    plt.xlabel(r'Mass of projectile $M$ [kg/mass of small particles] ', fontsize=20)
    plt.ylabel(r'Crater size', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    #plt.axis([min(factors), max(factors), 0.9*min(crater_size), 1.1*max(crater_size)])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('crater_size_vs_mass_of_projectile.png')

    
# ---------------------------------------------- #
# -------------------- MAIN -------------------- #
# ---------------------------------------------- #
if __name__ == "__main__":

    # --- Initialise elastisity constants.
    xi = 0.5

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