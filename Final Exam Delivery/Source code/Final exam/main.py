# Final Exam in Computational Physics
# start: 09.05.16 09:00
# end  : 12.05.16 23:00

import numpy as np
import copy
from matplotlib import pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
import time
import multiprocessing as mp
from heapq import heappush, heappop, heapify
from scipy import constants
from scipy import stats
from os.path import isfile

# --- Import my functions.
from myParticle import Particle
from myCollision import Collision 
   
# --- Function: Calculate time until collision with a wall.
#     input : Particle
#     output: time until collision with wall, a number representing
#             if it is horizontal (-2) or vertical wall (-1).
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
        return timex, -1
    else:
        return timey, -2

# --- Calculate time until collision between particles.
#     input : two Particles
#     output: time until collision between particles
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
        return dt

# --- Calculate new velocities for the involved particle(s)
#     in a collision.
#     input : the involved entities in the collision
#             - if particle-particle collision: two Particles 
#             - otherwise: a Particle and a number representing which wall.
#     output: time until collision between particles
def get_new_velocity(particle1, particle2, xi):
    if isinstance(particle2, Particle):
        dv  = particle2.velocity() - particle1.velocity()
        dx = particle2.position() - particle1.position()
        R  = particle1.radius + particle2.radius
        velocity1 = particle1.velocity() + ((1.0+xi)*particle2.mass/(particle1.mass+particle2.mass)*np.dot(dv,dx)/R**2)*dx
        velocity2 = particle2.velocity() - ((1.0+xi)*particle1.mass/(particle1.mass+particle2.mass)*np.dot(dv,dx)/R**2)*dx
        return velocity1, velocity2
    elif (particle2 == -1):
        return np.array([-xi*particle1.vx, xi*particle1.vy])
    elif (particle2 == -2):
        return np.array([xi*particle1.vx, -xi*particle1.vy])
    else:
        print('Something went wrong in the get_new_velocity-function.')

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
    if colcount[1] < 0:
        Ninv = 1
    else:
        Ninv = 2
    for i in range(Ninv):
        if colcount[i] != particles[involved[i]].colcount: 
            return False
    return True

# --- Function: Draws a uniform random number x in [xa, xb) 
#     and a uniform random number y in [yc, yd).
#     input : xa, xb, yc, yd
#     output: x, y
def draw_uniform_x_y(xa, xb, yc, yd):
    return np.random.uniform(xa, xb), np.random.uniform(yc, yd)

# --- Function: Calculate if a new particle overlaps one of 
#     the existing particles.
#     input : Particle, list of Particles.
#     output: boolean.
def is_overlap(particle, particles):
    for p in particles:
        distance = np.sqrt((particle.x - p.x)**2 + (particle.y - p.y)**2)
        if distance < (particle.radius + p.radius):
            return True
    return False

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

# --- Function: Calculates the average kinetic energy.
#     input : list of Particles.
#     output: average kinetic energy.
def get_average_kinetic_energy(particles):
    return sum(get_kinetic_energy_array(particles))/Nparticles

# --- Function: Calculates total kinetic energy.
#     input : list of Particles.
#     output: total kinetic energy.
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

# --- Function: Plot the distribution of speeds in a histogram.
#     input : list of speeds, filename of savefile.
#     output: none, but saves a figure of the histogram.
def plot_speed_distribution(speeds, filename):

    # --- Divide data into Nbins.
    Nbins = 101
    minimum = 0.0; maximum = 5.0
    hist, bins = np.histogram(speeds, bins=Nbins, range=(minimum, maximum), normed=True)
    bin_values = (bins[:-1] + bins[1:]) / 2

    # --- Plot histogram.
    histogram = plt.figure()
    plt.bar(bin_values, hist, width=bins[1]-bins[0], color='blue')

    plt.xlabel(r'Speed $v$ [m/s]', fontsize=30)
    plt.ylabel(r'$f(v)$', fontsize=30)
    plt.tick_params(axis='both', which='major', labelsize=25)
    plt.axis([minimum, maximum, 0.0, 1.2]) #max(hist)

    # --- Saving figure.
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
            dt = time_until_collision(particles[i], particles[j])
            if dt < float('inf'):
                temp = t0 + dt
                heappush(collisions, Collision(temp, np.array([i, j]), np.array([particles[i].colcount, particles[j].colcount])))

        # Calculate when particle will collide with a wall,
        # and store the collision times.
        dt, inv = time_until_wall_collision(particles[i])
        if dt < float('inf'):
            temp = t0 + dt
            heappush(collisions, Collision(temp, np.array([i, inv]), np.array([particles[i].colcount, -1])))

    avg_colcount = calculate_avg_colcount(particles)

    # --- Calculate initial total energy.
    initial_energy = get_total_kinetic_energy(particles)
    energy = initial_energy

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
    while(avg_colcount < stop_colcount and energy > 0.1*initial_energy):

        # --- Move all particles forward in time (straight lines, constant 
        #     velocity) until the earliest collision.
        for i in range (len(particles)):
            particles[i].update(collision.time - t0)

        # --- Set new time.
        t0 = collision.time

        # --- For the particle(s) involved in the collision,
        #     calculate new velocities and save the collision time.
        #     Increment colcount if it was a particle-particle collision.
        if collision.involved[1] > -1:
            inv1 = collision.involved[0]
            inv2 = collision.involved[1]
            if (t0 - particles[inv1].last) < tc or (t0 - particles[inv2].last) < tc:
                v1, v2 = get_new_velocity(particles[inv1], particles[inv2], 1.0)
            else:
                v1, v2 = get_new_velocity(particles[inv1], particles[inv2], xi)
            particles[inv1].set_velocity(v1)
            particles[inv2].set_velocity(v2)
            particles[inv1].last = particles[inv2].last = t0
            particles[inv1].increment_colcount()
            particles[inv2].increment_colcount()
        else:
            inv1 = collision.involved[0]
            if (t0 - particles[inv1].last) < tc:
                v1 = get_new_velocity(particles[inv1], collision.involved[1], 1.0)
            else:
                v1 = get_new_velocity(particles[inv1], collision.involved[1], xi)
            particles[inv1].last = t0
            particles[inv1].set_velocity(v1)
        
        # --- Update average colcount.
        avg_colcount = calculate_avg_colcount(particles)

        # --- Calculate total energy.
        energy = get_total_kinetic_energy(particles)

        if int(round(avg_colcount*len(particles))) % 1000 == 0:
            print('Ncollisions = ', int(round(avg_colcount*len(particles))), '. Time: ', t0,'. Conservation of energy: ', energy/initial_energy)   
            #plot_particles(particles, 'temporary_picture.png')     

        # --- Update priority queue for the particle(s) involved in the collision.
        for i in collision.involved:
            if i < 0:
                break

            # Calculate if and when the particle(s) involved in the collision 
            # will collide with all other particles, and store the collision times.
            for j in range(0, len(particles)):
                if j != i:
                    dt = time_until_collision(particles[i], particles[j])
                    if dt < float('inf'):
                        temp = t0 + dt
                        heappush(collisions, Collision(temp, np.array([i, j]), np.array([particles[i].colcount, particles[j].colcount])))

            # Calculate when the particle(s) involved in the collision
            # will collide with a wall, and store the collision times.
            dt, inv = time_until_wall_collision(particles[i])
            if dt < float('inf'):
                temp = t0 + dt
                heappush(collisions, Collision(temp, np.array([i, inv]), np.array([particles[i].colcount, -1])))

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

# --- Run simulation of a small and light particle with velocity
#     in the positive x-direction hitting one stationary, large
#     and heavy particle for different values of the impact parameter.
#     Calculate the scattering angle as a function of impact parameter.
def problem1():
    
    stop_colcount = 0.5     # in run(): loop only once
    start_time = 0.0

    # --- Run for different values of the impact parameter.
    print('Run for different values of the impact parameter.')
    Nb = 100
    angles = np.zeros(Nb)
    bs = np.linspace(0.0, 0.101, Nb)
    for i, b in enumerate(bs):

        # --- Initialise and run particles.
        particles = [Particle(0.5, 0.5, 0.0, 0.0, 1.e-1, 1.e6), Particle(0.1, 0.5+b, 1.0, 0.0, 1.e-3, 1.e0)]
        particles, time = run(particles, start_time, stop_colcount)
        
        # --- Calculate scattering angle.
        angle = np.arccos(np.dot(np.array([1.0,0.0]), particles[1].velocity())/1.0/np.sqrt(np.dot(particles[1].velocity(), particles[1].velocity())))
        angles[i] = angle*(180/np.pi)
    angles[-1] = 0.0

    # --- Plot scattering angle as a function of impact parameter.
    print('Plot scattering angle as a function of impact parameter.')
    plt.figure()
    plt.scatter(np.divide(bs, 0.101), angles, s=50, facecolors='none', edgecolors='b', label='Numerical')
    plt.plot(np.divide(bs, 0.101), 2*np.arccos(np.divide(bs, 0.101))*180/np.pi, 'r--', label='Exact')
    plt.axis([0.0, 1.0, 0.0, 180.0])
    plt.legend(loc='upper right')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlabel('Impact parameter $b$ [m/$R_{12}$]', fontsize=20); plt.ylabel('Scattering angle [deg]', fontsize=20)

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('impact_parameter.png')

def problem2(Nparticles):

    # --- Initialise particles.
    Nparticles = Nparticles
    fileid = '_Nparticles' + str(Nparticles)
    v0 = 1.0
    radius = 1.e-3
    mass = 1.0
    particles = []
    print('Initialise particles.')
    for i in range(Nparticles):
        x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles)):
            xx, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles.append(Particle(x, y, vx, vy, radius, mass))

    # --- Calculate speed of particles.
    speeds = get_speed_vector(particles)

    # --- Plot initial speed distribution.
    figname = 'speed_distribution' + fileid + '_initial.png'
    #plot_speed_distribution(speeds, figname)

    # --- Run system.
    stop_colcount = 100
    print('Start run.')
    particles, time = run(particles, 0.0, stop_colcount)

    # --- Calculate speed of particles.
    speeds = get_speed_vector(particles)

    # --- Appends speeds to file.
    filename = 'my_speeds' + fileid + '.dat'
    print('Save speeds to file: ', filename)
    fout = open(filename, 'ab'); np.savetxt(fout, speeds); fout.close()
    
    # --- Load speeds from file.
    speeds = np.loadtxt(filename)

    # --- Plot final speed distribution.
    figname = 'speed_distribution' + fileid + '_final.png'
    plot_speed_distribution(speeds, figname)

def problem3(Nparticles):

    # --- Initialise particles.
    Nparticles = Nparticles
    fileid = '_Nparticles' + str(Nparticles)
    v0 = 1.0
    radius = 1.e-3
    mass = 1.0
    particles1 = []
    print('Initialise particles of mass 1*mass.')
    for i in range(Nparticles/2):
        x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles1)):
            x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles1.append(Particle(x, y, vx, vy, radius, mass))
    print('Initialise particles of mass 4*mass.')
    for i in range(Nparticles/2):
        x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, 4*mass), particles1)):
            x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles1.append(Particle(x, y, vx, vy, radius, 4*mass))

    # --- Plot initial speed distribution.
    print('Plot initial speed distribution.')
    speeds_light = get_speed_vector(particles1[:Nparticles/2])
    speeds_heavy = get_speed_vector(particles1[Nparticles/2:])
    figname_light = 'mixture_light_speed_distribution_initial' + fileid + '.png'
    figname_heavy = 'mixture_heavy_speed_distribution_initial' + fileid + '.png'
    plot_speed_distribution(speeds_light, figname_light)
    plot_speed_distribution(speeds_heavy, figname_heavy)

    # --- Run system until average collision count is above 100.
    print('Run system until average collision count is above 100.')
    stop_colcount = 100
    particles1, time = run(particles1, 0.0, stop_colcount, False)

    # --- Calculate speed of particles.
    speeds_light = get_speed_vector(particles1[:Nparticles/2])
    speeds_heavy = get_speed_vector(particles1[Nparticles/2:])

    # --- Save speeds to file just in case.
    filename_light = 'my_light_speeds' + fileid + '.dat'
    filename_heavy = 'my_heavy_speeds' + fileid + '.dat'
    print('Save speeds to file: ', filename_light, 'and', filename_heavy)
    fout = open(filename_light, 'ab'); np.savetxt(fout, speeds_light); fout.close()
    fout = open(filename_heavy, 'ab'); np.savetxt(fout, speeds_heavy); fout.close()

    # --- Load final speeds from file.
    speeds_light = np.loadtxt(filename_light)
    speeds_heavy = np.loadtxt(filename_heavy)
    print('Number of speeds each: ', len(speeds_light), len(speeds_heavy))

    # --- Calculate average speed.
    print('Calculate average speed.')
    speeds_avg_light = np.sum(speeds_light)/(len(speeds_light))
    speeds_avg_heavy = np.sum(speeds_heavy)/(len(speeds_heavy))
    print('Light particles = ', speeds_avg_light, '   Heavy particles = ', speeds_avg_heavy)

    # --- Calculate average kinetic energy.
    print('Calculate average kinetic energy.')
    energy_avg_light = 0.5*mass*np.sum(np.power(speeds_light, 2))/(len(speeds_light))
    energy_avg_heavy = 0.5*4*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))
    print('Light particles = ', energy_avg_light, '   Heavy particles = ', energy_avg_heavy)

    # --- Plot final speed distribution.
    print('Plot final speed distribution.')
    figname_light = 'mixture_light_speed_distribution_final' + fileid + '.png'
    figname_heavy = 'mixture_heavy_speed_distribution_final' + fileid + '.png'
    plot_speed_distribution(speeds_light, figname_light)
    plot_speed_distribution(speeds_heavy, figname_heavy)

def problem4(Nparticles):

    # --- Initialisation of particles.
    Nparticles = Nparticles
    fileid = '_Nparticles' + str(Nparticles)
    v0 = 1.0
    radius = 1.e-3
    mass = 1.0
    particles = []
    print('Initialise particles of mass 1*mass.')
    for i in range(Nparticles/2):
        x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, mass), particles)):
            x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        theta = np.random.uniform(0.0, 2*np.pi)
        vx = v0*np.cos(theta)
        vy = v0*np.sin(theta)
        particles.append(Particle(x, y, vx, vy, radius, mass))
    print('Initialise particles of mass 4*mass.')
    for i in range(Nparticles/2):
        x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
        while (is_overlap(Particle(x, y, v0, v0, radius, 4*mass), particles)):
            x, y = draw_uniform_x_y(0.0+radius, 1.0-radius, 0.0+radius, 1.0-radius)
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
    print('Run system until average collision count reach 30.')
    for i, stop_colcount in enumerate(np.linspace(2/Nparticles, 20, Nsamples-1)):
        print(i+1)
        particles, times[i+1] = run(particles, times[i], stop_colcount, False)

        # --- Calculate speed of particles.
        speeds_light = get_speed_vector(particles[:Nparticles/2])
        speeds_heavy = get_speed_vector(particles[Nparticles/2:])

        # --- Calculate average kinetic energy.
        energy_avg_light[i+1] = 0.5*mass*np.sum(np.power(speeds_light, 2))/(len(speeds_light))
        energy_avg_heavy[i+1] = 0.5*4*mass*np.sum(np.power(speeds_heavy, 2))/(len(speeds_heavy))

    # --- Plot development of kinetic energy over time.
    print('Plot development of kinetic energy over time.')
    fig = plt.figure()
    plt.scatter(times, energy_avg_light, s=50, facecolors='none', edgecolors='b', label='Light particles')
    plt.scatter(times, energy_avg_heavy, s=50, facecolors='none', edgecolors='r', label='Heavy particles')
    plt.scatter(times, (energy_avg_light+energy_avg_heavy)/2, s=25, facecolors='none', edgecolors='g', label='All particles')
    plt.xlabel(r'Time $t$', fontsize=20)
    plt.ylabel(r'Kinetic energy $E_k$', fontsize=20)
    plt.legend(loc='upper right')
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.axis([min(times), max(times), 0.9*min(min(energy_avg_light),min(energy_avg_heavy)), 1.1*max(max(energy_avg_light), max(energy_avg_heavy))])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('kinetic_energy' + fileid + '.png')

def problem5(Nparticles, scan_across):

    # --- Initialise the wall of particles in (0.0, 1.0)x(0.0, 0.5).
    print('Initialising particles.')
    Nparticles = Nparticles
    fileid = '_Nparticles' + str(Nparticles)
    radius1 = np.sqrt(0.5*0.5/(Nparticles*np.pi)) # 4e-3
    mass1 = 1.0
    initial_particles = []
    posfilename = 'wall_of_particles' + fileid + '.out'
    if isfile(posfilename):
        positions = np.loadtxt(posfilename)
        for i in range(Nparticles):
            x = positions[i,0]
            y = positions[i,1]

            # --- Append Particle to list.
            initial_particles.append(Particle(x, y, 0.0, 0.0, radius1, mass1))
    else:
        for i in range(Nparticles):
            x, y = draw_uniform_x_y(0.0+radius1, 1.0-radius1, 0.0+radius1, 0.5-radius1)
            while (is_overlap(Particle(x, y, 0.0, 0.0, radius1, mass1), initial_particles)):
                x, y = draw_uniform_x_y(0.0+radius1, 1.0-radius1, 0.0+radius1, 0.5-radius1)
        
            # --- Save positions to file for use later and then append to initial particles.
            fout = open(posfilename,'ab'); np.savetxt(fout, np.atleast_2d(np.array([x, y]))); fout.close()

            # --- Append Particle to list.
            initial_particles.append(Particle(x, y, 0.0, 0.0, radius1, mass1))

    # --- Append the projectile with larger mass and larger radius.
    radius2 = 5*radius1
    mass2 = 25*mass1
    v0 = 5.0
    initial_particles.append(Particle(0.5, 0.75, 0.0, -v0, radius2, mass2))

    print('Nparticles = ', len(initial_particles))

    # --- Plot initial positions of all particles.
    figname = 'crater_formation_initial_' + scan_across + fileid + '.png'
    #splot_particles(initial_particles, figname)

    # --- Run the simulation several times to get a
    #     parameter scan of the mass of the projectile.
    Nfactors = 10
    crater_size = np.zeros(Nfactors)
    if scan_across == 'mass':
        factors = np.linspace(1e-3, 35, Nfactors)
    elif scan_across == 'radius':
        factors = np.linspace(1e-3, 30, Nfactors)
    elif scan_across == 'velocity':
        factors = np.linspace(1e-3, 10, Nfactors)
    
    filename = 'crater_size_vs_' + scan_across + fileid + '.out'

    # MULTIPROCESSING START
    #number_of_cores = mp.cpu_count()
    queue = mp.Queue()
    sub_process = [None] * Nfactors
    for i, f in enumerate(factors):

        # --- Initialise particles to start values.
        particles = copy.deepcopy(initial_particles)

        # --- Set the variable for which the parameter scan is done for.
        if scan_across == 'mass':
            print(i, 'Run for mass of projectile ', f, ' times bigger than the small mass.')
            particles[-1].mass = f*mass1
        elif scan_across == 'radius':
            print(i, ': Run for radius of projectile ', f, ' times bigger than the small radius.')
            particles[-1].radius = f*radius1
        elif scan_across == 'velocity':
            print(i, ': Run for velocity of projectile ', f)
            particles[-1].vy = -f

        # --- Run parameter scan.
        figname = 'crater_formation_final_' + scan_across + fileid + '_' + str(i) + '.png'
        sub_process[i] = mp.Process(target=parameter_scan, args=(initial_particles, particles, f, Nparticles, filename, figname, queue))
        sub_process[i].start()
    for i in range(Nfactors):
        crater_size[i] = queue.get()
    for i in range(Nfactors):
        sub_process[i].join()
    # MULTIPROCESSING END

    # --- Load data from file.
    data        = np.loadtxt(filename)
    factors     = data[:,0]
    crater_size = data[:,1]

    # --- Plot size of crater as a function of scan parameter.
    fig = plt.figure()
    plt.scatter(factors, crater_size, s=50, facecolors='none', edgecolors='b')
    if scan_across == 'mass':
        plt.xlabel(r'Mass of projectile $M$/mass of small particles', fontsize=20)
    elif scan_across == 'radius':
        plt.xlabel(r'Radius of projectile $R$/radius of small particles', fontsize=20)
    elif scan_across == 'velocity':
        plt.xlabel(r'velocity of projectile $v_0$', fontsize=20)
    plt.ylabel(r'Crater size', fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.axis([min(factors), max(factors), 0.0, 1.1*max(crater_size)])

    # --- Saving figure.
    plt.tight_layout()
    plt.savefig('crater_size_vs_' + scan_across + fileid + '.png')

# --- Function to multiprocess.
def parameter_scan(initial_particles, particles, f, Nparticles, filename, figname, queue):

    # --- Run simulation until only 10 % of the initial energy remains.
    stop_colcount = float('inf'); start_time = 0.0
    particles, time = run(particles, start_time, stop_colcount)

    # --- Plot final positions of all particles.
    #plot_particles(particles, figname)

    # --- Calculate size of crater.
    crater_size = get_crater_size(initial_particles[:Nparticles], particles[:Nparticles])

    # --- Save crater size and corresponding mass factor to file (just in case).
    print('Save to file: ', filename)
    fout = open(filename,'ab'); np.savetxt(fout, np.atleast_2d(np.array([f, crater_size]))); fout.close()

    queue.put(crater_size)
    
# ---------------------------------------------- #
# -------------------- MAIN -------------------- #
# ---------------------------------------------- #
if __name__ == "__main__":

    # --- Set start time of simulation.
    tic = time.clock()

    # --- Initialise elastisity constants.
    xi = 1.0
    print('xi = ', xi)

    # --- Set duration of contact.
    tc = 1e-5

    # --- PROBLEM 1: Calculate angle as a function of impact parameter.
    #problem1() DONEDONEDONEDONEDONEDONEDONE

    # --- PROBLEM 2: Speed distribution with one particle type.
    Nparticles = 1000
    #problem2(Nparticles) DONEDONEDONEDONEDONEDONEDONE

    # --- PROBLEM 3: Speed distribution with two different particles.
    Nparticles = 1000
    #problem3(Nparticles)

    # --- PROBLEM 4: Speed distribution with two different particles.
    xi = 0.8
    print('Problem 4: xi = ', xi)
    Nparticles = 1000
    #problem4(Nparticles) DONEDONEDONEDONEDONEDONEDONE

    # --- PROBLEM 5: Crater formation.
    xi = 0.5
    print('Problem 5: xi = ', xi)
    Nparticles  = 1000
    scan_across = 'radius' # or 'mass' or 'radius'
    fileid      = '_Nparticles' + str(Nparticles)
    #problem5(Nparticles, scan_across)

    # --- Print how long the simulation took.
    seconds = time.clock() - tic
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    print('Time = %f seconds (%02d:%02d:%02d).' % (seconds, h, m, s))