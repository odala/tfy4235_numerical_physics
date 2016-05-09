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

def is_overlap(particle, particles):
    for p in particles:
        distance = np.sqrt((particle.x - p.x)**2 + (particle.y - p.y)**2)
        if distance < particle.radius + p.radius:
            return True
    return False

def plot_speed_distribution(particles, filename):
    
    # --- Calculate speed of particles.
    Nparticles = len(particles)
    speeds = np.zeros(Nparticles)
    for i in range(Nparticles):
        speeds[i] = np.sqrt(particles[i].vx**2 + particles[i].vy**2)

    # --- Save speeds to file just in case.
    fout = open('my_speeds.dat', 'ab')
    np.savetxt(fout, speeds)
    fout.close()

    # --- Divide data into Nbins.
    Nbins = 101
    minimum = -1.0; maximum = 3.0
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

def calculate_avg_colcount(particles):
    avg = 0.0
    for p in particles:
        avg += p.colcount
    return avg/len(particles)

def is_colcountOK(involved, particles, colcount):
    for i in range(len(involved)):
        if colcount[i] != particles[involved[i]].colcount:
            return False
    return True

# ------------------------------------------------------ #
# ------------------------ MAIN ------------------------ #
# ------------------------------------------------------ #
def main(particles, plotter=False):

    # --- Initialise square box.
    '''xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0'''

    # --- Set start time.
    t0 = 0.0

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
    print('First collision at t = ', collision.time)

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
    while(collisions and avg_colcount < 100):
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

    return particles

if __name__ == "__main__":

    # --- Initialise elastisity constants.
    xi = 1.0

    # --- PROBLEM 2: speed distribution.
    # --- Initialisation of particles.
    Nparticles = 2000
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
    plot_speed_distribution(particles1+particles2+particles3, 'speed_distribution_initial.png')
    print('Starting 1.')
    particles1 = main(particles1, False)
    print('Starting 2.')
    particles2 = main(particles2, False)
    print('Starting 3.')
    particles3 = main(particles3, False)
    plot_speed_distribution(particles1+particles2+particles3, 'speed_distribution_final.png')

    # --- PROBLEM 1: Calculate angle as a function of impact parameter.
    '''
    # in main: run loop only once
    Nb = 1000
    angles = np.zeros(Nb)
    i=0
    bs = np.linspace(0.0, 0.101, Nb)
    for i in range(Nb):
        particles = [Particle(0.5, 0.5, 0.0, 0.0, 1.e-1, 1.e6), Particle(0.1, 0.5+bs[i], 1.0, 0.0, 1.e-3, 1.e0)]
        particles = main(particles, False)
        angle = np.arccos(np.dot(np.array([1.0,0.0]), particles[1].velocity())/1.0/np.sqrt(np.dot(particles[1].velocity(), particles[1].velocity())))
        angles[i] = angle*(180/np.pi)
    angles[-1] = 0.0

    plt.figure()
    plt.plot(np.divide(bs, 0.101), angles, 'bo', label='empty')
    plt.axis([0.0, 1.0, 0.0, 180.0])
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlabel('Impact parameter $b$ [m/$R_{12}$]', fontsize=20); plt.ylabel('Scattering angle [deg]', fontsize=20)    
    plt.savefig('impact_parameter_02.png')
    plt.show()'''
 

    # --- Test of heapq.
    items = [(3, "Clear drains"), (4, "Feed cat"), (5, "Make tea"), (1, "Solve RC tasks"), (2, "Tax return")]
    heappush(items, (3, "test"))
    heapify(items)

    while items:
        print(heappop(items))