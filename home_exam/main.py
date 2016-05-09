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


# ------------------------------------------------------ #
# ------------------------ MAIN ------------------------ #
# ------------------------------------------------------ #
def main(b, plotter=False):

    # --- Initialise square box.
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0

    # --- Set start time.
    t0 = 0.0

    # --- Initialisation of particles and collision priority heap.
    particles = [Particle(0.5, 0.5, 0.0, 0.0, 1.e-1, 1.e6), Particle(0.1, 0.5+b, 1.0, 0.0, 1.e-3, 1.e0)]
    print('Start1: ', particles[0].position(), particles[0].velocity(), particles[0].radius, particles[0].mass)
    print('Start2: ', particles[1].position(), particles[1].velocity(), particles[1].radius, particles[1].mass)
    collisions = []
    for i in range (len(particles)):

        # Calculate if and when particle will collide with all other particles,
        # and store the collision times.
        for j in range(i+1, len(particles)):
            #print(i, j)
            dt, v1, v2 = time_until_collision(particles[i], particles[j])
            if dt < float('inf'):
                temp = t0 + dt
                heappush(collisions, (temp, np.array([i, j]), np.array([particles[i].colcount, particles[j].colcount]), np.array([v1, v2])))

        # Calculate when particle will collide with a wall,
        # and store the collision times.
        dt, v1 = time_until_wall_collision(particles[i])
        if dt < float('inf'):
            temp = t0 + dt
            heappush(collisions, (temp, np.array([i]), np.array([particles[i].colcount]), np.array([v1])))

    # --- Calculate initial total energy.
    energy0 = 0.0
    for p in particles:
        energy0 += 0.5*p.mass*(p.vx**2 + p.vy**2)

    # --- Identify the earliest collision.
    t, involved, colcount, new_velocities = heappop(collisions)
    #print('First collision at t = ', t)

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
    #while(collisions):
    for loop in range(1):

        # --- Move all particles forward in time (straight lines, constant 
        #     velocity) until the earliest collision.
        for i in range (len(particles)):
            particles[i].update(t - t0)

        # --- For the particle(s) involved in the collision, calculate new 
        #     velocities.
        for i in range(len(involved)):
            particles[involved[i]].increment_colcount()
            particles[involved[i]].set_velocity(new_velocities[i])

        # --- Calculate total energy.
        energy = 0.0
        for p in particles:
            energy += 0.5*p.mass*(p.vx**2 + p.vy**2)
        #print('Conservation of energy: ', energy/energy0)

        # --- Set new time.
        t0 = t

        #print('1: ', particles[0].position(), particles[0].velocity())
        #print('2: ', particles[1].position(), particles[1].velocity())
        #print(particles[0].colcount, particles[1].colcount)

        for i in involved:

            # Calculate if and when the particle(s) involved in the collision 
            # will collide with all other particles, and store the collision times.
            for j in range(0, len(particles)):
                if j != i:
                    dt, v1, v2 = time_until_collision(particles[i], particles[j])
                    if dt < float('inf'):
                        temp = t0 + dt
                        heappush(collisions, (temp, np.array([i, j]), np.array([particles[i].colcount, particles[j].colcount]), np.array([v1, v2])))

            # Calculate when the particle(s) involved in the collision
            # will collide with a wall, and store the collision times.
            dt, v1 = time_until_wall_collision(particles[i])
            if dt < float('inf'):
                temp = t0 + dt
                heappush(collisions, (temp, np.array([i]), np.array([particles[i].colcount]), np.array([v1])))

        # --- Identify the earliest collision.
        t, involved, colcount, new_velocities = heappop(collisions)
        for i in range(len(involved)):
            while colcount[i] != particles[involved[i]].colcount:
                #print('INVALID!')
                t, involved, colcount, new_velocities = heappop(collisions)
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

    #while collisions:
    #    print(heappop(collisions))

if __name__ == "__main__":

    # --- Initialise elastisity constants.
    xi = 1.0

    main(b, True)