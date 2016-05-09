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

class Particle():

    particle_count = 0

    def __init__(self):
        self.x, self.y   = 0.5, 0.5
        self.vx, self.vy = 1.0, 0.0
        self.radius      = 1.e-1
        self.mass        = 1.0
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
        return self.colcount-1

    # --- Settere.
    def set_velocity(self, v):
        self.vx = v[0]
        self.vy = v[1]

    # --- Functions.
    def update(self, dt):
        self.x = self.x + dt*self.vx
        self.y = self.y + dt*self.vy        
   
# --- Calculate time until collision with a wall.
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
        return timex, np.array([-xi_wall*particle.vx, xi_wall*particle.vy])
    else:
        return timey, np.array([xi_wall*particle.vx, -xi_wall*particle.vy])

# --- Calculate time until collision between particles.
#     input : two Particles
#     output: time until collision between particles, corresponding new velocity after collision.
def time_until_collision(particle1, particle2):
    dx = np.array([particle2.x - particle1.x, particle2.y - particle1.y])     # particle1.position() - particle2.position()
    dv = np.array([particle2.vx - particle1.vx, particle2.vy - particle1.vy]) # particle1.velocity() - particle2.velocity()
    R  = particle1.radius + particle2.radius
    d  = (np.dot(dv,dx))**2 - (np.dot(dv,dv))**2 * ((np.dot(dx,dx))**2 - R**2)
    if np.dot(dv,dx) >= 0:
        return float('inf')
    elif d <= 0:
        return float('inf')
    else:
        velocity1 = particle1.velocity() + ((1+xi)*particle2.mass()/(particle1.mass()+particle2.mass())*np.dot(dv,dx)/R**2)*dx
        velocity2 = particle2.velocity() - ((1+xi)*particle1.mass()/(particle1.mass()+particle2.mass())*np.dot(dv,dx)/R**2)*dx
        return -(np.dot(dv,dx) + np.sqrt(d)) / (np.dot(dv,dv)), velocity1, velocity2

def main():

    # --- Initialise square box.
    xmin = 0.0
    xmax = 1.0
    ymin = 0.0
    ymax = 1.0

    # --- Initialisation of particles and collision priority heap.
    particles = [Particle()]
    collisions = []
    for i in range (len(particles)):

        # Calculate if and when particle will collide with all other particles,
        # and store the collision times.
        for j in range(i+1, len(particles)):
            t, v1, v2 = time_until_collision(particles[i], particles[j])
            if t < float('inf'):
                heappush(collisions, (t, np.array([i, j]), np.array([particles[i].increment_colcount(), particles[j].increment_colcount()]), np.array([v1, v2])))

        # Calculate when particle will collide with a wall,
        # and store the collision times.
        t, v1 = time_until_wall_collision(particles[i])
        heappush(collisions, (t, np.array([i]), np.array([particles[i].increment_colcount()]), np.array([v1])))

    # --- Identify the earliest collision.
    dt, involved, coll_count, new_velocities = heappop(collisions)
    print(dt)
    print(involved)
    print(coll_count)
    print(new_velocities)
    
    # --- Loop:
    #while(collisions):
    for i in range(5):
        # --- Move all particles forward in time (straight lines, constant 
        #     velocity) until the earliest collision.
        for i in range (len(particles)):
            particles[i].update(dt)

        # --- For the particle(s) involved in the collision, calculate new 
        #     velocities.
        for i in range(len(involved)):
            particles[involved[i]].set_velocity(new_velocities[i])

        for i in involved:

            # Calculate if and when the particle(s) involved in the collision 
            # will collide with all other particles, and store the collision times.
            for j in range(0, len(particles)):
                if j != i:
                    t, v1, v2 = time_until_collision(particles[i], particles[j])
                    if t < float('inf'):
                        heappush(collisions, (t, np.array([i, j]), np.array([particles[i].increment_colcount(), particles[j].increment_colcount()]), np.array([v1, v2])))

            # Calculate when the particle(s) involved in the collision
            # will collide with a wall, and store the collision times.
            t, v1 = time_until_wall_collision(particles[i])
            heappush(collisions, (t, np.array([i]), np.array([particles[i].increment_colcount()]), np.array([v1])))

        # --- Identify the earliest collision.
        dt, involved, coll_count, new_velocities = heappop(collisions)
        print(dt)
        print(involved)
        print(coll_count)
        print(new_velocities)

        # --- Having resolved the collision, identify the new earliest collision,
        #     and check if it is still valid (if the particle(s) involved have 
        #     collided with something since the collision was calculated, it is 
        #     no longer valid). If the collision is invalid, discard it and move
        #     to the next earliest. Once a valid collision is found, repeat the 
        #     steps in the loop.



    while collisions:
        print(heappop(collisions))

if __name__ == "__main__":

    # --- Initialise elastisity constants.
    xi_wall = 1.0
    xi = 1.0

    main()