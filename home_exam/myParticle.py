import numpy as np

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
        self.last        = -1.0
        Particle.particle_count += 1

    def __init__(self, x, y, vx, vy, r, m):
        self.x, self.y   = x, y
        self.vx, self.vy = vx, vy
        self.radius      = r
        self.mass        = m
        self.colcount    = 0
        self.id          = Particle.particle_count
        self.last        = -1.0
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