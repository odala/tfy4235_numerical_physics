# --- Class: Collision
class Collision():
    def __init__(self, time, involved, colcount):  
        self.time = time
        self.involved = involved
        self.colcount = colcount

    def __cmp__(self, other):
        return -1.0*(self.time < other.time) + 1.0*(self.time >= other.time)