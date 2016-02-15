import numpy as np
import matplotlib.pyplot as plt

# --- Load data
minloc  = np.loadtxt('minloc.dat', delimiter='\n')

# --- Get number of steps
N       = minloc.size

# --- Ask user to give number of species
prompt = 'What is the number of species? (note use a even number)';
nSpecies = input(prompt);

print nSpecies, N

# --- Create an array of time step (1,2,3,...,N)
time    = np.arange(0, N)
test    = np.arange(0,3)
print len(minloc)
print len(time)

# --- Plot time vs location of minimum
# --- comment/uncomment the one you wish
# --- With a linear scale
plt.plot(minloc, time)
plt.xlim(0, nSpecies)
plt.ylim(0, N)
plt.xlabel('Mutation activity at species index')
plt.ylabel('Time (steps)')
plt.show()
plt.savefig('evolution.png')

# --- With a log scale on the y-axis
plt.semilogy(minloc, time)
plt.xlim(0, nSpecies)
plt.ylim(0, N)
plt.xlabel('Mutation activity at species index')
plt.ylabel('Time (steps)')
plt.show()
plt.savefig('evolution_logy.png')


# --- Create an array of distances between to successive minimum location
#distance = zeros(N-1,1);

# --- Loop over steps
#for i = 1:N-1
   #d = abs(minloc(i+1)-minloc(i));
   # --- Remember periodic boundary conditions
   #if (d > Nspecies/2)
    #  d = Nspecies - d;
   #end
   #distance(i) = d;
#end

# --- Construct histogram of distances
#     We disregard the first part of the evolution
#     so that we reach steady state (this is kind of arbitrary here)
#Nstart = floor(N/10);
#frequency = zeros(Nspecies/2+1,1);
#for i = Nstart:N-1
#   d = distance(i);
#   frequency(d+1) = frequency(d+1) + 1;
#end

# --- Normalization
#frequency = frequency./(N-Nstart);

# --- Plot histogram of distances between consecutives minima.
#figure(2)
#loglog(0:Nspecies/2,frequency,'color','k')
#axis([0 Nspecies/2 1E-10 1])
#xlabel('Distance between consecutive minima.')
#ylabel('Frequency of occurence.')