import genSemigroups as gs

for gamma in range(2,12,2):
    for N in range(3*gamma+2, 3*gamma+12,2):
        for d in range(2*gamma+3, 2*gamma+13,2):
            sg = gs.staircaseSG(N,gamma,d)
            sg.showStats()
