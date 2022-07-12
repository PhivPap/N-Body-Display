import random


# CFG
universe_file = "in/a9.tsv"
bodies = 2000
universe = ((-8.888e15, 8.888e15), (-5e15, 5e15))   # m (x, y)
mass_magnitude_range = (10, 34)             # kg   |    translates to [1e+10, 1e+25]
vel_range = (-5e4, 5e4)                     # m / s

# LOGIC
def exporand(mang_low, mang_high):
    magnitude = random.uniform(mang_low, mang_high)
    rgn = random.uniform(1,10)
    rgn *= 10**magnitude
    return rgn

t = '\t'
def new_body(num):
    return {
        "id"    : str(num),
        "mass"  : str(exporand(mass_magnitude_range[0], mass_magnitude_range[1])),
        "x"     : str(random.uniform(universe[0][0], universe[0][1])),
        "y"     : str(random.uniform(universe[1][0], universe[1][1])),
        # "vel_x": "0",
        # "vel_y": "0"
        "vel_x" : str(random.uniform(vel_range[0], vel_range[1])),
        "vel_y" : str(random.uniform(vel_range[0], vel_range[1]))
    }

f = open(universe_file, "w")
f.write("id" + t + "mass" + t + "x" + t + "y" + t + "vel_x" + t + "vel_y" + "\n")

for i in range(bodies):
    b = new_body(i)
    f.write(b["id"] + t + b["mass"] + t + b["x"] + t + b["y"] + t + b["vel_x"] + t + b["vel_y"] + "\n")