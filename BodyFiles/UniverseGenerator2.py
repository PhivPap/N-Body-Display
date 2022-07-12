import random
import math

# CFG

G = 6.67e-11
t = '\t'
universe_file = "in/a13.tsv"
orbit_variance = 0.1
# mmr: mass mang range
# ormr: orbit range magnitude range
black_holes = {
    "count": 1,
    "mmr": (36, 36),
    "pos_range": (-5e17, 5e17),
    "vel_range": (-1e6, 1e6),
    "list": []
}
stars = {
    "count": 200,
    "mmr": (30, 33),
    "pos_range": (-1e15, 1e15),
    "vel_range": (-1e4, 1e4),
    "ormr": (15.5, 16.5),
    "list": []
}
planets = {
    "count": 400,
    "mmr": (21, 28),
    "pos_range": (-1e15, 1e15),
    "vel_range": (-1e4, 1e4),
    "ormr": (10.7, 13),
    "list": []
}
satelites = {
    "count": 600,
    "mmr": (18, 20),
    "pos_range": (-1e15, 1e15),
    "vel_range": (-1e3, 1e3),
    "ormr": (8.5, 9),
    "list": []
}

# LOGIC
def exporand(magn_range):
    mang_low, mang_high = magn_range
    magnitude = random.uniform(mang_low, mang_high)
    rgn = random.uniform(1,10)
    rgn *= 10**magnitude
    return rgn

def new_body_in_orbit(id, stats, parent):
    orbit_radius = exporand(stats["ormr"])
    mass = exporand(stats["mmr"])
    dx = random.uniform(-orbit_radius, orbit_radius)
    dy2 = abs(dx**2 - orbit_radius**2)
    dy = math.sqrt(dy2)
    if bool(random.getrandbits(1)):
        dy = -dy
    vel = math.sqrt(G * (mass + parent["mass"]) / orbit_radius)
    vel += vel * random.uniform(-orbit_variance, orbit_variance)
    theta = (parent["x"] - dx) / (parent["y"] - dy)
    # print(f"theta = {theta}")
    angle = math.atan(theta)
    # print(f"angle = {angle}")
    dvx = math.cos(-angle) * vel
    dvy = math.sin(-angle) * vel
    # print(dvx, " :: ", dvy)
    if dx > 0:
        dvy = -abs(dvy)
    else:
        dvy = abs(dvy)
    if dy > 0:
        dvx = abs(dvx)
    else:
        dvx = -abs(dvx)
    return {
        "id": id,
        "mass": mass,
        "x": parent["x"] + dx,
        "y": parent["y"] + dy,
        "vel_x": parent["vel_x"] + dvx,
        "vel_y": parent["vel_y"] + dvy
    }

def new_body(id, stats):
    return {
        "id": id,
        "mass": exporand(stats["mmr"]),
        "x": random.uniform(stats["pos_range"][0], stats["pos_range"][1]),
        "y": random.uniform(stats["pos_range"][0], stats["pos_range"][1]),
        "vel_x" : random.uniform(stats["vel_range"][0], stats["vel_range"][1]),
        "vel_y": random.uniform(stats["vel_range"][0], stats["vel_range"][1])
    }


f = open(universe_file, "w")
f.write("id" + t + "mass" + t + "x" + t + "y" + t + "vel_x" + t + "vel_y" + "\n")

counter = 0

if black_holes["count"] == 1:
    black_holes["list"].append({ "id": 0, "mass": exporand(black_holes["mmr"]), "x": 0, "y": 0, "vel_x": 0, "vel_y": 0 })
    counter += 1

else:
    for _ in range(black_holes["count"]):
        black_holes["list"].append(new_body(counter, black_holes))
        counter += 1

if stars["count"] == 1 and black_holes["count"] == 0:
    stars["list"].append({ "id": 0, "mass": exporand(stars["mmr"]), "x": 0, "y": 0, "vel_x": 0, "vel_y": 0 })
else:
    for _ in range(stars["count"]):
        stars["list"].append(new_body_in_orbit(counter, stars, random.choice(black_holes["list"])) if black_holes["count"] > 0 else new_body(counter, stars))
        counter += 1

for _ in range(planets["count"]):
    planets["list"].append(new_body_in_orbit(counter, planets, random.choice(stars["list"])) if stars["count"] > 0 else new_body(counter, planets))
    counter += 1

for _ in range(satelites["count"]):
    satelites["list"].append(new_body_in_orbit(counter, satelites, random.choice(planets["list"])) if planets["count"] > 0 else new_body(counter, satelites))
    counter += 1

bodies = black_holes["list"] + stars["list"] + planets["list"] + satelites["list"]
for b in bodies:
    f.write(str(b["id"]) + t + str(b["mass"]) + t + str(b["x"]) + t + str(b["y"]) + t + str(b["vel_x"]) + t + str(b["vel_y"]) + "\n")

f.close()



# for i in range(bodies):
#     b = new_body(i)
#     f.write(b["id"] + t + b["mass"] + t + b["x"] + t + b["y"] + t + b["vel_x"] + t + b["vel_y"] + "\n")


# def new_body(num):
#     return {
#         "id"    : str(num),
#         "mass"  : str(exporand(mass_magnitude_range)),
#         "x"     : str(random.uniform(universe[0][0], universe[0][1])),
#         "y"     : str(random.uniform(universe[1][0], universe[1][1])),
#         # "vel_x": "0",
#         # "vel_y": "0"
#         "vel_x" : str(random.uniform(vel_range[0], vel_range[1])),
#         "vel_y" : str(random.uniform(vel_range[0], vel_range[1]))
#     }