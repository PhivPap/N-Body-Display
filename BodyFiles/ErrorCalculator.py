import csv

# CFG
reference_filepath = "out/b.tsv"
approx_filepath = "out/a.tsv"

# LOGIC
def parse_bodies(filepath):
    bodies = []
    with open(filepath) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        next(rd) # skips tsv header row
        for row in rd:
            body = [float(attr) for attr in row[1:6]]
            bodies.append(body)
    return bodies


ref_bodies = parse_bodies(reference_filepath)
approx_bodies = parse_bodies(approx_filepath)
ref_b_count = len(ref_bodies)
apx_b_count = len(approx_bodies)

if ref_b_count != apx_b_count:
    print("rip ->", ref_b_count, "ref bodies VS", apx_b_count, "approx bodies.")
    exit(1)

invalid_bodies = 0
max_x_error_body = max_y_error_body = max_vel_y_error_body = max_vel_x_error_body = -1
max_x_error_prcnt = max_y_error_prcnt = max_vel_y_error_prcnt = max_vel_x_error_prcnt = 0
avg_x_error_prcnt = avg_y_error_prcnt = avg_vel_y_error_prcnt = avg_vel_x_error_prcnt = 0

# max_p_error_prcnt = max_v_error_prcnt = 0
# avg_p_error_prcnt = avg_v_error_prcnt = 0

for i in range(ref_b_count):
    ref_body = ref_bodies[i]
    approx_body = approx_bodies[i]

    # Take a look below. vvvvvvvvvvvvvvvv
    # ref_body[0] - approx_body[0] - mass
    # ref_body[1] - approx_body[1] - x
    # ref_body[2] - approx_body[2] - y
    # ref_body[3] - approx_body[3] - vel_x
    # ref_body[4] - approx_body[4] - vel_y
    # ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    if(ref_body[0] != approx_body[0]):  # if masses are not the same: problem :)
        print("Body mass does not match. Bye!")
        exit(1)
    
    x_diff = ref_body[1] - approx_body[1]
    y_diff = ref_body[2] - approx_body[2]
    vel_x_diff = ref_body[3] - approx_body[3]
    vel_y_diff = ref_body[4] - approx_body[4]

    x_error_prcnt = abs(x_diff / ref_body[1]) * 100
    y_error_prcnt = abs(y_diff / ref_body[2]) * 100
    vel_x_error_prcnt = abs(vel_x_diff /  ref_body[3]) * 100
    vel_y_error_prcnt = abs(vel_y_diff /  ref_body[4]) * 100

    if x_error_prcnt > max_x_error_prcnt:
        max_x_error_prcnt = x_error_prcnt
        max_x_error_body = i

    if y_error_prcnt > max_y_error_prcnt:
        max_y_error_prcnt = y_error_prcnt
        max_y_error_body = i

    if vel_x_error_prcnt > max_vel_x_error_prcnt:
        max_vel_x_error_prcnt = vel_x_error_prcnt
        max_vel_x_error_body = i

    if vel_y_error_prcnt > max_vel_y_error_prcnt:
        max_vel_y_error_prcnt = vel_y_error_prcnt
        max_vel_y_error_body = i

    if x_error_prcnt > 1 or y_error_prcnt > 1 or vel_x_error_prcnt > 1 or vel_y_error_prcnt > 1:
        invalid_bodies += 1
    
    avg_x_error_prcnt += x_error_prcnt
    avg_y_error_prcnt += y_error_prcnt
    avg_vel_x_error_prcnt += vel_x_error_prcnt
    avg_vel_y_error_prcnt += vel_y_error_prcnt


avg_x_error_prcnt = round(avg_x_error_prcnt / ref_b_count, 2)
avg_y_error_prcnt =     round(avg_y_error_prcnt / ref_b_count, 2)
avg_vel_x_error_prcnt = round(avg_vel_x_error_prcnt / ref_b_count, 2)
avg_vel_y_error_prcnt = round(avg_vel_y_error_prcnt / ref_b_count, 2)

max_x_error_prcnt  = round(max_x_error_prcnt, 2) 
max_y_error_prcnt = round(max_y_error_prcnt, 2)
max_vel_x_error_prcnt = round(max_vel_x_error_prcnt, 2)
max_vel_y_error_prcnt = round(max_vel_y_error_prcnt, 2)


print("\nInvalid bodies (error > 1%):", invalid_bodies, "\n")
print("AVG ERROR (%):\n  Position: [x] ", avg_x_error_prcnt, "%, [y] ", avg_y_error_prcnt, "%\n  Velocity: [x] ", avg_vel_x_error_prcnt, "%, [y] ", avg_vel_y_error_prcnt, "%\n", sep="")
print("MAX ERROR (%):\n  Position: [x] ", max_x_error_prcnt, "% (body: ",  max_x_error_body, "), [y] " ,max_y_error_prcnt, "% (body: ", max_y_error_body, ")", sep="")
print("  Velocity: [x] ", max_vel_x_error_prcnt, "% (body: ",max_vel_x_error_body, "), [y] ", max_vel_y_error_prcnt, "% (body: ", max_vel_y_error_body, ")\n", sep="")