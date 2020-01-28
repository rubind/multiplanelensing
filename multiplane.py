import numpy as np
import decimal
decimal.getcontext().prec = 50
D = decimal.Decimal

def mydot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def myadd(a, b):
    return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]

def mysub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]

def myscalarmult(a, b):
    return [a[0]*b, a[1]*b, a[2]*b]

def myscalardiv(a, b):
    return [a[0]/b, a[1]/b, a[2]/b]

def myarray(vals):
    return [D(item) for item in vals]

def mypathintegral(xyz_start, xyz_fin, steps, masses_xyz, masses):
    assert steps % 2 == 1
    
    tvals = [D(item) for item in np.linspace(0, 1, steps)]
    
    xyzs = []
    for tval in tvals:
        xyzs.append(myadd(myscalarmult(xyz_start, D(tval)), myscalarmult(xyz_fin, D(1 - tval))))

    total_integral = D(0.)
    dxyz3 = mydot(mysub(xyzs[1], xyzs[0]), mysub(xyzs[1], xyzs[0])).sqrt() / D(3)

    for i in range(steps):
        if i == 0 or i == steps - 1:
            simpson_weight = D(1)
        else:
            simpson_weight = D(2) + D(2*(i % 2))
        

        for j in range(len(masses)):
            total_integral += masses[j] * D(9.5738379e-8) * simpson_weight/ mydot(mysub(xyzs[i], masses_xyz[j]),
                                                                                  mysub(xyzs[i], masses_xyz[j])).sqrt()
    total_integral *= dxyz3
    
    return total_integral


def get_travel_time(starting_xyz, unit_vector, target_xyz):
    # Solve[D[(t1 - (s1 + u1*t))^2 + (t2 - (s2 + u2*t))^2 + (t3 - (s3 + u3*t))^2, t] == 0, t]
    travel_time = mydot(mysub(target_xyz, starting_xyz), unit_vector)
    travel_time = travel_time/mydot(unit_vector, unit_vector)

    print("travel_time", travel_time)
    assert travel_time > 0
    return travel_time

def get_point_of_closest_approach(starting_xyz, unit_vector, target_xyz):
    travel_time = get_travel_time(starting_xyz, unit_vector, target_xyz)
    
    return myadd(starting_xyz, myscalarmult(unit_vector, travel_time))

def get_new_unit_vector(unit_vector, mass_xyz, mass, point_of_closest_approach):
    diff_vector = mysub(mass_xyz, point_of_closest_approach)
    impact_squared = mydot(diff_vector, diff_vector)
    impact = impact_squared.sqrt()
    
    defl_angle = D(1.91476758e-7) * mass/impact

    new_unit = myadd(unit_vector, myscalarmult(diff_vector, defl_angle/impact))
    new_unit = myscalardiv(new_unit, mydot(new_unit, new_unit).sqrt())

    return new_unit


def get_diff_time(masses_xyz, masses, source_xyz_init, Shapiro_steps = 1001):
    """Observer at origin, distances in Mpc, masses in 1e12 solar masses
    masses_xyz will be sorted by descending LoS distance
    """

    travel_times = [get_travel_time(myarray(source_xyz_init), myscalarmult(myarray(source_xyz_init), D(-1)), myarray(item)) for item in masses_xyz]
    print("travel times, neglecting lensing:", travel_times)
    inds = np.argsort(travel_times)
    print("inds for sorting", inds)

    
    masses_xyz = [myarray(masses_xyz[ind]) for ind in inds]
    masses = myarray([masses[ind] for ind in inds])
    source_xyz_init = myarray(source_xyz_init)

    source_unit_vector_init = myscalarmult(source_xyz_init, D(-1))

    for iteration in range(6):
        print("\n\n\nInteration %i\n\n\n" % iteration)
        source_xyz = [myscalarmult(source_xyz_init, D(1))]

        source_unit_vector = myscalarmult(source_unit_vector_init, D(1))
        source_unit_vector = myscalardiv(source_unit_vector, mydot(source_unit_vector, source_unit_vector).sqrt())

        for i in range(len(masses)):
            print("source_unit_vector", source_unit_vector)
            point_of_closest_approach = get_point_of_closest_approach(source_xyz[i], source_unit_vector, masses_xyz[i])
            print(point_of_closest_approach, point_of_closest_approach)

            source_unit_vector = get_new_unit_vector(source_unit_vector, masses_xyz[i], masses[i], point_of_closest_approach)
            print("source_unit_vector", source_unit_vector)

            source_xyz.append(point_of_closest_approach)


        point_of_closest_approach = get_point_of_closest_approach(source_xyz[-1], source_unit_vector, myarray([0, 0, 0]))
        print("point_of_closest_approach to observer", point_of_closest_approach)
        source_xyz.append(point_of_closest_approach)

        for this_xyz in source_xyz:
            print("source_xyz", this_xyz)
        source_unit_vector_init = myadd(source_unit_vector_init, myscalarmult(point_of_closest_approach, -1))

    assert len(source_xyz) - 2 == len(masses)


    direct_path = mydot(mysub(source_xyz[0],  source_xyz[-1]), mysub(source_xyz[0], source_xyz[-1])).sqrt()
    total_shapiro_direct = mypathintegral(source_xyz[0], source_xyz[-1], steps = Shapiro_steps, masses_xyz = masses_xyz, masses = masses)
    
    total_path = D(0)
    total_shapiro_indirect = D(0)
    
    for i in range(len(source_xyz) - 1):
        total_path += mydot(mysub(source_xyz[i], source_xyz[i+1]), mysub(source_xyz[i], source_xyz[i+1])).sqrt()
        total_shapiro_indirect += mypathintegral(source_xyz[i], source_xyz[i+1], steps = Shapiro_steps, masses_xyz = masses_xyz, masses = masses)
        
    print(direct_path, total_path)

    diff_path = total_path - direct_path
    print(diff_path)

    diff_time = D(1.02927125e14) * diff_path

    print(diff_time)
    shapiro_delay = D(1.02927125e14)*(total_shapiro_indirect - total_shapiro_direct)
    
    print("Shapiro", shapiro_delay, "masses", masses, "just path length", diff_time, "ratio", shapiro_delay/diff_time, "abs delay", D(1.02927125e14)*total_shapiro_direct)
    
    return float(diff_time), float(shapiro_delay)

def one_mass_delay(mass, impact, dLS, dL):
    return 0.5*(1.91476758e-7*mass/impact)**2. * 1.02927125e14 * dLS*dL/(dLS + dL)

if __name__== "__main__":
    print("Running tests!")

    diff_time, NA = get_diff_time(masses_xyz = [[10, 10, 0]], masses = [1], source_xyz_init = [20., 0., 0])
    diff_time, NA = get_diff_time(masses_xyz = [[10, 100, 0]], masses = [1], source_xyz_init = [20., 0., 0])

    
    diff_time, NA = get_diff_time(masses_xyz = [[10, 0.1, 0]], masses = [1], source_xyz_init = [20., 0., 0])
    assert abs(diff_time - one_mass_delay(1., 0.1, 10., 10.)) < 1
    
    diff_time, NA = get_diff_time(masses_xyz = [[10, 1., 0]], masses = [1], source_xyz_init = [20., 0., 0])
    assert abs(diff_time - one_mass_delay(1., 1., 10., 10.)) < 0.01

    diff_time, NA = get_diff_time(masses_xyz = [[10, 0, 1.]], masses = [1], source_xyz_init = [20., 0., 0])
    assert abs(diff_time - one_mass_delay(1., 1., 10., 10.)) < 0.01

    diff_time, NA = get_diff_time(masses_xyz = [[0, 10, 1.]], masses = [1], source_xyz_init = [0., 20., 0])
    assert abs(diff_time - one_mass_delay(1., 1., 10., 10.)) < 0.01

    diff_time, NA = get_diff_time(masses_xyz = [[1., 1., 10.]], masses = [1], source_xyz_init = [0., 0., 20.])
    assert abs(diff_time - one_mass_delay(1., 1.41421356237, 10., 10.)) < 0.01

    diff_time, NA = get_diff_time(masses_xyz = [[10, 1., 0]], masses = [2], source_xyz_init = [20., 0., 0])
    assert abs(diff_time - one_mass_delay(2., 1., 10., 10.)) < 0.01
    
    diff_time, NA = get_diff_time(masses_xyz = [[10, 1., 0]], masses = [2], source_xyz_init = [40., 0., 0])
    assert abs(diff_time - one_mass_delay(2., 1., 30., 10.)) < 0.01

    get_diff_time(masses_xyz = [[1, 1., 0], [10, 1., 0]], masses = [2, 1], source_xyz_init = [40., 0., 0])
    get_diff_time(masses_xyz = [[10, 1., 0], [1, 1., 0]], masses = [1, 2], source_xyz_init = [40., 0., 0])
