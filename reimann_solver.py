import numpy                #numpy is a library for array operations akin to MATLAB
import matplotlib.pyplot as plt    #matplotlib is 2D plotting library

# Reimann solver: different ICs on x<0 and x>0. At time t=0+ the gases interact
# and we can see the resulting conditions. No boundaries so if t is too high all
# you will see is a flat line.

# Start with the function output, input the values for initial conditions and run
# to produce a graph at the desired time. Variables should be non-dimensional

# nx = number of x-steps, increasing this value increases resolution and 
# computation time

# t=time, change this to see a different moment in the expansion

def output(nx=100, t=0.04, dim=2):
    # input desired variables
    gamma = 1.4
    
    input_density_right = 1
    input_density_left = 1

    input_velocity_right = -2
    input_velocity_left = 2

    input_pressure_right = 0.4
    input_pressure_left = 0.4
    
    # Create array of x values and and array for outputs
    x = numpy.linspace(-1,1,nx)
    u = numpy.ones(nx)

    # Prepare the left and right initial conditions vectors
    leftconds = numpy.array([input_density_left, input_velocity_left, input_pressure_left])
    rightconds = numpy.array([input_density_right, input_velocity_right, input_pressure_right])

    # Find the desired output for each x value
    for n in range(nx):
        w = reimann_solver(gammaout, rightconds, leftconds, x[n] / t)
        u[n] = w[dim]
   
    # Plot the graph
    plt.plot(x,u)
    plt.show()


def reimann_solver(gamma,conditions_right,conditions_left,xt):
    #rmannsol function takes primitive variables across a discontinuity to find
    #the distribution of those variables after a time t

    #   conditions=[density,velocity,pressure]
    #   conditions_right=variables on right of discontinuity
    #   conditions_left=variables on left of discontinuity

    # initializing the error
    error = 5

    #Set the accuracy of the P* solution
    acc = 0.001

    #Calculate speed of sound based on initial conditions on both sides
    speed_right=(gamma*conditions_right[2]/conditions_right[0])**0.5
    speed_left=(gamma*conditions_left[2]/conditions_left[0])**0.5

    pressure_guess = conditions_left[2] * conditions_right[0] * speed_right + conditions_right[2] * conditions_left[0] * speed_left+ (conditions_left[1] - conditions_right[1]) * conditions_right[0] * conditions_left[0] * speed_right * speed_left
    #Guess the inital pressure using the acoustic wave approximation

    pressure_guess = pressure_guess/(conditions_right[0] * speed_right + conditions_left[0] * speed_left)

    #In case of negative values set pressure small
    if pressure_guess < 0:
        pressure_guess = 0.001

    #counter
    count = 0

    while error > acc:
        # decide the type of wave (exp fan or shock)
        if pressure_guess <= conditions_right[2]:
            # right-facing expansion fan
            f_press_right, df_press_right = exp_fan(
                                                    conditions_right[0],
                                                    pressure_guess,
                                                    conditions_right[2],
                                                    gamma,
                                                    speed_right
                                                    )
        else:
            # right-facing shock wave
            f_press_right, df_press_right = shock(
                                                    conditions_right[0],
                                                    pressure_guess,
                                                    conditions_right[2],
                                                    gamma
                                                    )
        if pressure_guess <= conditions_left[2]:
            # left facing expansion fan
            f_press_left, df_press_left = exp_fan(
                                                    conditions_left[0],
                                                    pressure_guess,
                                                    conditions_left[2],
                                                    gamma,
                                                    speed_left)
        else:
            f_press_left, df_press_left = shock(
                                                    conditions_left[0],
                                                    pressure_guess,
                                                    conditions_left[2],
                                                    gamma
                                                    )
        
        # Use Newton-Raphson to find P*: There should be a wave of some 
        # kind of wave on both sides and some pressure P* between them) 
        # Using the values of f(P*) and df(P*) calculated above we can
        # find a value for P* and iterate until it stops changing
        # i.e. error < accuracy
        
        pstar = pressure_guess - (
                                    (
                                    f_press_right + f_press_left +
                                    conditions_right[1] - conditions_left[1]
                                    )
                                    
                                    /
                                    
                                    (
                                    df_press_left + df_press_right
                                    )
                                 )

        # Update the error
        error = abs( (pstar - pressure_guess) / pstar)
        # print(error)

        # update pressure_guess
        pressure_guess = pstar

        # Update the counter
        count += 1
        # print(count)

        # Timeout check
        if count >= 200:
            print ('timeout')
            break
    
    # Calculate the conditions on the right/left sides to find the
    # velocity of the contact surface u*
    if pstar <= conditions_right[2]:
        # right-facing expansion fan
        f_press_right, df_press_right = exp_fan(
                                                conditions_right[0],
                                                pstar,
                                                conditions_right[2],
                                                gamma,
                                                speed_right
                                                )
    else:
        # right-facing shock wave
        f_press_right, df_press_right = shock(
                                                conditions_right[0],
                                                pstar,
                                                conditions_right[2],
                                                gamma
                                                )
    if pstar <= conditions_left[2]:
        # left facing expansion fan
        f_press_left, df_press_left = exp_fan(
                                                conditions_left[0],
                                                pstar,
                                                conditions_left[2],
                                                gamma,
                                                speed_left)
    else:
        # left-facing shock wave
        f_press_left, df_press_left = shock(
                                                conditions_left[0],
                                                pstar,
                                                conditions_left[2],
                                                gamma
                                                )    
    # Calculate u*
    ustar =   0.5 * (conditions_left[1] + conditions_right[1])
    + 0.5 * (f_press_right - f_press_left)

    # We find the conditions based on the ratio of x / t for the solution
    if xt >= ustar:
        if pstar <= conditions_right[2]:
            # Region behind right-facing expansion fan
            density, rspeed = kfan(
                                    conditions_right[0],
                                    pstar,
                                    conditions_right[2],
                                    gamma,
                                    speed_right
                                 )
            # Find the x / t of the head and tail of the fan
            head_fan_right = conditions_right[1] + speed_right
            tail_fan_right = ustar + rspeed
            if xt <= head_fan_right:
                if xt <= tail_fan_right:
                    result = numpy.array([density, ustar, pstar])
                else:
                    # In this situation we need the conditions inside 
                    # the expansion fan
                    result = numpy.zeros(3)
                    result[0] = conditions_right[0]*(2 / (gamma + 1) - 
                                (gamma - 1) * 
                                (conditions_right[1] - xt) / 
                                ((gamma + 1) * speed_right))**(2 / (gamma - 1))

                    result[1] = ( 2 / (gamma + 1))*(-speed_right + (gamma - 1)
                                * conditions_right[1] / 2
                                + xt)

                    result[2] = conditions_right[2] * (2 / (gamma + 1) - 
                                (gamma - 1) * (conditions_right[1] - xt) 
                                / ((gamma + 1) * speed_right))**((2 * gamma) / (gamma - 1))

            else:
                result = conditions_right

        else:
            # Find values of shock wave on the right side
            density = kshock(conditions_right[0], pstar, conditions_right[2], gamma)
            shockfront =    conditions_right[1] + shockspeed(pstar, conditions_right[2], speed_right, gamma)
            if xt >= shockfront:
                result = conditions_right
            else:
                result = numpy.array([density, ustar, pstar])

    else:
        # xt is on the left side of the contact surface
        if pstar <= conditions_left[2]:
            # Region behind left-facing expansion fan
            ldensity1, speed = kfan(
                                    conditions_left[0],
                                    pstar,
                                    conditions_left[2],
                                    gamma,
                                    speed_left
                                 )
            # Find the x / t of the head and tail of the fan
            head_fan_left = conditions_left[1] - speed_left
            tail_fan_left = ustar - speed
            if xt >= head_fan_left:
                if xt >= tail_fan_left:
                    result = numpy.array([ldensity1, ustar, pstar])
                
                else:
                    # In this situation we need the conditions inside 
                    # the expansion fan
                    result = numpy.zeros(3)
                    result[0] = conditions_left[0]*(2 / (gamma + 1) + 
                                (gamma - 1) * 
                                (conditions_left[1] - xt) / 
                                ((gamma + 1) * speed_left)) **(2 / (gamma - 1))

                    result[1] = ( 2 / (gamma + 1)) *(speed_left + (gamma - 1)* conditions_left[1] / 2
                                + xt)

                    result[2] = conditions_left[2] * (2 / (gamma + 1) + 
                                (gamma - 1) * (conditions_left[1] - xt) 
                                / ((gamma + 1) * speed_left)) **((2 * gamma) / (gamma - 1))

            else:
                result = conditions_left

        else:
            # Find values of shock wave on the left side
            ldensity2 = kshock(conditions_left[0], pstar, conditions_left[2], gamma)
            shockfront =    conditions_left[1] - shockspeed(pstar, conditions_left[2], speed_left, gamma)
            if xt < shockfront:
                result = conditions_left
            else:
                result = numpy.array([ldensity2, ustar, pstar])
   
    # just in case it falls to zero at any point
    if result[2] == 0:
        result[2] = 0.001
    # And thats all, return whatever output its supposed to
    return result

def shock(rhok,p,pk,gamma):
    # %         calculates f(P*) behind a shock wave
    
    fp = (p - pk) * ((2 / (rhok * (gamma + 1))) / (p + pk * (gamma - 1) / (gamma + 1))) **0.5
    dfp=(1 - (p - pk) / (2 * (p + pk * (gamma - 1) / (gamma + 1)))) * ((2 / (rhok * (gamma + 1))) / (p + pk * (gamma - 1) / (gamma+1))) ** 0.5;
    return fp, dfp

def exp_fan(rhok3,p3,pk4,gamma,ck):
    # %         calculates f(P*) behind expansion fan
    
    fp =    ((p3 / pk4) ** ((gamma - 1) / (2 * gamma)) - 1) * (2 * ck) / (gamma - 1)
    dfp =   (1 / (rhok3 * ck)) * (p3 / pk4) ** (-1 * (gamma + 1) / (2 * gamma))
    return fp, dfp

def kshock(rhok1,pstar,pk1,gamma):
    # %         Calculates the parameters behind a shock
    
    density = rhok1 * ((pstar / pk1) + (gamma - 1) / (gamma + 1)) / ((gamma - 1) * pstar / ((gamma + 1) * pk1) + 1);
    return density

def kfan(rhok2,p1,pk2,gamma,ck1):
    #Calculates the parameters behind an expansion fan
    
    pr = p1 / pk2
    density = rhok2 * pr ** (1 / gamma)
    speed = ck1 * pr ** ((gamma - 1) / (2 * gamma))
    return density, speed

def shockspeed(p2,pk3,ck2,gamma):
    # Calculates the velocity of a shock wave in the fixed reference
    # frame
    
    shockfront = ck2 * ((gamma + 1) * p2 / (2 * gamma * pk3)
                 + (gamma - 1) / (2 * gamma)) ** 0.5
    return shockfront