import numpy as np              #numpy is a library for array operations akin to MATLAB
import matplotlib.pyplot as plt    #matplotlib is 2D plotting library

#  W is primitive variables [Density, Velocity, Pressure]
#  U is conserved variables []
#  F is flux variables []

def muscl(tmax = 5, Xtrue = 100):
    # Start with the number of steps and the
    # minimum and maximum extent of the domain
    
    X = Xtrue + 4
    xmin = -1
    xmax = 1

    # N = 2 * X
    # Set dx
    dx = (xmax - xmin) / (Xtrue)
    t = 0
    n = 2
    w = np.zeros(4, X)
    w[0, : ] = 1

    w[2, 0.1 * Xtrue + 2 : 0.25 * Xtrue + 2] = 2
    w[3, 2 : X] = np.linspace(xmin, xmax, Xtrue)
    w_old = w

    plt.axis([xmin, xmax, -2, 2])
    plt.ion()

    while t < tmax


        # Set dt using CFL condition (max value of u*dt/dx = 1)
        dt = dx / np.amax(w[2, : ])
        
        if dt < 0.000001:
            dt = 0.000001
        
        if (t + dt) > tmax:
            dt = tmax - t
            t = tmax
        else
            t = t + dt

        for i in range(Xtrue)
            ii = i + 2

            # We start with finding predicted values for the
            # primitive variables. Using the function we can
            # get predicted values and limited slopes
            # (important to stop superfluous fluctuations)


            w_pred_dw, slope_dw = pred(
            	                        w_old[ : , ii],
            	                        w_old[ : , ii + 1],
            	                        w_old[ : , ii + 2],
            	                        dt,
            	                        gamma
            	                        )
            
            w_pred, slope = pred(
            	                w_old[ : , ii -1],
            	                w_old[ : , ii],
            	                w_old[ : , ii + 1],
            	                dt,
            	                gamma)
            
            w_pred_uw, slope_uw = pred(
            	                        w_old[ : , ii - 2],
            	                        w_old[ : , ii - 1],
            	                        w_old[ : , ii],
            	                        dt,
            	                        gamma
            	                        )

            # Now using those values to correct predictions

            dw_left =   0.5 *
                        (w_old[0 : 2, ii] + 
                        w_pred[0 : 2]) + 
                        (0.5 * abs(w_old[2, ii + 1] - 
                        w_old[3 , ii])) *
                        slope

            dw_right =   0.5 *
                        (w_old[0 : 2, ii + 1] + 
                        w_pred_dw[0 : 2]) - 
                        (0.5 * abs(w_old[3, ii] - 
                        w_old[3 , ii - 1])) *
                        slope_dw

            uw_left =   0.5 *
                        (w_old[0 : 2, ii - 1] + 
                        w_pred_uw[0 : 2]) + 
                        (0.5 * abs(w_old[3, ii] - 
                        w_old[3 , ii - 1])) *
                        slope_uw

            uw_right =  0.5 *
                        (w_old[0 : 2, ii - 1] + 
                        w_pred_[0 : 2]) - 
                        (0.5 * abs(w_old[3, ii] - 
                        w_old[3 , ii - 1])) *
                        slope

            # Using the reimann solver we can get the conditions
            # between the predicted values

            notflux_uw = rs.reimann_solver(uw_left, uw_right, gamma, 0)
            notflux_dw = rs.reimann_solver(dw_left, dw_right, gamma, 0)
             
            # Convert primitive variables into flux variables
            flux_dw = w2f(notflux_dw)
            flux_uw = w2f(notflux_uw)

            u_uw = w2u(w_old[0 : 2, ii], gamma)
            u_current = u_uw - (dt / dx) * (flux_dw - flux_uw)

            doubleu = u2w(u_current, g)
            w[0 : 2, ii] = doubleu[0 : 2]

        w[0 : 2, 0 ] = w[0 : 2, 3]
        w[0 : 2, 1 ] = w[0 : 2, 2]
        w[1, 0 : 1] = -w[1, 0 : 1]
        w[0 : 2, X - 1] = w[0 : 2, X - 4]
        w[0 : 2, X - 2] = w[0 : 2, X - 3]
        w[1, X - 2 : X - 1] = -w[1, X - 2 : X - 1]

        w_old = w
        plt.plot(w[3, : ], w[2, : ])
        plt.pause(0.05)

def w2u(win,g):
    uout = np.zeros(3)
    uout[0] = win[0]
    uout[1] = win[1] * win[1]
    uout[2] = win[2] / (g - 1) + 0.5 * win[2] * win[1] ** 2
    
def u2w(uin,g):
    wout[0] = uin[0]
    wout[1] = uin[1] / uin[0]
    wout[2] = (g - 1) * (uin[2] - 0.5 * wout[0] * wout[1] ** 2)

def w2f(win,g):
    fout[0] = win[0] * win[1]
    fout[1] = win[0] * win[1] ** 2 + win[2]
    fout[2] = ((win[2] / (g - 1) + 0.5 * win[0] * win[1] ** 2) + win[2]) * win[1]

def u2f(uin,g):
    fout[0] = uin[1]
    fout[1] = 0.5 * (3 - g) * (uin[1] ** 2) / uin[0] + (g - 1) * uin[2]
    fout[2] = 0.5 * (1 - g) * (uin[1] **3)/(uin[0] ** 2) + g * uin[2] * uin[1] / uin[0]

def pred(w0,w1,w2,dt,g):
    # %         choose slope limiter beta is superbee for b=2, minmod for
    # %         b=1,intermediate dissipation for 1<b<2
    dw = beta(w1 - w0, w2 - w1, 1.9);
    w2l = w1[0 : 2] + 0.5 * dw * abs(w2[3] - w1[3])
    w0r = w1[0 : 2] - 0.5 * dw * abs(w1[3] - w0[3])
    fpl = w2f(w2l, g)
    fpr = w2f(w0r, g)
    fpr[3] = w1[3] - 0.5 * abs(w1[3] - w0[3])
    fpl[3] = w1[3] + 0.5 * abs(w2[3] - w1[3])
    
    up = w2u(w1[0 : 2], g) - (dt / abs(fpl[3] - fpr[3])) * (fpl[0 : 2] - fpr[0 : 2])
    
    dw=dw
    wp = u2w(up, g)
    wp[3] = w1[3]