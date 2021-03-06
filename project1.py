# author: Nicholas Prayogo
# simulating crank rocker mechanism with coupler point to produce a path similar to a figure 8
# Figure 8 theory: Lemniscate of Gerono https://www.mathcurve.com/courbes2d.gb/gerono/gerono.shtml

import numpy as np
import matplotlib.pyplot as plt
from math import *
from sklearn.metrics import mean_squared_error, max_error
import pandas as pd
import time
from tqdm import tqdm
import warnings
from scipy import integrate
from scipy.optimize import minimize
warnings.filterwarnings("ignore")

def calculate_coupler_position(r1, r2, r3, r4, theta2, r6, r7):
    h1 = r1 / r2
    h3 = r1 / r4
    h5 = (r1**2 + r2**2 - r3**2 + r4**2) / (2 * r2 * r4)
    h2 = r1 / r3
    h4 = (-r1**2 - r2**2 - r3**3 + r4**2) / (2 * r2 * r3)
    b = -2 * sin(theta2)
    d = -h1 + (1 - h3) * cos(theta2) + h5
    e = h1 - (1 + h3) * cos(theta2) + h5
    a = -h1 + (1 + h2) * cos(theta2) + h4
    c = h1 + (1 - h2) * cos(theta2) + h4

    # 2 possible vals of theta3, take 1
    # print((b**2 - 4*a*c))
    if (b**2 - 4*a*c)>=0:
        theta3 = 2 * atan((-b - sqrt(b**2 - 4 * a * c)) / (2 * d))
        theta3_2 = 2 * atan((-b + sqrt(b**2 - 4 * a * c)) / (2 * d))
        # print(theta3)
        if theta3 == None or theta3_2 == None:
            trajectory_possible = False
            return None, None, None, None,  trajectory_possible
            # return None, None,  trajectory_possible

    else:
        trajectory_possible = False
        return None, None,None, None, trajectory_possible
        # return None, None,  trajectory_possible

    r8 = r2 * cos(theta2) + r6 * cos(theta3) - r7 * sin(theta3)
    r9 = r2 * sin(theta2) + r6 * sin(theta3) + r7 * cos(theta3)
    # print(f"theta3: {theta3}, theta32: {theta3_2}")
    r8_2 = r2 * cos(theta2) + r6 * cos(theta3_2) - r7 * sin(theta3_2)
    r9_2 = r2 * sin(theta2) + r6 * sin(theta3_2) + r7 * cos(theta3_2)

    # r_coupler = sqrt(r8**2 + r9**2)
    # theta_coupler = atan(r9 / r8)

    # x_coupler = r8
    # y_coupler = r9

    trajectory_possible = True
    return r8, r9, r8_2, r9_2,  trajectory_possible
    # return r8, r9, trajectory_possible

def constraint1():
    # h1 = r1 / r2
    # h3 = r1 / r4
    # h5 = (r1**2 + r2**2 - r3**2 + r4**2) / (2 * r2 * r4)
    # h2 = r1 / r3
    # h4 = (-r1**2 - r2**2 - r3**3 + r4**2) / (2 * r2 * r3)
    # b = -2 * sin(theta2)
    # d = -h1 + (1 - h3) * cos(theta2) + h5
    # e = h1 - (1 + h3) * cos(theta2) + h5
    # a = -h1 + (1 + h2) * cos(theta2) + h4
    # c = h1 + (1 - h2) * cos(theta2) + h4
    # (b**2 - 4*a*c)
    return

def calculate_error(r1, r2, r3, r4, theta2):
    r7 = 0
    r6 = r3/2
    h1 = r1 / r2
    h3 = r1 / r4
    h5 = (r1**2 + r2**2 - r3**2 + r4**2) / (2 * r2 * r4)
    h2 = r1 / r3
    h4 = (-r1**2 - r2**2 - r3**3 + r4**2) / (2 * r2 * r3)
    b = -2 * sin(theta2)
    d = -h1 + (1 - h3) * cos(theta2) + h5
    e = h1 - (1 + h3) * cos(theta2) + h5
    a = -h1 + (1 + h2) * cos(theta2) + h4
    c = h1 + (1 - h2) * cos(theta2) + h4
    theta3 = 2 * atan((-b - sqrt(b**2 - 4 * a * c)) / (2 * d))
    x_coupler = r2 * cos(theta2) + r6 * cos(theta3) - r7 * sin(theta3)
    y_coupler = r2 * sin(theta2) + r6 * sin(theta3) + r7 * cos(theta3)
    x_target = x_coupler
    y_target = sqrt((x_target**2 - x_target**4))
    error = abs(y_target-y_coupler)**2
    return error

def objective_function(params):
    r1, r2, r3, r4 = params
    try:
        area_error = integrate.quad(calculate_error, 0, np.pi/2, args= (r1,r2,r3,r4,)) # area from 0 to pi/2
    except:
        return 1000000000
    return area_error

if __name__ == "__main__":
    # input link theta
    theta2_range = np.linspace(0, 2 * np.pi, 100)

    # unit in inches
    # s + l < p + q
    # for crank rocker, r2 has to be s
    # let r1 be p, r4 be q, and r3 be l
    resolution = 1000
    min_val = 1.1
    max_val = 10
    # r2_candidates = np.linspace(min_val,max_val, 2)
    r2_candidates = [1.0]
    # r2_candidates = [1.0, 5.0]
    # r3_candidates = np.linspace(4,10,  resolution)
    #
    # r4_candidates = np.linspace(10,15 ,resolution)
    # r1_candidates = np.linspace(3,10,  resolution)
    # r7_candidates = np.linspace(4,10, resolution)

    # r2_candidates = np.linspace(min_val,max_val, resolution)
    # r3_candidates = np.linspace(min_val,max_val,  resolution)
    # r4_candidates = np.linspace(min_val,max_val,  resolution)
    # r1_candidates = np.linspace(min_val,max_val,  resolution)
    # r6_candidates = np.linspace(min_val,max_val,  resolution)
    # r7_candidates = np.linspace(min_val,max_val, 2)

    # 0.3333222225925802, 1e-05, 0.3333222225925802, 0.3333222225925802
 #    solval = [0.2568888888888889,
 # 0.00001,
 # 0.2568888888888889,
 # 0.22033333333333333,
 # 0.03755555555555556]
    # r1_candidates = [solval[0]]
    # r2_candidates = [solval[1]]
    # r3_candidates = [solval[2]]
    # r4_candidates = [solval[3]]
    # r7_candidates = [solval[4]]

    r3_candidates = [1.599]
    r1_candidates = [1.59]
    r4_candidates = [1.554]
    r7_candidates = [1.1]
    experiment_results = []
    a_given = 4
    best_rmse = 99999999
    start = time.process_time()

    new_method = False

    if new_method:
        optimal_result = minimize(objective_function, [1, 1, 1, 1], method="BFGS")
        print(optimal_result)
        print("Time taken: ", time.process_time() - start)

    else:
        counter = 0
        for r2 in r2_candidates:
            # print("r2:", r2)
            for r3 in r3_candidates:
                for r4 in r4_candidates:
                    for r1 in r1_candidates:
                        for r7 in r7_candidates:
                            # only simulate for crank rocker dimensions
                            if ((r2 + r3) < (r4 + r1)) and ((r2+r1+r4) > r3) and (r3>=r2 and r3>=r1 and r3>=r4):
                                # assume r6 is half of r3 (coupler point in the middle)
                                # r6 = r3/2


                                r6_candidates = np.linspace(1,r3,resolution)
                                r1 = np.round(r1,3)
                                r2 = np.round(r2,3)
                                r3 = np.round(r3,3)
                                r4 = np.round(r4,3)
                                r7 = np.round(r7,3)

                                # for r6 in r6_candidates:
                                r6 = 1.0
                                r6 = np.round(r6,3)
                                coupler_trajectory_polar = []
                                coupler_trajectory_cartesian_x = []
                                coupler_trajectory_cartesian_y = []
                                target_trajectory_polar = []
                                target_trajectory_cartesian_x = []
                                target_trajectory_cartesian_y = []
                                trajectory_possible = True
                                for theta2 in theta2_range:
                                    x_coupler, y_coupler,x_coupler2, y_coupler2, trajectory_possible = calculate_coupler_position(
                                        r1, r2, r3, r4, theta2, r6, r7)
                                    # x_coupler, y_coupler,trajectory_possible = calculate_coupler_position(
                                        # r1, r2, r3, r4, theta2, r6, r7)
                                    if trajectory_possible:
                                        # ensure position within domain of curve
                                        # print(x_coupler)
                                        if abs(x_coupler)<=1 :
                                            x_target = x_coupler
                                            # equation converted from polar to cartesian
                                            y_target = sqrt((x_target**2 - x_target**4))
                                            coupler_trajectory_cartesian_x.append(x_coupler)
                                            coupler_trajectory_cartesian_y.append(y_coupler)
                                            target_trajectory_cartesian_x.append(x_target)
                                            target_trajectory_cartesian_y.append(y_target)

                                        if abs(x_coupler2)<=1:
                                            x_target2 = x_coupler2
                                            y_target_2 = - sqrt((x_target2**2 - x_target2**4))
                                            coupler_trajectory_cartesian_x.append(x_coupler2)
                                            coupler_trajectory_cartesian_y.append(y_coupler2)
                                            target_trajectory_cartesian_x.append(x_target2)
                                            target_trajectory_cartesian_y.append(y_target_2)

                                            # use this for polar
                                            # if theta_touse < np.pi/2
                                            #     theta_touse = theta_coupler
                                                # r_target = sqrt(r6**2 * cos(2 * theta_touse) * (1 / cos(theta_touse))**4)
                                                # x_target = r_target*cos(theta_touse)
                                                # y_target = r_target*sin(theta_touse)
                                                # coupler_trajectory_polar.append(r_coupler)
                                                # coupler_trajectory_polar.append(r_coupler)
                                                # target_trajectory_polar.append(r_target)
                                try:
                                    r6 = r3/2
                                    plt.scatter(target_trajectory_cartesian_x, target_trajectory_cartesian_y, c='r')
                                    plt.scatter(coupler_trajectory_cartesian_x, coupler_trajectory_cartesian_y, c= 'b')
                                    plt.xlabel("X-coordinate (inches)")
                                    plt.ylabel("y-coordinate (inches)")
                                    plt.legend(["Target trajectory", "Coupler trajectory"])
                                    # plt.title(f"r1:{r1}, r2:{r2}, r3:{r3}, r4:{r4}, r6:{r6}, r7:{r7}, rmse:{np.round(rmse,3)}, max_error:{np.round(maxerror,3)}")
                                    plt.show()
                                    if len(target_trajectory_cartesian_x) >=80:
                                        # maxy = max(target_trajectory_cartesian_y)
                                        rmse = sqrt(mean_squared_error((target_trajectory_cartesian_y), (coupler_trajectory_cartesian_y)))
                                        # print(len(target_trajectory_polar),len(coupler_trajectory))
                                        maxerror = max_error(target_trajectory_cartesian_y, coupler_trajectory_cartesian_y)
                                        if rmse < best_rmse:
                                            counter+=1
                                            # if counter% 5 == 0:
                                            plt.scatter(target_trajectory_cartesian_x, target_trajectory_cartesian_y, c='r')
                                            plt.scatter(coupler_trajectory_cartesian_x, coupler_trajectory_cartesian_y, c= 'b')
                                            plt.xlabel("X-coordinate (inches)")
                                            plt.ylabel("y-coordinate (inches)")
                                            plt.legend(["Target trajectory", "Coupler trajectory"])
                                            plt.title(f"r1:{r1}, r2:{r2}, r3:{r3}, r4:{r4}, r6:{r6}, r7:{r7}, rmse:{np.round(rmse,3)}, max_error:{np.round(maxerror,3)}")
                                            plt.show()

                                            best_rmse = rmse
                                            best_trajectory_cartesian = [coupler_trajectory_cartesian_x, coupler_trajectory_cartesian_y]
                                            best_trajectory_polar = coupler_trajectory_polar
                                            best_target_trajectory = [target_trajectory_cartesian_x, target_trajectory_cartesian_y]
                                            result = {"dimensions": {"r1":r1, "r2":r2, "r3":r3, "r4":r4, "r6":r6, "r7":r7, }, "rmse":rmse, "max_error":maxerror}
                                            experiment_results.append(result)
                                            print("Best Result: ", result)
                                except:
                                    continue

        print("Best RMSE: ", best_rmse)
        print("Best Result: ", result)
        print("Time taken: ", time.process_time() - start)

        plt.scatter(best_target_trajectory[0], best_target_trajectory[1], c='r')
        plt.scatter(best_trajectory_cartesian[0], best_trajectory_cartesian[1], c='b')
        plt.xlabel("X-coordinate (inches)")
        plt.ylabel("y-coordinate (inches)")
        plt.show()
