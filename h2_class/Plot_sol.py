# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 20:30:29 2023

@author: s1859154
"""
import matplotlib.pyplot as plt

def plot_sol(m, p):
    plt.figure()
    all_loc = p.loc
    W_col = ['blue' for i in range(p.n)]
    for k in p.K:
        plt.scatter(p.loc_supply[k][0], p.loc_supply[k][1], c='orange')
    for g in p.G_data:
        plt.scatter(p.loc_gas[g][0], p.loc_gas[g][1], c='blue')
    for i in p.I:
        if m.x[i, i].solution_value > 0.9:
            plt.scatter(p.loc[i][0], p.loc[i][1], c='green')
            plt.annotate(i, (p.loc[i][0]+2, p.loc[i][1]))
        else:
            plt.scatter(p.loc[i][0], p.loc[i][1], c='red')
            plt.annotate(i, (p.loc[i][0]+2, p.loc[i][1]))
    for i in p.I:
        for k in p.K:
            if m.T[i, k].solution_value > 0.9:
                plt.plot([p.loc[i][0], p.loc_supply[k][0]], [
                         p.loc[i][1], p.loc_supply[k][1]], c='orange')
    for i in p.I:
        for g in p.M:
            if m.u[g, i].solution_value > 0.9:
                if g >= len(p.I):
                    plt.plot([p.loc[i][0], p.loc_gas[g-p.n][0]],
                             [p.loc[i][1], p.loc_gas[g-p.n][1]], c='blue')
                else:
                    plt.plot([p.loc[i][0], p.loc[g][0]], [
                             p.loc[i][1], p.loc[g][1]], c='blue')
    plt.axis([p.width_2, p.width, p.width_2, p.width])
    plt.grid()
    fig = plt.gcf()

    fig.set_size_inches(8, 8)