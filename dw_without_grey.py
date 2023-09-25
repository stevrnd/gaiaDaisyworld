# -*- coding: utf-8 -*-
"""

Module  : CMP-6013Y - CMP Third Year Project
File    : daisyworld.py
Date    : Thursday 26 November 2020
Desc.   : Modelling the Gaia hypothesis through Daisyworld with a genetic algorithm introduced
History : 26/11/2020 - v1.0 - Created project file, created and appropriately named empty functions
          29/11/2020 - v1.1 - Added random number generator and radiation factor function
          13/12/2020 - v1.2 - Created dictionary for initial conditions of Daisyworld
          14/12/2020 - v1.3 - Daisy class added
          17/12/2020 - v1.4 - Updated solar_factor method
          19/12/2020 - v1.5 - Added basic method calc_opt_temp to set an optimum temperature for each daisy
          23/12/2020 - v1.6 -

"""
import math
import random

import numpy as np
import matplotlib.pyplot as plt
from point_without_grey import Point

__author__ = "Steven Diep"
__maintainer__ = "Steven Diep"
__email__ = "steven_diep@hotmail.co.uk"
__status__ = "Prototype"  # "Development" "Prototype" "Production"


class Daisyworld:
    def __init__(self, x_dim, y_dim, luminosities, init_pop):
        self.num_b = 0  # Number of black daisies
        self.num_w = 0  # Number of white daisies
        self.x_dim = x_dim
        self.y_dim = y_dim
        self.luminosities = luminosities
        self.init_pop = init_pop
        self.points = dict()
        self.generation = 0

        a_d = self.calc_avg_albedo()
        for x in range(self.x_dim):
            for y in range(self.y_dim):
                point = Point(x_coord=x,
                              y_coord=y)
                point.calc_temp(a_d, self.luminosities[0])
                self.points[(x, y)] = point

        i = 0
        while i < self.init_pop:
            x = random.randint(0, self.x_dim-1)
            y = random.randint(0.2*self.y_dim, 0.8*self.y_dim-1)
            daisy = self.points.get((x, y))
            daisy.allocate_nutrients()
            daisy.randomise_age()
            daisy.grow_daisy()
            if daisy.colour == Point.black:
                self.num_b += 1
            elif daisy.colour == Point.white:
                self.num_w += 1
            i += 1

    def calc_avg_albedo(self):
        """Calculates average albedo of Daisyworld

        :rtype: float
        :return: Average albedo with a range of 0 - 1
        """
        num_points = self.x_dim * self.y_dim
        area_b = self.num_b/num_points
        area_w = self.num_w/num_points
        u_area = 1 - (area_b + area_w)
        return Point.black*area_b + Point.white*area_w + Point.ground*u_area

    def find_diffuse_temp(self, neighbours, point):
        # Finds average temperature of point
        num_neighbours = len(neighbours) + 1
        total_temp = 0
        for neighbour in neighbours:
            total_temp += self.points.get(neighbour).local_temp
        return (total_temp + point)/num_neighbours

    def run(self):
        avg_temps = []
        avg_albedos = []
        num_black = []
        num_white = []
        for lumen in self.luminosities:
            avg_albedo_per_cycle = []
            avg_temps_per_cycle = []
            t = 0
            # Goes through 5 cycles before rise in luminosity
            while t < 5:
                a_d = self.calc_avg_albedo()
                temp_map = []
                mature_daisies = []
                for x in range(self.x_dim):
                    y_map = []
                    for y in range(self.y_dim):
                        point = self.points.get((x, y))
                        point.calc_temp(a_d, lumen)
                        neighbours = point.find_neighbours()
                        smoothed_temp = self.find_diffuse_temp(neighbours, point.local_temp)
                        y_map.append(smoothed_temp)
                        # If it is a daisy, grow
                        if not point.check_pos():
                            beta = point.beta_y(smoothed_temp)
                            alive = point.grow(beta)
                            if not alive:
                                if point.colour == Point.black:
                                    self.num_b -= 1
                                elif point.colour == Point.white:
                                    self.num_w -= 1
                                point.colour = Point.ground
                            if point.age is None:
                                pass
                            elif point.age > Point.maturity_age and point.nutrients > Point.req_resource:
                                mature_daisies.append((x, y))
                    temp_map.append(y_map)
                # Randomise list of mature daisies so daisies closer to 0x0 will
                # not get an advantage in selection process
                random.shuffle(mature_daisies)
                # Selection phase happens
                if not mature_daisies:
                    continue
                # Now we have a set of mature daisies
                self.generation += 1
                # Begin selection process
                index = 0
                while index < len(mature_daisies):
                    # Sexually reproducing should be more favourable than clonally reproducing,
                    # however, daisy must be a certain age and have enough nutrients
                    best_fitness = 0
                    best_mate = []
                    # Perform a search for best mate available within a range
                    # Certainly not the best method of selection but I am short on time
                    # Work out fitness in terms of locality and resources
                    for i in range(1, len(mature_daisies)):
                        x1, y1 = mature_daisies[index]
                        x2, y2 = mature_daisies[i]
                        # Works out Euclidean distance between daisies
                        locality = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
                        fitness = self.points.get(mature_daisies[i]).nutrients / locality
                        # Daisies that are too far are not counted
                        # Wanted locality to also be a gene trait but not enough time
                        if locality <= 7 and fitness > best_fitness:
                            best_fitness = fitness
                            best_mate.append(mature_daisies[index + 1])
                    # If no mates are found for daisy in list, clonally reproduce
                    if not best_mate:
                        clones = random.randint(3, 7)
                        # Daisy can fail to have offspring
                        if clones == 0:
                            self.points.get(mature_daisies[index]).nutrients -= Point.clonal_cost
                            mature_daisies.remove(mature_daisies[index])
                            continue
                        possible_points = self.points.get(mature_daisies[index]).possible_points()
                        for clone in range(clones):
                            while True:
                                # Randomly allocated position for child daisy
                                allocated_position = random.choice(possible_points)
                                child_daisy = self.points.get(allocated_position)
                                # Checks if daisy is already present on point, if present then remove coordinate
                                if not child_daisy.check_pos():
                                    # All possible areas are occupied by daisies
                                    # This daisy died from overcrowding
                                    break
                                elif child_daisy.check_pos():
                                    # Set coordinates of child
                                    child_genes = self.points.get(mature_daisies[index]).genes
                                    child_daisy.grow_daisy(child_genes)
                                    child_daisy.age = 0
                                    child_daisy.allocate_nutrients()
                                    # Asexual reproduction does not introduce enough variety
                                    # to planet, therefore, higher mutation rate for selfing is increased from 1% to 5%
                                    child_daisy.mutate_high()
                                    if child_daisy.colour == Point.black:
                                        self.num_b += 1
                                    elif child_daisy.colour == Point.white:
                                        self.num_w += 1
                                    break
                        self.points.get(mature_daisies[index]).nutrients -= Point.clonal_cost
                        mature_daisies.remove(mature_daisies[index])
                    else:
                        # Reproduce here, can produce between 0 - 3 children
                        daisy = self.points.get(mature_daisies[index])
                        # Final element of the best_mate list is always most fit daisy
                        daisy_mate = self.points.get(best_mate[-1]).genes
                        children = random.randint(3, 7)
                        if children == 0:
                            self.points.get(mature_daisies[index]).nutrients -= Point.sexual_cost
                            mature_daisies.remove(mature_daisies[index])
                            self.points.get(best_mate[-1]).nutrients -= Point.sexual_cost
                            mature_daisies.remove(best_mate[-1])
                            continue
                        # Get midpoint between parents
                        x1, y1 = mature_daisies[index]
                        x2, y2 = best_mate[-1]
                        x_mid = math.floor((x1 + x2)/2)
                        y_mid = math.floor((y1 + y1)/2)
                        # Find range around midpoint for dispersal
                        possible_points = self.points.get((x_mid, y_mid)).possible_points()
                        for child in range(children):
                            while True:
                                # Randomly allocated position for child daisy
                                allocated_position = random.choice(possible_points)
                                child_daisy = self.points.get(allocated_position)
                                if not child_daisy.check_pos():
                                    # All possible areas are occupied by daisies
                                    # This daisy died from overcrowding
                                    break
                                elif child_daisy.check_pos():
                                    # Set coordinates of child
                                    child_genes = daisy.s_reproduce(daisy_mate)
                                    child_daisy.grow_daisy(child_genes)
                                    child_daisy.age = 0
                                    child_daisy.allocate_nutrients()
                                    child_daisy.mutate_low()
                                    if child_daisy.colour == Point.black:
                                        self.num_b += 1
                                    elif child_daisy.colour == Point.white:
                                        self.num_w += 1
                                    break
                        self.points.get(mature_daisies[index]).nutrients -= Point.sexual_cost
                        mature_daisies.remove(mature_daisies[index])
                        self.points.get(best_mate[-1]).nutrients -= Point.sexual_cost
                        mature_daisies.remove(best_mate[-1])
                t += 1
                num_points = self.x_dim * self.y_dim
                total_temp = np.sum(temp_map)
                avg_planet_temp = total_temp/num_points
                avg_albedo_per_cycle.append(self.calc_avg_albedo())
                avg_temps_per_cycle.append(avg_planet_temp)
            avg_albedo = sum(avg_albedo_per_cycle) / len(avg_albedo_per_cycle)
            avg_albedos.append(avg_albedo)
            avg_temp = sum(avg_temps_per_cycle) / len(avg_temps_per_cycle)
            avg_temps.append(avg_temp)
            num_black.append(self.num_b)
            num_white.append(self.num_w)
        plt.plot(self.luminosities, avg_temps, 'b')
        plt.title('Temperature over luminosity')
        plt.xlabel('Solar Luminosity')
        plt.ylabel('Temperature (°C)')
        plt.show()

        plt.plot(self.luminosities, avg_albedos, 'b')
        plt.title('Average albedo over luminosity')
        plt.xlabel('Solar Luminosity')
        plt.ylabel('Albedo')
        plt.show()

        plt.plot(self.luminosities, num_black, 'b', label='Black daisies')
        plt.plot(self.luminosities, num_white, 'g', label='White daisies')
        plt.legend(loc='upper right')
        plt.title('Number of daisies over luminosity')
        plt.xlabel('Solar Luminosity')
        plt.ylabel('Count')
        plt.show()