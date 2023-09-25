import math
import random


class Point:
    total_daisies = 0
    alive_daisies = 0

    age_of_death = 30
    maturity_age = 23

    req_resource = 20
    sexual_cost = 15
    clonal_cost = 12

    mutation_rate_low = 0.01
    mutation_rate_high = 0.05
    gene_length = 5

    colours = {
        "black": 0.25,
        "white": 0.75,
        "ground": 0.5
    }
    black = 0.25
    white = 0.75
    ground = 0.5

    x_dimension = 50
    y_dimension = 50

    flux = 1050  # Rate of energy received in Watts per metre**2.
    # Note: Value used is smaller than observed constant to simulate a younger star.

    def __init__(self, x_coord, y_coord):
        # Positional attributes of daisy or daisies
        self.x = x_coord
        self.y = y_coord
        self.coordinates = [(x_coord, y_coord)]  # Location of daisy on a 50x50 grid

        # Point attributes
        self.colour = Point.ground
        self.local_temp = None
        self.age = None  # Daisy age
        self.nutrients = None  # Number of accumulated nutrients
        self.genes = None
        self.opt_temp = None

    def grow_daisy(self, genes_list=None):
        Point.total_daisies += 1
        Point.alive_daisies += 1
        # Check if daisy exists, this is to reuse this method for new daisies
        self.genes = genes_list
        if genes_list is None:
            self.genes = [None] * Point.gene_length
            for i in range(Point.gene_length):
                self.genes[i] = random.randint(1, 10) / 10
        self.expressed_colour()
        self.expressed_opt_temp()

    def expressed_colour(self):
        prob_list = []
        total_val = self.genes[0] + self.genes[1]
        # Up to three colours
        for i in range(2):
            prob_list.append(self.genes[i] / total_val)
        res = self.pick_one(prob_list)
        self.colour = list(Point.colours.values())[res]

    def expressed_opt_temp(self):
        # The expressed genes would be closest to local temperature
        allele_A = self.local_temp * (self.colour + self.genes[3])
        allele_B = self.local_temp * (self.colour + self.genes[4])
        delta_A = abs(self.local_temp - allele_A)
        delta_B = abs(self.local_temp - allele_B)
        if delta_B > delta_A:
            self.opt_temp = allele_A
        else:
            self.opt_temp = allele_B

    @staticmethod
    def pick_one(probabilities):
        index = 0
        r = random.randint(0, 10) / 10
        while r > 0:
            r = r - probabilities[index]
            index += 1
            if index >= len(probabilities):
                break
        index -= 1
        return index

    def solar_factor(self):
        """Generates number between 0.8 - 1.2 based on y-coordinate between North and South pole and rounded to 2 d.p.

        :param int y_coord: Specified y-coordinate on map

        :rtype: float
        :return: Multiplier
        """
        return round(1.2 - 0.00064 * (self.y - 25) ** 2, 2)

    def check_pos(self):
        """Checks if daisy is present on a point

        :param int x: X coordinate of daisy
        :param int y: Y coordinate of daisy

        :rtype: int
        :return: 0 for not present, 1 for present
        """
        return self.colour == Point.ground

    def calc_temp(self, a_d, lumen):
        """Works out temperature at a point

        :param double flux: Solar flux constant
        :param double lumen: Solar luminosity
        :param double a_d: Planetary albedo

        :rtype: double
        :return: Average planetary temperature
        """
        sigma = 5.67037e-8  # Stefan-Boltzmann constant
        abs_zero = 273.15  # Used to calculate temperature in celsius
        q = 20

        temp_d = (((self.solar_factor() * Point.flux * lumen * (1 - a_d)) / sigma) ** 0.25) - abs_zero
        self.local_temp = q * (a_d - self.colour) + temp_d

    def find_neighbours(self):
        # Six points going out from the chosen point
        # Return all valid positions
        delta = [(0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1)]
        neighbours = []
        # Unpack change x and y
        for dx, dy in delta:
            x = self.x + dx
            y = self.y + dy
            if self.is_valid_point(x, y):
                neighbours.append((x, y))
        return neighbours

    def possible_points(self):
        delta = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2),
                 (2, 3), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (5, 0), (-1, 0), (-1, 1), (-1, 2), (-1, 3), (-1, 4),
                 (-2, 0), (-2, 1), (-2, 2), (-2, 3), (-3, 0), (-3, 1), (-3, 2), (-4, 0), (-4, 1), (-5, 0), (0, -1),
                 (1, -1), (2, -1), (3, -1), (4, -1), (0, -2), (1, -2), (2, -2), (3, -2), (0, -3), (1, -3), (2, -3),
                 (0, -4), (1, -4), (0, -5), (-1, -1), (-2, -1), (-3, -1), (-4, -1), (-1, -2), (-2, -2), (-3, -2),
                 (-1, -3), (-2, -3), (-1, -4)]
        points = []
        # Unpack change x and y
        for dx, dy in delta:
            x = self.x + dx
            y = self.y + dy
            if self.is_valid_point(x, y):
                points.append((x, y))
        return points

    @staticmethod
    def is_valid_point(x, y):
        # Returns boolean value
        return -1 < x < Point.x_dimension and -1 < y < Point.y_dimension

    def beta_y(self, temp_y, c=0.003265):
        """Daisy growth rate function

        :param double temp_y: Local temperature experienced by the daisies
        :param double opt_temp: Optimum temperature for black and white daisies to thrive at
        :param double c: Determines quadratic width, allowing growth to range from 5 to 40 degrees celsius on a negative parabolic curve

        :rtype: double
        :return: A value between 0 - 1, 0 indicating no growth, 1 indicating maximal growth
        """
        min_growth = round(self.opt_temp - math.sqrt(1 / c), 1)
        max_growth = round(self.opt_temp + math.sqrt(1 / c), 1)
        if min_growth <= temp_y <= max_growth:
            return 1 - c * (self.opt_temp - temp_y) ** 2
        else:
            return 0

    def grow(self, beta):
        if self.age >= Point.age_of_death:
            # Now an empty patch, bare ground
            self.age = None
            self.nutrients = None
            self.genes = None
            Point.alive_daisies -= 1
            return False
        else:
            self.nutrients += 5 * beta
            self.age += 1
            return True

    def s_reproduce(self, partner):
        # Produces progeny
        genes = [None] * Point.gene_length
        crossover = random.randint(0, Point.gene_length)
        for i in range(Point.gene_length):
            if i > crossover:
                genes[i] = self.genes[i]
            else:
                genes[i] = partner[i]
        return genes

    def mutate_low(self):
        # For those wondering, there are two variables for mutation rate, high and low. It is set to low because of the
        # large population size we have of daisies, variety does not need to be pushed for when handling with large
        # population as it's size makes up for it
        for i in range(Point.gene_length):
            rand = random.randint(0, 10) / 10
            if rand < Point.mutation_rate_low:
                self.genes[i] = random.randint(1, 10) / 10

    def mutate_high(self):
        for i in range(Point.gene_length):
            rand = random.randint(1, 10) / 10
            if rand < Point.mutation_rate_high:
                self.genes[i] = random.randint(1, 10) / 10

    def allocate_nutrients(self):
        self.nutrients = random.randint(2, 5)

    def randomise_age(self):
        self.age = random.randint(0, 15)

    def __str__(self):
        return "Coordinates: " + str(self.coordinates) + ", Colour: " + str(self.colour) + \
               ", Temperature: " + str(self.local_temp) + ", Age: " + str(self.age) + \
               ", Nutrients: " + str(self.nutrients) + ", Genes: " + str(self.genes) + \
               ", Optimum temperature: " + str(self.opt_temp)

    def __del__(self):
        return "Deleted daisy at coordinates: " + str(self.coordinates)