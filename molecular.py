"""
molecular.py is a unit of Kasuga computational package responsible for all operations with molecular geometries such as:
    1. Reading, writing and processing .cif files containing crystal packing data.
    2. Reading, writing and processing .xyz, .mol2 and .pdb files containing single molecule and molecular cluster data
    3. Reading and processing .out (Gaussian), .cube and .fchk files containing computational properties for single
        molecules and molecular clusters.
    4. Symmetry operations with single molecules and molecular clusters, molecular cluster generation.
    5. Modification of a molecular geometry: changing bond lengths, angles and dihedral angles, changing atomic and
        molecular properties.

As a rule of thumb, if you need to process and modify molecular geometries or read and store any property (such as
electron density) of a single molecule - it is most likely handled by molecular.py unit. Any other, more advanced
processing down the line (such as integration, IC and ISC calculations) is handled by other modules.

⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀⣀⣀⣀⣤⣤⣤⣤⣤⣤⣤⣤⣤⣀⣀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣤⣴⣶⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣶⣤⣀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⣴⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⣶⣿⢿⡿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⢿⣿⢿⡿⣿⡿⣿⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣾⣿⣯⣟⣯⣿⣿⣽⣾⣿⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⣡⢺⡆⠙⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣿⣷⣻⣯⣿⣷⣆⠀⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠟⢹⣿⣿⡿⢋⠞⠃⠀⣷⠀⠈⠻⣿⣿⣿⣿⠛⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣇⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠀⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣟⣡⣿⣿⡿⠧⠤⣾⡿⢋⠤⠊⠀⠀⠀⢸⡄⠀⠀⠈⠻⣿⣿⡦⠬⠽⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡆⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⠟⠁⠀⣸⣿⣯⣤⣤⣼⣏⡠⠎⠀⠀⠀⠀⠀⠈⢷⠀⠀⠀⠀⣈⣻⣷⣤⣤⣬⣻⣧⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⣼⣿⣿⣿⣿⣿⣿⡿⠋⢀⣴⢾⢏⣿⣿⣿⣿⣿⣮⣅⠀⠀⠀⠀⠀⠀⠀⠘⠙⠀⠀⠀⣩⣵⣿⣿⣿⣷⣮⡙⢳⣦⠻⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⣿⡿⠃⢠⣿⠁⢠⣿⣿⣿⣿⣿⣿⣿⣿⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⣿⣿⣿⣿⣿⣿⣿⣷⠀⢻⣧⠈⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⣿⣿⣿⣿⣿⣿⡇⠀⢹⣿⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⠀⢸⣿⠀⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢀⣿⣿⣿⣿⣿⣿⣷⠀⠀⠁⠀⠈⢷⣤⣼⣿⣿⣿⠿⠗⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠧⢤⣾⣿⣿⣿⣿⠟⠀⠘⠉⠀⣿⣿⣿⣿⣿⣿⣿⡟⠀⠀⠀
⠀⠀⠀⠀⠀⠀⢰⣿⣿⣿⣿⣿⣿⣿⠆⠀⠀⠀⠀⠀⠈⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢠⣶⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⣿⣿⣿⣿⣿⣿⣿⣿⣿⡆
⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣄⣀⣀⠀⣁⢐⣀⣠⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣸⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀
⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣾⠿⠟⠉⠋⠉⠋⠙⠛⢿⡇⠀⠀⠀⠀⠀⠀⠀⣀⢀⣠⣴⣿⣿⣿⣿⣿⣿⣿⣿⣿⠁⠀⠀
⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣤⣀⣀⠀⠀⠀⠀⠀⠙⣦⣀⡆⠈⡄⠀⣘⣠⡞⠁⠀⠀⢀⡀⣄⣶⣶⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀
⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣶⣶⣀⣰⣀⢈⠙⠛⠛⠛⠛⠋⠉⣀⣶⣰⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⣸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⢛⣧⠉⠉⠙⠒⠛⠛⠛⠚⠛⠉⠉⢰⡿⠿⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣿⡿⠟⠋⠁⠰⣯⠟⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠄⢺⣡⡆⠀⠈⠻⢿⡛⢹⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀
⠀⠀⠀⠀⠀⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠟⠊⠉⠀⠀⠀⠀⠀⢹⡀⠀⠈⠑⠲⢤⣄⣂⣈⣐⣀⣀⣠⣟⣀⠤⠔⠒⠀⠀⠉⠛⢿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠘⣿⣿⣿⣿⣿⣿⣿⣯⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⢷⡒⣒⠲⣒⢒⢲⣶⡒⠖⡲⢶⡞⠁⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⠘⣿⣿⣿⣿⣿⣿⣿⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣷⢠⣣⠔⡊⠖⢫⠘⡴⣱⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣿⣿⣿⣿⣿⣿⣿⡟⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⠀⠀⣨⣿⣇⢻⣿⣿⣿⣿⡢⣀⠀⠀⠀⠀⠀⣀⣀⣠⡀⠀⠈⣧⠳⡘⠴⡉⢆⠽⣿⠋⠀⠀⠀⠀⣀⡀⠀⠀⠀⠀⣼⣿⣿⣿⣿⢛⣿⠿⡄⠀⠀⠀⠀⠀
⠀⠀⠀⠀⠀⢄⡾⣣⢁⡛⢣⢉⡻⡿⢿⣷⢦⣉⠓⢒⣊⣭⡭⢀⡏⠀⠀⠀⠈⢳⣱⡢⠑⡌⡿⠁⠀⠀⠀⠀⢀⡇⠨⣭⣖⠢⣼⣿⣿⠿⠋⣄⣾⠏⡒⣽⣦⠀⠀⠀⠀
⠀⠀⠀⠀⣠⡟⢚⢄⡒⢌⢂⠣⡑⠤⡃⢍⠣⡌⢍⣃⢲⣏⡀⠸⠤⣀⣀⠀⠀⠀⠹⣦⠣⣼⠃⠀⢀⣀⡤⠤⠚⢁⣀⣷⠄⢿⡟⢫⢁⠎⡱⢌⠡⢃⠜⣘⡹⣇⠀⠀⠀
⠀⠀⠀⣰⠏⣷⡃⢆⡘⠤⢃⣌⡇⠖⡩⢌⠢⣑⠊⡔⢢⢒⠩⣉⠓⡒⠦⠭⠭⢥⣒⣚⣻⣗⣈⠭⠥⢖⢒⡚⠩⡍⡒⢤⡉⠆⡍⢒⠌⡒⣡⢊⡴⢃⠜⣠⢑⡼⠀⠀⠀
⠀⠀⣰⢫⢰⣉⣿⡌⣴⣁⢣⣼⡏⡜⢰⣈⢱⢠⢩⣐⢡⢊⢱⢠⡍⣰⢉⡜⣤⢓⢸⣏⣀⣀⣿⡘⢢⢡⡒⣌⡑⣆⡑⢢⣌⢱⣈⣬⡘⣤⡁⣎⡟⡌⢢⢡⣾⢃⠀⠀⠀
⠀⢰⡟⡏⠷⡙⢎⣿⡶⡉⢇⣿⡷⢙⠳⡉⠇⡎⢳⠉⢇⢋⠞⢳⢁⠳⠞⡁⢇⠞⡸⢹⡟⣿⡏⢇⡎⠷⡙⢎⡹⢈⡹⢃⠞⡰⣉⠆⢳⠆⢳⡘⣿⡁⡏⣾⠇⢏⠀⠀⠀
⢀⡾⣥⢋⠴⡁⠆⣿⣷⣱⣿⣻⡇⢎⠰⡁⢎⠰⡁⢎⠢⢌⡘⠄⡎⠱⡨⢑⠌⡢⠅⣿⡇⣻⣏⠦⡘⠰⡉⠆⡔⢣⢐⡡⢊⠔⡰⢨⡁⢎⡡⢘⣿⣷⣼⡟⣌⠢⠀⠀⠀
⢸⢇⠹⢌⠢⢅⠓⢾⣟⡿⣷⣿⠇⡊⢔⠡⡊⢔⠡⡊⠔⢢⠘⡰⢈⡱⢐⡉⢆⠡⢍⣿⢦⢙⣿⢂⡍⠱⡘⠌⡔⢡⢂⠒⡡⢊⠔⡡⡘⢄⡒⡡⣿⣿⣿⠒⡄⢣⠀⠀⠀
⣿⢂⠣⢌⠣⡘⢌⠺⣿⣽⢻⣿⢀⠓⡌⢒⠡⢊⠒⡡⡉⢆⡑⢢⢁⠦⢡⠘⠤⣉⢺⣟⢸⡌⣿⡆⡌⢱⠈⡜⢠⠃⡌⢒⠡⢊⠔⡡⡘⢄⠆⡱⣿⢯⡽⢃⠜⡠⠀⠀⠀
⢠⡇⢎⢂⠧⢌⠱⡈⠜⣿⣞⣯⡗⢨⡘⠤⢃⠥⢃⡱⠤⡑⢢⠘⡄⡊⠔⡡⢊⠔⡤⢻⡏⡤⡇⢿⡧⡘⢄⠣⡘⢄⠣⡘⠤⡉⢆⠚⠤⡑⢊⡔⢡⣿⣻⡏⠴⡘⠰⠀⠀⠀
⣼⠓⡌⡒⢌⣂⠣⡉⢆⡹⢯⣷⡏⠔⡨⢒⡉⠆⠥⢂⡱⢈⠆⡱⠠⢅⠓⡄⢣⠘⠤⣿⡇⢇⡹⢸⣷⢁⠎⡰⢁⠎⡰⢁⡒⢡⠊⡜⢠⠃⡥⠘⠤⣿⣳⡟⠤⡑⢃⠀⠀⠀
"""

import os
import kasuga_io
import numpy as np
from constants import element_weight
from constants import covalent_radius
from constants import HM2Hall
from constants import SymOpsHall


def log_notation_to_float(s: str):
    if s.find("E") != -1:
        cut = s.split("E")
    elif s.find("D") != -1:
        cut = s.split("D")
    else:
        return None
    return float(cut[0]) * (10 ** float(cut[1]))


class Vector:
    """
    Vector is a base class containing coordinates represented as a np.array(3)
    """

    def __init__(self):
        self.coord = np.zeros(3)

    def transform(self, matrix: np.array):
        """
        General method for vector coordinate transformation using a transformation matrix.
        :param matrix: transformation 3x3 matrix (np.array)
        """
        self.coord = np.matmul(self.coord, matrix)

    def invert(self, inv_coord=np.zeros(3)):
        """
        Invert vector around arbitrary center.
        :param inv_coord: inversion center coordinates (np.array, [0,0,0] by default)
        """
        shift = inv_coord - self.coord
        self.coord += 2 * shift

    def mirror(self, normal=np.array([1, 0, 0]), point=np.zeros(3)):
        """
        Mirror vector in an arbitrary mirror plane
        :param normal: normal vector of a mirror plane (np.array, [1,0,0] by default)
        :param point: arbitrary point that belongs to a mirror plane (np.array, [0,0,0] by default)
        :return:
        """
        # normalize n just to be safe
        n = normal / np.linalg.norm(normal)
        # plane equation coefficients
        a = normal[0]
        b = normal[1]
        c = normal[2]
        d = -1 * (a * point[0] + b * point[1] + c * point[2])
        # distance between a point (our vector) and a mirror plane
        distance = abs(a * self.coord[0] + b * self.coord[1] + c * self.coord[3] + d) / np.sqrt(
            a ** 2 + b ** 2 + c ** 2)
        # provided normal vector can either face the same direction as mirrored point or the opposite
        test1 = self.coord + 2 * n * distance
        test2 = self.coord - 2 * n * distance
        distance_test1 = abs(a * test1[0] + b * test1[1] + c * test1[3] + d) / np.sqrt(a ** 2 + b ** 2 + c ** 2)
        distance_test2 = abs(a * test2[0] + b * test2[1] + c * test2[3] + d) / np.sqrt(a ** 2 + b ** 2 + c ** 2)
        if distance_test1 < distance_test2:
            # same direction
            self.coord = test1
        else:
            # opposite direction
            self.coord = test2

    def xyz_mirror(self, plane="xy", plane_point=np.zeros(3)):
        """
        Simplified mirror method to reflect vector in xy, xz and yz planes
        :param plane: string representing one of default planes: "xy", "xz", "yz" or "ab", "ac", "bc"
        :param plane_point: arbitrary point that belongs to a mirror plane (np.array, [0,0,0] by default)
        :return:
        """
        match plane:
            case "xy" | "ab":
                self.mirror(np.array([0, 0, 1]), plane_point)
            case "xz" | "ac":
                self.mirror(np.array([0, 1, 0]), plane_point)
            case "yz" | "bc":
                self.mirror(np.array([1, 0, 0]), plane_point)

    def rotate(self, angle: float, axis_vector: np.array, axis_point=np.zeros(3)):
        """
        Rotate Vector around arbitrary axis.
        :param angle: angle of rotation (in degrees)
        :param axis_point: point of origin for rotation axis (np.array) (np.array, [0,0,0] by default)
        :param axis_vector:
        """
        axis_diff = self.coord - axis_point
        axis_vector = axis_vector / np.linalg.norm(axis_vector)
        angle_rad = np.deg2rad(angle)
        matrix = np.zeros((3, 3))
        matrix[0, 0] = np.cos(angle_rad) + axis_vector[0] ** 2 * (1 - np.cos(angle_rad))
        matrix[0, 1] = axis_vector[0] * axis_vector[1] * (1 - np.cos(angle_rad)) - axis_vector[2] * np.sin(angle_rad)
        matrix[0, 2] = axis_vector[0] * axis_vector[0] * (1 - np.cos(angle_rad)) + axis_vector[1] * np.sin(angle_rad)
        matrix[1, 1] = np.cos(angle_rad) + axis_vector[1] ** 2 * (1 - np.cos(angle_rad))
        matrix[1, 2] = axis_vector[1] * axis_vector[2] * (1 - np.cos(angle_rad)) - axis_vector[0] * np.sin(angle_rad)
        matrix[2, 2] = np.cos(angle_rad) + axis_vector[2] ** 2 * (1 - np.cos(angle_rad))
        matrix[1, 0] = matrix[0, 1]
        matrix[2, 0] = matrix[0, 2]
        axis_diff_rotated = np.matmul(axis_diff, matrix)
        axis_translation = axis_diff - axis_diff_rotated
        self.coord += axis_translation

    def improper_rotate(self, angle: float, axis_vector=np.zeros(3), point=np.zeros(3)):
        self.rotate(angle, axis_vector, point)
        self.mirror(axis_vector, point)

    def screw_axis(self, angle: float, axis_vector=np.zeros(3), point=np.zeros(3), translation_vector=np.zeros(3)):
        self.rotate(angle, axis_vector, point)
        self.coord += translation_vector

    def glide_plane(self, normal=np.array([1, 0, 0]), point=np.zeros(3), translation_vector=np.zeros(3)):
        self.mirror(normal, point)
        self.coord += translation_vector

    def distance(self, v2):
        delta = self.coord - v2.coord
        return np.linalg.norm(delta)

    def distance_rough(self, v2):
        delta = self.coord - v2.coord
        return max(float(abs(delta[0])), float(abs(delta[1])), float(abs(delta[2])))


class ConnectivityGraph:

    def __init__(self, size):
        self.size = size
        self.nodes = np.zeros((size, size))

    def flood_fill_search(self, startpoint: int, excluded=()):
        inside = np.zeros(self.size)  # points that are connected to the starting one
        checked = np.zeros(self.size)  # points that we have already checked
        while True:
            for i in range(self.size):
                if self.nodes[startpoint, i] == 1 or self.nodes[i, startpoint] == 1:
                    if i not in excluded:
                        inside[i] = 1  # include all points that are connected to current start point
            checked[startpoint] = 1  # so, we've checked current start point
            count_inside = 0
            count_checked = 0
            for i in range(self.size):
                if checked[i] == 1:
                    count_checked += 1
                if inside[i] == 1:
                    count_inside += 1
                if checked[i] != 0 and inside[i] == 1:
                    startpoint = i
            if count_inside == count_checked:
                return inside

    def subsets_connected(self, subset_one=np.zeros(1), subset_two=np.zeros(1)):
        for i1 in range(subset_one.size):
            for i2 in range(subset_two.size):
                if self.nodes[i1, i2] == 1 or self.nodes[i2, i1] == 1:
                    return True
        return False


class Atom(Vector):

    def assign_weight(self):
        if self.symbol in element_weight:
            self.weight = element_weight[self.symbol]
        else:
            kasuga_io.quit_with_error(f'Unrecognized {self.symbol} atom encountered!')

    def __init__(self, symbol=""):
        self.weight = 0.0  # Atomic weight
        self.charge = 0.0
        self.symbol = ""  # Chemical symbol of an atom
        super().__init__()
        if symbol != "":
            self.assign_weight()

    def __eq__(self, other):
        distance = np.linalg.norm(self.coord - other.coord)
        if distance < 0.001 and self.symbol == other.symbol:
            return True
        else:
            return False

    def __ne__(self, other):
        distance = np.linalg.norm(self.coord - other.coord)
        if distance < 0.001 and self.symbol == other.symbol:
            return False
        else:
            return True

    def connected(self, b, simplified=False, cutoff=0.025):
        if simplified:
            d = Vector.distance_rough(self, b)
        else:
            d = Vector.distance(self, b)
        if d <= cutoff:
            return True
        else:
            return False


class Molecule:

    def __init__(self):
        self.atoms = []
        self.mass_center = None
        self.molecular_formula = None
        self.connectivity_graph = None
        self.inertia_eigenvectors = None
        self.inertia_eigenvalues = None
        self.inertia_vector_x, self.inertia_vector_y, self.inertia_vector_z = None, None, None
        self.symmetrized = False
        self.dipole_moment = Vector()
        self.quadrupole_moment = np.zeros((3, 3))
        self.point_group = ""

    def add_atom(self, atom: str, v: Vector):
        if atom in element_weight:
            new_atom = Atom()
            new_atom.coord = v
            new_atom.symbol = atom
            new_atom.assign_weight()
            self.atoms.append(new_atom)
        else:
            kasuga_io.quit_with_error(f'add_atom method error: unknown atom {atom}')

    def rebuild_connectivity(self):
        self.connectivity_graph = ConnectivityGraph(len(self.atoms))
        for i1 in range(len(self.atoms) - 1):
            for i2 in range(len(self.atoms) - 1):
                diff = self.atoms[i1].coord - self.atoms[i2].coord
                vdw_sum = covalent_radius[self.atoms[i1].symbol] + covalent_radius[self.atoms[i2].symbol]
                if abs(np.linalg.norm(diff) - vdw_sum) < 0.05:
                    self.connectivity_graph.nodes[i1, i2] = 1
                    self.connectivity_graph.nodes[i2, i1] = 1

    def separate_molecules(self):
        checked = []
        vectors = []
        flag = True
        while flag:
            for i in range(self.connectivity_graph.size - 1):
                checked.append(False)
            for i in range(self.connectivity_graph.size - 1):
                if not checked:
                    vector = self.connectivity_graph.flood_fill_search(startpoint=i)
                    break
            for i in range(self.connectivity_graph.size - 1):
                if vector[i] == 1:
                    checked[i] = True
                else:
                    checked[i] = False
            vectors.append(vector)
            check_count = 0
            for i in range(self.connectivity_graph.size - 1):
                if checked:
                    check_count += 1
            if check_count == self.connectivity_graph.size:
                flag = False
        molecules = []
        for v in vectors:
            new_mol = Molecule()
            for i in range(self.connectivity_graph.size - 1):
                if v[i] == 1:
                    new_mol.atoms.append(self.atoms[i])
            molecules.append(new_mol)
        if len(molecules) == 1:
            return None
        return molecules

    def __add__(self, other):
        for i in other.atoms:
            self.atoms.append(i)
        self.rebuild_connectivity()

    def __sub__(self, other):
        for_deletion = []
        for i1 in other.atoms:
            for num, i2 in enumerate(self.atoms):
                if i1 == i2:
                    for_deletion.append(self.atoms[num])
        for i1 in for_deletion:
            if i1 in self.atoms:
                self.atoms.remove(i1)
        self.rebuild_connectivity()

    def __eq__(self, other):
        diff = 0.0
        # Simple tests first to potentially save the hassle
        if self.get_molecular_formula() != self.get_molecular_formula():
            return False
        # For symmetry cloned molecules it's safe to assume that the order of atoms is still the same
        # But generally it's not always the case, especially if molecules originate from different sources
        for i1 in range(len(self.atoms)):
            current_min = 1000
            for i2 in range(len(other.atoms)):
                # We look for the closest atom with the same symbol
                if self.atoms[i1].symbol == other.atoms[i2].symbol:
                    current_diff = self.atoms[i1].distance(other.atoms[i2])
                    if current_diff < current_min:
                        current_min = current_diff
            diff += current_min
        if diff < 0.05:
            return True
        else:
            return False

    def __ne__(self, other):
        diff = 0.0
        if self.get_molecular_formula() != self.get_molecular_formula():
            return False
        for i1 in range(len(self.atoms)):
            current_min = 1000
            for i2 in range(len(other.atoms)):
                if self.atoms[i1].symbol == other.atoms[i2].symbol:
                    current_diff = self.atoms[i1].distance(other.atoms[i2])
                    if current_diff < current_min:
                        current_min = current_diff
            diff += current_min
        if diff < 0.05:
            if diff < 0.05:
                return False
            else:
                return True

    def is_connected(self, other):
        for i1 in self.atoms:
            for i2 in other.atoms:
                if i1.connected(i2):
                    return True
        return False

    def get_mass_center(self):
        if self.mass_center is None:
            self.mass_center = Vector()
        mass = 0.0
        mass_vector = np.zeros((3, 1))
        for atom in self.atoms:
            mass += element_weight[atom.symbol]
            mass_vector += atom.coord * element_weight[atom.symbol]
        self.mass_center = mass_vector / mass
        return self.mass_center

    def get_molecular_formula(self):
        if self.molecular_formula is None:
            self.molecular_formula = ""
        atom_list = []
        count_list = []
        for i, atom in enumerate(self.atoms):
            if atom.symbol not in atom_list:
                atom_list.append(atom.symbol)
        atom_list.sort()
        for i1, atom_in_formula in enumerate(atom_list):
            count_list.append(0)
            for i2, atom_in_list in enumerate(self.atoms):
                if atom_in_formula == atom_in_list:
                    count_list[i1] += 1
        result = ""
        for i, atom in enumerate(atom_list):
            result += (atom + str(count_list[i]))
        self.molecular_formula = result
        return self.molecular_formula

    def get_connectivity_matrix(self):
        if self.connectivity_graph is None:
            self.connectivity_graph = ConnectivityGraph(len(self.atoms))
        for i1 in range(len(self.atoms)):
            for i2 in range(i1, len(self.atoms)):
                if self.atoms[i1].connected(self.atoms[i2]):
                    self.connectivity_graph.nodes[i1, i2] = 1
                    self.connectivity_graph.nodes[i2, i1] = 1
        return self.connectivity_graph.nodes

    def get_inertia_vectors(self):
        if self.inertia_eigenvalues is None or self.inertia_eigenvectors is None:
            self.inertia_eigenvectors = np.zeros((3, 3))
            self.inertia_eigenvalues = np.zeros(3)
        # First, we translate origin to mass center
        if self.mass_center is None:
            mass_center = self.get_mass_center()
        else:
            mass_center = self.mass_center
        atoms_mc_system = []
        for a in self.atoms:
            new_a = Atom()
            new_a.coord = a.coord - mass_center
            new_a.symbol = a.symbol
            atoms_mc_system.append(new_a)
        # Calculate inertia tensor
        for a in atoms_mc_system:
            # xx
            self.inertia_eigenvectors[0, 0] += element_weight[a.symbol] * (a.coord[1, 0] ** 2 + a.coord[2, 0] ** 2)
            # xy
            self.inertia_eigenvectors[0, 1] += -1 * element_weight[a.symbol] * a.coord[0, 0] * a.coord[1, 0]
            # xz
            self.inertia_eigenvectors[0, 2] += -1 * element_weight[a.symbol] * a.coord[0, 0] * a.coord[2, 0]
            # yy
            self.inertia_eigenvectors[1, 1] += element_weight[a.symbol] * (a.coord[0, 0] ** 2 + a.coord[2, 0] ** 2)
            # yz
            self.inertia_eigenvectors[1, 2] += -1 * element_weight[a.symbol] * a.coord[1, 0] * a.coord[2, 0]
            # zz
            self.inertia_eigenvectors[2, 2] += element_weight[a.symbol] * (a.coord[0, 0] ** 2 + a.coord[1, 0] ** 2)
        self.inertia_eigenvectors[1, 0] = self.inertia_eigenvectors[0, 1]
        self.inertia_eigenvectors[2, 0] = self.inertia_eigenvectors[0, 2]
        self.inertia_eigenvectors[2, 1] = self.inertia_eigenvectors[1, 2]
        # Calculate eigenvalues and eigenvectors of the inertia tensor
        self.inertia_eigenvalues, self.inertia_eigenvectors = np.linalg.eig(self.inertia_eigenvectors)
        # Assign eigenvectors to Cartesian axis: highest for Z, lowest for X
        internal_e = self.inertia_eigenvalues
        if self.inertia_vector_x is None or self.inertia_vector_y is None or self.inertia_vector_z is None:
            self.inertia_vector_x = Vector()
            self.inertia_vector_y = Vector()
            self.inertia_vector_z = Vector()
        for i in range(3):
            index = np.where(internal_e == internal_e.max())
            if self.inertia_vector_z.coord != np.zeros(3):
                self.inertia_vector_z = (self.inertia_eigenvalues[index[0] - 1:index[0], :] /
                                         np.linalg.norm(self.inertia_eigenvalues[index[0] - 1:index[0], :]))
            elif self.inertia_vector_y.coord != np.zeros(3):
                self.inertia_vector_y.coord = (self.inertia_eigenvalues[index[0] - 1:index[0], :] /
                                               np.linalg.norm(self.inertia_eigenvalues[index[0] - 1:index[0], :]))
            elif self.inertia_vector_x.coord != np.zeros(3):
                self.inertia_vector_x.coord = (self.inertia_eigenvalues[index[0] - 1:index[0], :] /
                                               np.linalg.norm(self.inertia_eigenvalues[index[0] - 1:index[0], :]))
            internal_e[index[0], 0] = -1.0
        return self.inertia_vector_x, self.inertia_vector_y, self.inertia_vector_z

    def match_rotation_to(self, other):
        # Extract principal axes for each molecule
        # Static stays in place, rotated is transformed
        rotated_x, rotated_y, rotated_z = self.get_inertia_tensor()
        static_x, static_y, static_z = other.get_inertia_tensor()
        # Since the vectors are stored in ((3,1)) shape we need transpose them first
        rotated_x = np.transpose(rotated_x)
        rotated_y = np.transpose(rotated_y)
        rotated_z = np.transpose(rotated_z)
        static_x = np.transpose(static_x)
        static_y = np.transpose(static_y)
        static_z = np.transpose(static_z)
        # Create matrices that rotate standard coordinate system 0 to molecule's principal axes
        rot_mat_rotated = np.array([[rotated_x], [rotated_y], [rotated_z]])
        rot_mat_static = np.array([[static_x], [static_y], [static_z]])
        # Combine two rotations: Rotated -> 0  (inverted 0 -> Rotated) and 0 -> Static
        final_rotation = np.matmul(np.linalg.inv(rot_mat_rotated), rot_mat_static)
        # We translate Rotated to 0 system, perform rotation, and translate it back
        mass_center = self.get_mass_center()
        for i in range(len(self.atoms)):
            self.atoms[i].coord -= mass_center
            self.atoms[i].coord = np.transpose(np.matmul(np.transpose(self.atoms[i].coord), final_rotation))
            self.atoms[i].coord += mass_center

    def change_bond(self, bond: tuple, delta: float):
        first_fragment = self.connectivity_graph.flood_fill_search(bond[0], (bond[1]))
        second_fragment = self.connectivity_graph.flood_fill_search(bond[1], (bond[0]))
        translation_vector = self.atoms[bond[0]].coord - self.atoms[bond[1]].coord
        translation_vector = np.linalg.norm(translation_vector)
        if self.connectivity_graph.subsets_connected(first_fragment, second_fragment):
            self.atoms[bond[0]].coord += delta * translation_vector / 2
            self.atoms[bond[1]].coord -= delta * translation_vector / 2
        else:
            for i in first_fragment:
                self.atoms[i].coord += delta * translation_vector / 2
            for i in second_fragment:
                self.atoms[i].coord -= delta * translation_vector / 2

    def change_angle(self, angle: tuple, delta: float):
        first_fragment = self.connectivity_graph.flood_fill_search(angle[0], (angle[1]))
        second_fragment = self.connectivity_graph.flood_fill_search(angle[2], (angle[1]))
        v1 = self.atoms[angle[0]].coord - self.atoms[angle[1]].coord
        v2 = self.atoms[angle[2]].coord - self.atoms[angle[1]].coord
        rotation_vector = np.cross(v1, v2)
        if self.connectivity_graph.subsets_connected(first_fragment, second_fragment):
            self.atoms[angle[0]].rotate(delta / 2, rotation_vector, self.atoms[angle[1]])
            self.atoms[angle[0]].rotate(-1 * delta / 2, rotation_vector, self.atoms[angle[1]])
        else:
            for i in first_fragment:
                self.atoms[i].rotate(delta / 2, rotation_vector, self.atoms[angle[1]])
            for i in second_fragment:
                self.atoms[i].rotate(-1 * delta / 2, rotation_vector, self.atoms[angle[1]])

    def change_dihedral(self, dihedral: tuple, delta: float):
        f_fragment = self.connectivity_graph.flood_fill_search(dihedral[0], (dihedral[1], dihedral[2], dihedral[3]))
        s_fragment = self.connectivity_graph.flood_fill_search(dihedral[3], (dihedral[0], dihedral[1], dihedral[2]))
        rotation_vector = self.atoms[dihedral[1]].coord - self.atoms[dihedral[2]].coord
        if self.connectivity_graph.subsets_connected(f_fragment, s_fragment):
            self.atoms[dihedral[0]].rotate(delta / 2, rotation_vector, self.atoms[dihedral[1]])
            self.atoms[dihedral[3]].rotate(-1 * delta / 2, rotation_vector, self.atoms[dihedral[1]])
        else:
            for i in f_fragment:
                self.atoms[i].rotate(delta / 2, rotation_vector, self.atoms[dihedral[1]])
            for i in s_fragment:
                self.atoms[i].rotate(-1 * delta / 2, rotation_vector, self.atoms[dihedral[1]])


class GaussianFile:

    class SCF:

        def __init__(self):
            self.RMS_DP_criteria = 0.0
            self.Max_DP_criteria = 0.0
            self.DeltaE_criteria = 0.0
            self.energy = []
            self.delta_energy = []
            self.final_energy = []
            self.RMSDP = []
            self.MaxDP = []
            self.OVMax = []
            self.occ_eigenvalues = []
            self.virt_eigenvalues = []

    def __init__(self):
        self.path = ""
        self.file_raw_contents = []
        self.__start_end = []
        self.version = ""
        self.revision_date = ""
        self.execution_date = ""
        self.link0_instructions = {}
        self.calculation_instructions = []
        self.iops_instructions = []
        self.calculation_title = ""
        self.geometries = []
        self.scf_stages = []

    def link1(self):

        def get_title(lines):
            for n, line in enumerate(lines):
                if line[1:4] == "***":
                    first = lines[n + 1].strip().split()
                    second = lines[n + 1].strip().split()
                    self.version = f'{first[0]} {first[1]}'
                    self.revision_date = first[2]
                    self.execution_date = second[0]
                    return n + 2

        def get_link0(start_pos: int, lines):
            for i1 in range(start_pos, len(lines)):
                a = lines[i1].strip()
                if a[0] == "-":
                    return i1 + 1
                else:
                    splitted = a.split("=")
                    self.link0_instructions[splitted[0].strip()[1:]] = splitted[1].strip()

        def get_calculation_instructions(start_pos: int, lines):
            for i1 in range(start_pos, len(lines)):
                a = lines[i1].strip()
                if a[0] == "-":
                    return i1 + 1
                else:
                    self.calculation_instructions.append(a)

        extracted_lines = self.file_raw_contents[self.__start_end[0]: self.__start_end[1]]
        last_line = get_title(extracted_lines)
        last_line = get_link0(last_line, extracted_lines)
        last_line = get_calculation_instructions(last_line, extracted_lines)
        for i in range(last_line, len(extracted_lines)):
            self.iops_instructions.append(extracted_lines[i].strip())

    def link101(self):
        extracted_lines = self.file_raw_contents[self.__start_end[0]: self.__start_end[1]]
        self.calculation_title = extracted_lines[1].strip()
        for num, s in enumerate(extracted_lines):
            if s[0] == " ":
                return None
            new_atom = Atom()
            line = s.split()
            new_atom.symbol = line[0]
            new_atom.coord[0] = float(line[1])
            new_atom.coord[1] = float(line[2])
            new_atom.coord[2] = float(line[3])
            self.geometries[0].atoms.append(new_atom)

    def link103(self):
        extracted_lines = self.file_raw_contents[self.__start_end[0]: self.__start_end[1]]
        if extracted_lines[3].strip() == "Initialization pass.":
            return None

    def link202(self):
        extracted_lines = self.file_raw_contents[self.__start_end[0]: self.__start_end[1]]
        new_molecule = Molecule()
        if extracted_lines[0].strip() == "Symmetry turned off by external request.":
            new_molecule.symmetrized = False
        else:
            new_molecule.symmetrized = True
        a = extracted_lines[4].split()
        new_molecule.point_group = a[1]
        i = 0
        check_line = extracted_lines[i + 10]
        while check_line.strip()[0] != "-":
            a = check_line.split()
            atom_index = int(a[1])
            new_atom = Atom()
            new_atom.symbol = list(element_weight)[atom_index]
            new_atom.assign_weight()
            new_atom.coord[0] = float(a[3])
            new_atom.coord[1] = float(a[4])
            new_atom.coord[2] = float(a[5])
            new_molecule.atoms.append(new_atom)
            i += 1
            check_line = extracted_lines[i + 10]
        self.geometries.append(new_molecule)

    def link502(self):
        extracted_lines = self.file_raw_contents[self.__start_end[0]: self.__start_end[1]]
        new_scf = GaussianFile.SCF()
        for num, line in enumerate(extracted_lines):
            str1 = line.split()
            if str1[0] == "Requested":
                str2 = line.split("=")
                value = log_notation_to_float(str2[1][:9])
                str2split = str2[0].split()
                if str2split[len(str2split) - 1] == "matrix":
                    if str2split[len(str2split) - 3] == "RMS":
                        new_scf.RMS_DP_criteria = value
                    elif str2split[len(str2split) - 3] == "MAX":
                        new_scf.Max_DP_criteria = value
                elif str2split[len(str2split) - 1] == "energy":
                    new_scf.DeltaE_criteria = value
            if new_scf.DeltaE_criteria != 0.0 and new_scf.RMS_DP_criteria != 0.0 and new_scf.Max_DP_criteria != 0.0:
                break
        first_space = 0
        for num, line in enumerate(extracted_lines):
            if line == " ":
                first_space = num
                break
        for num, line in enumerate(extracted_lines):
            line_split = line.split()
            if line_split[0] == "E=":
                new_scf.energy.append(float(line_split[1]))
                break
        for i in range(first_space + 1, len(extracted_lines)):
            if extracted_lines[i] == " ":
                line = extracted_lines[i - 1]
                line_split = line.split()
                for num, l in enumerate(line_split):
                    l_split = l.split("=")
                    if l_split[0] == "RMSDP":
                        new_scf.RMSDP.append(log_notation_to_float(l_split[1]))
                    elif l_split[0] == "MaxDP":
                        new_scf.MaxDP.append(log_notation_to_float(l_split[1]))
                    elif l_split[0] == "DE":
                        new_scf.delta_energy.append(log_notation_to_float(l_split[1]))
                    elif l_split[0] == "OVMax":
                        new_scf.OVMax.append(log_notation_to_float(l_split[1]))
        self.scf_stages.append(new_scf)

    def link601(self):
        extracted_lines = self.file_raw_contents[self.__start_end[0]: self.__start_end[1]]
        mulliken_charges_pos = 0
        for num, line in enumerate(extracted_lines):
            if line.strip() == "Mulliken charges:":
                mulliken_charges_pos = num
                break
            line_split = line.split("--")
            a = line_split[0]
            b = line_split[1]
            a_split = a.split()
            value = 0.0
            if a_split[1] == "occ.":
                b_split = b.split(".")
                for i in range(len(b_split)):
                    if i % 2 == 0:
                        value = float(b_split[i])
                    else:
                        value += float(10 ** (-1 * len(b_split[i]))) * float(b_split[i])
                self.scf_stages[len(self.scf_stages) - 1].occ_eigenvalues.append(value)
            elif a_split[1] == "virt.":
                b_split = b.split(".")
                for i in range(len(b_split)):
                    if i % 2 == 0:
                        value = float(b_split[i])
                    else:
                        value += float(10 ** (-1 * len(b_split[i]))) * float(b_split[i])
                self.scf_stages[len(self.scf_stages) - 1].virt_eigenvalues.append(value)
        for i in range(mulliken_charges_pos + 1, len(extracted_lines) - 1):
            line_split = extracted_lines[i]
            if line_split[0] == "Sum":
                break
            self.geometries[len(self.geometries) - 1].atoms[int(line_split[0])].charge = float(line_split[2])
        for num, line in enumerate(extracted_lines):
            line_split = line.split()
            if line_split[0] == "Dipole":
                dipole_line = extracted_lines[num + 1].split()
                self.geometries[len(self.geometries) - 1].dipole_moment[0] = float(dipole_line[1])
                self.geometries[len(self.geometries) - 1].dipole_moment[1] = float(dipole_line[3])
                self.geometries[len(self.geometries) - 1].dipole_moment[2] = float(dipole_line[5])
            if line_split[0] == "Quadrupole":
                quadrupole_line = extracted_lines[num + 1].split() + extracted_lines[num + 2].split()
                self.geometries[len(self.geometries) - 1].quadrupole_moment[0, 0] = float(quadrupole_line[1])
                self.geometries[len(self.geometries) - 1].quadrupole_moment[1, 1] = float(quadrupole_line[3])
                self.geometries[len(self.geometries) - 1].quadrupole_moment[2, 2] = float(quadrupole_line[5])
                self.geometries[len(self.geometries) - 1].quadrupole_moment[0, 1] = float(quadrupole_line[7])
                self.geometries[len(self.geometries) - 1].quadrupole_moment[0, 2] = float(quadrupole_line[9])
                self.geometries[len(self.geometries) - 1].quadrupole_moment[1, 2] = float(quadrupole_line[11])
                self.geometries[len(self.geometries) - 1].quadrupole_moment[1, 0] = (
                    self.geometries[len(self.geometries) - 1].quadrupole_moment)[0, 1]
                self.geometries[len(self.geometries) - 1].quadrupole_moment[2, 0] = (
                    self.geometries[len(self.geometries) - 1].quadrupole_moment)[0, 2]
                self.geometries[len(self.geometries) - 1].quadrupole_moment[2, 1] = (
                    self.geometries[len(self.geometries) - 1].quadrupole_moment)[2, 1]

    def read(self, file_path=""):
        links_dict = {
            "L1": self.link1(),
            "L101": self.link101(),
            "L103": self.link103(),
            "L202": self.link202(),
            "L502": self.link502(),
            "L601": self.link601(),
        }
        self.path = file_path
        try:
            if "\\" not in file_path:  # Try to open file in the same directory
                file = open(os.path.join(os.getcwd(), file_path), "r")
            else:
                file = open(file_path, "r")
            self.file_raw_contents = file.readlines()
            file.close()
        except OSError:
            kasuga_io.quit_with_error(f'Can`t open: {file_path}')
        self.__start_end = [0, 0]
        for num, line in enumerate(self.file_raw_contents):
            split_line = line.split()
            if "Link 1" in split_line and "Entering" in split_line:
                self.__start_end[0] = num
            if "Link 1" in split_line and "Leaving" in split_line:
                self.__start_end[1] = num
            if self.__start_end != [0, 0]:
                self.link1()
                break
        self.__start_end = [0, 0]
        for num, line in enumerate(self.file_raw_contents):
            split_line = line.split()
            detected_link = ""
            if split_line[0] == "(Enter":
                self.__start_end[0] = num
                link_split1 = split_line[1].split("/")
                link_split2 = link_split1[len(link_split1) - 1].split(".")
                detected_link = link_split2[0]
                detected_link[0] = "L"
            if "Leave" in split_line and "Link" in split_line:
                self.__start_end[1] = num
                if detected_link in links_dict:
                    links_dict[detected_link]
                    self.__start_end = [0, 0]


class GaussianCube:

    def __init__(self):
        self.num_atoms = 0
        self.steps = np.zeros((3, 1), dtype=int)
        self.volume = np.zeros((3, 3))
        self.origin = np.zeros((3))
        self.voxels = np.zeros((3, 3, 3))
        self.grid = []
        self.dv = 0

    def read(self, path=""):
        try:
            file = open(path, "r")
        except OSError:
            print("Could not open the CUBE file at: ", path)
            exit(-1)
        contents = file.readlines()
        file.close()
        words = contents[2].split()
        self.num_atoms = int(words[0])
        self.origin = np.array([float(words[1]), float(words[2]), float(words[3])])
        for i in range(3):
            words = contents[i + 3].split()
            self.steps[i, 0] = int(words[0])
            self.volume[i, 0] = float(words[1])
            self.volume[i, 1] = float(words[2])
            self.volume[i, 2] = float(words[3])
        # self.molecule = Molecule()
        # for i in range(6, 6 + self.num_atoms):
        #     words = contents[i].split()
        #     self.molecule.atom_label.append(dic.nrelements[int(words[0])])
        #     self.molecule.atom_coord[i:i + 1] = np.array([float(words[2]), float(words[3]), float(words[4])])
        self.voxels = np.zeros((self.steps[0, 0], self.steps[1, 0], self.steps[2, 0]))
        sep_contents = []
        for i1 in range(6 + self.num_atoms, len(contents)):
            words = contents[i1].split()
            for i2 in range(len(words)):
                sep_contents.append(float(words[i2]))
        for x in range(self.steps[0, 0]):
            for y in range(self.steps[1, 0]):
                for z in range(self.steps[2, 0]):
                    n = z + y * self.steps[2, 0] + x * self.steps[1, 0] * self.steps[2, 0]
                    self.voxels[x, y, z] = sep_contents[n]
        self.dv = self.volume[0, 0] * self.volume[1, 1] * self.volume[2, 2]
        for x in range(self.steps[0, 0]):
            for y in range(self.steps[1, 0]):
                for z in range(self.steps[2, 0]):
                    temp = self. origin + np.array([x * self.volume[0, 0], y * self.volume[1, 1], z * self.volume[2, 2]])
                    self.grid.append(temp)


class CifFile:
    # Parser and processor for .cif files according to CIF v1.1 standard

    def __init__(self):
        self.transform_matrix = None
        self.tags = {}  # Single fields from CIF
        self.loops = []  # Looped fields from CIF
        # Cell parameters
        self.cell_length_a = float
        self.cell_length_b = float
        self.cell_length_c = float
        self.cell_angle_alpha = float
        self.cell_angle_beta = float
        self.cell_angle_gamma = float
        # Cartesian translation vectors
        self.translation_a = Vector()
        self.translation_b = Vector()
        self.translation_c = Vector()
        # Transformation matrix from abc-system to Cartesian
        self.coord_transform_matrix = np.zeros((3, 3))
        # Asymmetric unit of a primitive cell
        self.as_unit = Molecule()
        self.as_molecules = []
        self.xyz_eq = []

    @staticmethod
    def parse_xyz_eq(eq: list):
        for num, s in enumerate(eq):
            eq[num] = s.strip()
        transformation_matrix = np.zeros((3, 3))
        translation_vector = np.zeros(3)
        c = ["x", "y", "z"]
        for i1 in range(3):
            for i2 in range(3):
                if c[i2] in eq[i1]:
                    if eq[i1][0] == "-":
                        transformation_matrix[i2, i1] = -1
                        eq[i1] = eq[i1][2:]
                    else:
                        transformation_matrix[i2, i1] = 1
                        eq[i1] = eq[i1][1:]
        for i in range(3):
            translation = 0
            if eq[i] != "":
                split = eq[i].split("/")
                translation = float(split[0]) / float(split[1])
            translation_vector[i] = translation
        return translation_vector, transformation_matrix

    @staticmethod
    def parse_line(line):
        line = line.strip()
        if line == "?" or line == ".":
            return ""
        if line[0] == "'" and line[len(line) - 1] == "'":
            return line[1:len(line) - 1]
        split = line.split()
        out = []
        for s in split:
            if s[0] == "'" and s[len(s) - 1] == "'":
                out.append(s[1:len(s) - 1])
                continue
            else:
                try:
                    pre = int(s)
                    out.append(pre)
                except ValueError:
                    try:
                        pre = float(s)
                        out.append(pre)
                    except ValueError:
                        if isinstance(s, str):
                            out.append(s)
                        elif "(" in s and ")" in s:
                            temp = float(s[0:s.find("(")])
                            for i in range(3):
                                l1 = len(s[s.find(".") + 1: s.find("(")])
                                l2 = len(s[s.find("(") + 1: s.find(")")])
                                temp += float(s[s.find("(") + 1:s.find(")")]) * 10 ** (-1 * (l1 + (i + 1) * l2))
                            out.append(temp)
                        else:
                            kasuga_io.quit_with_error(f'Unrecognized variable in "{line}"!')
        if len(out) == 1:
            return out[0]
        else:
            return out

    def read_raw(self, file_path):
        file_contents = []
        parsed_loop = []
        loop_tags = []
        loop_contents = []
        try:
            if "\\" not in file_path:  # Try to open file in the same directory
                file = open(os.path.join(os.getcwd(), file_path), "r")
            else:
                file = open(file_path, "r")
            file_contents = file.readlines()
            file.close()
        except OSError:
            kasuga_io.quit_with_error(f'Can`t open: {file_path}')

        for i in range(len(file_contents)):  # Initial survey for any Shelxl data to be expunged
            if "_shelx_res_file" in file_contents[i]:
                file_contents = file_contents[:i]  # Slice away everything below first Shelxl tag
                break

        in_loop = False
        loop_parsed = False

        for index in range(len(file_contents)):
            if "loop_" in file_contents[index]:
                in_loop = True
                parsed_loop = []
                loop_tags = []
                loop_contents = []
                continue

            if in_loop and loop_parsed:
                if file_contents[index].strip() == "":
                    in_loop = False
                    loop_parsed = False
                    continue
                else:
                    continue

            if in_loop:
                start_index = index
                reading_tags = True

                for i in range(start_index, len(file_contents)):
                    if file_contents[i].strip()[0] == "_":
                        split = file_contents[i].split()
                        loop_tags.append(split[0][1:])
                    else:
                        reading_tags = False
                        start_index = i
                        break

                if not reading_tags:
                    for i in range(start_index, len(file_contents)):
                        if file_contents[i].strip() == "":
                            if len(loop_contents) % len(loop_tags) != 0:
                                kasuga_io.quit_with_error(f'Faulty loop block around "{file_contents[index]}" '
                                                          f'and "{file_contents[i]}"! '
                                                          f'Please verify "{file_path}" integrity.')
                            else:
                                for i1 in range(len(loop_contents) // len(loop_tags)):
                                    d = {}
                                    for i2 in range(len(loop_tags)):
                                        d[loop_tags[i2]] = loop_contents[i1 + i2]
                                    parsed_loop.append(d)
                                self.loops.append(parsed_loop)
                                loop_parsed = True
                                break
                        else:
                            loop_pre_contents = self.parse_line(file_contents[i])
                            if isinstance(loop_pre_contents, list):
                                for i2 in loop_pre_contents:
                                    loop_contents.append(i2)
                            else:
                                loop_contents.append(loop_pre_contents)

            if file_contents[index][0] == "_" and not in_loop:
                split = file_contents[index].split()
                tag_content = ""

                cd_block_encountered = False
                if len(split) > 1:
                    for i in range(1, len(split)):
                        tag_content += split[i]
                elif ";" in file_contents[index + 1]:
                    cd_block_encountered = True
                    ind = index + 2
                    if file_contents[ind] == ";":
                        kasuga_io.quit_with_error(f'Faulty tag ;-; block encountered around "{file_contents[index]}"! '
                                                  f'Please verify "{file_path}" integrity.')
                    else:
                        for i in range(ind, len(file_contents)):
                            if ";" in file_contents[i]:
                                break
                            else:
                                tag_content += file_contents[i]
                else:
                    tag_content = file_contents[index + 1]

                if tag_content == "" or split[0] == "_":
                    kasuga_io.quit_with_error(f'Faulty CIF tag encountered around "{file_contents[index]}"!'
                                              f' Please verify "{file_path}" integrity.')
                else:
                    if cd_block_encountered:
                        self.tags[split[0][1:]] = tag_content
                    else:
                        self.tags[split[0][1:]] = self.parse_line(tag_content)

    def parse_raw(self):
        # Primitive cell dimensions
        self.cell_length_a = self.tags['cell_length_a']
        self.cell_length_b = self.tags['cell_length_b']
        self.cell_length_c = self.tags['cell_length_c']

        # Primitive cell angles
        self.cell_angle_alpha = self.tags['cell_angle_alpha']  # Between c and b
        self.cell_angle_beta = self.tags['cell_angle_beta']  # Between c and a
        self.cell_angle_gamma = self.tags['cell_angle_gamma']  # Between a and b

        # Generate transformation matrix from abc to Cartesian
        cosa = np.cos(np.deg2rad(self.cell_angle_alpha))
        cosb = np.cos(np.deg2rad(self.cell_angle_beta))
        cosg = np.cos(np.deg2rad(self.cell_angle_gamma))
        sing = np.sin(np.deg2rad(self.cell_angle_gamma))
        volume = np.sqrt(1.0 - cosa ** 2.0 - cosb ** 2.0 - cosg ** 2.0 + 2.0 * cosa * cosb * cosg)
        self.transform_matrix = np.array([[self.cell_length_a, self.cell_length_b * cosg, self.cell_length_c * cosb],
                                         [0, self.cell_length_b * sing, self.cell_length_c * (cosa - cosb * cosg) / sing],
                                         [0, 0, self.cell_length_c * volume / sing]])

        # Translation vectors in cartesian coordinates
        # Most of the symmetry operations are performed in the abc-system for the sake of simplicity
        # Yet, for some processing down the line we might need cartesian vectors as well
        self.translation_a.coord[0] = self.cell_length_a  # We assume that X-axis is aligned with a-axis
        self.translation_b.coord = np.matmul(np.array([0, 1, 0]), self.coord_transform_matrix)
        self.translation_c.coord = np.matmul(np.array([0, 0, 1]), self.coord_transform_matrix)

        # Extract fractional coordinates from CIF loops
        found_as = False
        for i1 in range(len(self.loops)):
            if "atom_site_label" in self.loops[i1][0]:
                if found_as:
                    kasuga_io.quit_with_error(f'Duplicated asymmetric units in CIF file!')
                else:
                    found_as = True
                for i2 in range(len(self.loops[i1])):
                    a = Atom()
                    a.symbol = self.loops[i1][i2]['atom_site_type_symbol']
                    a.assign_weight()
                    a.coord[0] = self.loops[i1][i2]['atom_site_fract_x']
                    a.coord[1] = self.loops[i1][i2]['atom_site_fract_y']
                    a.coord[2] = self.loops[i1][i2]['atom_site_fract_z']
                    self.as_unit.atoms.append(a)
        self.xyz_eq = SymOpsHall[HM2Hall[self.tags["_symmetry_space_group_name_H-M"]]]

    def build_as_molecules(self):
        mol_to_add = []
        self.as_unit.rebuild_connectivity()
        self.as_molecules = self.as_unit.separate_molecules()
        for s in self.xyz_eq:
            vector, matrix = self.parse_xyz_eq(s)
            for mol in self.as_molecules:
                new_molecule = Molecule()
                new_molecule.atoms = mol.atoms
                for i in range(len(new_molecule.atoms) - 1):
                    new_molecule.atoms[i] = np.matmul(new_molecule.atoms[i], matrix)
                    new_molecule.atoms[i] += vector
                    if mol != new_molecule:
                        flag = False
                        for i2 in range(len(self.as_molecules) - 1):
                            if new_molecule.is_connected(self.as_molecules[i2]):
                                flag = True
                        if flag:
                            mol_to_add.append(new_molecule)
        for m in mol_to_add:
            self.as_molecules.append(m)
        whole = Molecule()
        for m in self.as_molecules:
            whole += m
        whole.rebuild_connectivity()
        self.as_molecules = whole.separate_molecules()


# Cluster is an array of molecules either natively generated from an associated CifFile or appended through other means
class Cluster:

    def __init__(self):
        self.molecules = []
        self.cif = CifFile()

    def molecule_is_inside_pc(self, cartesian=True, mol=Molecule()):
        mc = mol.get_mass_center()
        if cartesian:
            if mc[0] <= self.cif.cell_length_a and mc[1] <= self.cif.cell_length_b and mc[2] <= self.cif.cell_length_c:
                return True
        else:
            if mc[0] <= 1 and mc[1] <= 1 and mc[2] <= 1:
                return True
        return False

    def molecule_is_in_cluster(self, m: Molecule):
        for mol in self.molecules:
            if m == mol:
                return True
        return False

    def init_cif(self, cif_file: CifFile):
        self.cif = cif_file
        self.cif.build_as_molecules()
        self.molecules = self.cif.as_molecules

    def build_primitive_cell(self):
        new_molecule = Molecule()
        for m in self.molecules:
            for xyz in self.cif.xyz_eq:
                vector, matrix = CifFile.parse_xyz_eq(xyz)
                new_molecule.atoms = m.atoms
                for a in range(len(new_molecule.atoms) - 1):
                    new_molecule.atoms[a] = np.matmul(new_molecule.atoms[a], matrix)
                    new_molecule.atoms[a] += vector
                flag = True
                for m in self.molecules:
                    if new_molecule == m:
                        flag = False
                if flag:
                    for m1 in self.molecules:
                        if new_molecule.is_connected(m1):
                            m1 += new_molecule
                        else:
                            self.molecules.append(new_molecule)
        for_deletion = []
        for m1 in self.molecules:
            for m2 in self.molecules:
                if m1.is_connected(m2):
                    m2 += m1
                    if m1 not in for_deletion:
                        for_deletion.append(m1)
            if not self.molecule_is_inside_pc(m1):
                if m1 not in for_deletion:
                    for_deletion.append(m1)
        for d in for_deletion:
            self.molecules.remove(d)

    def transform_to_cartesian(self):
        for num1, mol in enumerate(self.molecules):
            for num2, a in enumerate(mol.atoms):
                mol.atoms.coord[num2] = np.matmul(a.coord, self.cif.transform_matrix)
            self.molecules[num1].atoms = mol.atoms

    def multiply(self, direction="x", count=1):
        translation_vector =  np.zeros(3)
        xyz = ["x", "y", "z"]
        abc = ["a", "b", "c"]
        if direction in xyz:
            translation_vector[xyz.index(direction)] = 1
        elif direction in xyz:
            translation_vector[abc.index(direction)] = 1
        elif type(direction) is np.ndarray:
            translation_vector = direction
        else:
            kasuga_io.quit_with_error(f'vector in multiply cluster routine: {direction}')
        cloned_molecules = []
        for mol in self.molecules:
            for i in range(count)
                new_molecule = Molecule()
                new_molecule.atoms = mol.atoms
                new_molecule
                if not self.molecule_is_in_cluster(new_molecule):

