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
"""

import os
import kasuga_io
import numpy as np

# Atomic weights for each respected atom
element_weight = {
    'H': 1.0075,
    'D': 2.01410178,
    'He': 4.002,
    'Li': 6.9675,
    'Be': 9.012,
    'B': 10.8135,
    'C': 12.0106,
    'N': 14.0065,
    'O': 15.999,
    'F': 18.998,
    'Ne': 20.1797,
    'Na': 22.989,
    'Mg': 24.3050,
    'Al': 26.981,
    'Si': 28.085,
    'P': 30.973,
    'S': 32.0675,
    'Cl': 35.4515,
    'Ar': 39.948,
    'K': 39.0983,
    'Ca': 40.078,
    'Sc': 44.955,
    'Ti': 47.867,
    'V': 50.9415,
    'Cr': 51.9961,
    'Mn': 54.938,
    'Fe': 55.845,
    'Co': 58.933,
    'Ni': 58.6934,
    'Cu': 63.546,
    'Zn': 65.38,
    'Ga': 69.723,
    'Ge': 72.63,
    'As': 74.921,
    'Se': 78.96,
    'Br': 79.904,
    'Kr': 83.798,
    'Rb': 85.4678,
    'Sr': 87.62,
    'Y': 88.905,
    'Zr': 91.224,
    'Nb': 92.906,
    'Mo': 95.96,
    'Tc': 98,
    'Ru': 101.07,
    'Rh': 102.905,
    'Pd': 106.42,
    'Ag': 107.8682,
    'Cd': 112.411,
    'In': 114.818,
    'Sn': 118.710,
    'Sb': 121.760,
    'Te': 127.60,
    'I': 126.904,
    'Xe': 131.293,
    'Cs': 132.905,
    'Ba': 137.327,
    'La': 138.905,
    'Ce': 140.116,
    'Pr': 140.907,
    'Nd': 144.242,
    'Pm': 145,
    'Sm': 150.36,
    'Eu': 151.964,
    'Gd': 157.25,
    'Tb': 158.925,
    'Dy': 162.500,
    'Ho': 164.930,
    'Er': 167.259,
    'Tm': 168.934,
    'Yb': 173.054,
    'Lu': 174.9668,
    'Hf': 178.49,
    'Ta': 180.947,
    'W': 183.84,
    'Re': 186.207,
    'Os': 190.23,
    'Ir': 192.217,
    'Pt': 195.084,
    'Au': 196.966,
    'Hg': 200.59,
    'Tl': 204.3835,
    'Pb': 207.2,
    'Bi': 208.980,
    'Po': 209,
    'At': 210,
    'Rn': 222,
    'Fr': 223,
    'Ra': 226,
    'Ac': 227,
    'Th': 232.038,
    'Pa': 231.035,
    'U': 238.028,
    'Np': 237,
    'Pu': 244,
    'Am': 243,
    'Cm': 247,
    'Bk': 247,
    'Cf': 251,
    'Es': 252,
    'Fm': 257,
    'Md': 258,
    'No': 259,
    'Lr': 262,
    'Rf': 267,
    'Db': 268,
    'Sg': 271,
    'Bh': 272,
    'Hs': 270,
    'Mt': 276,
    'Ds': 281,
    'Rg': 280,
    'Cn': 285
}

covalent_radius = {
    'H': 0.32,
    'D': 0.32,
    'Ne': 0.71,
    'F': 0.72,
    'O': 0.73,
    'N': 0.75,
    'C': 0.77,
    'B': 0.82,
    'Be': 0.90,
    'He': 0.93,
    'Ar': 0.98,
    'Cl': 0.99,
    'S': 1.02,
    'P': 1.06,
    'Si': 1.11,
    'Kr': 1.12,
    'Br': 1.14,
    'Ni': 1.15,
    'Se': 1.16,
    'Co': 1.16,
    'Cu': 1.17,
    'Fe': 1.17,
    'Mn': 1.17,
    'Al': 1.18,
    'Cr': 1.18,
    'As': 1.20,
    'Ge': 1.22,
    'V': 1.22,
    'Li': 1.23,
    'Rh': 1.25,
    'Ru': 1.25,
    'Zn': 1.25,
    'Ga': 1.26,
    'Os': 1.26,
    'Ir': 1.27,
    'Tc': 1.27,
    'Re': 1.28,
    'Pd': 1.28,
    'W': 1.30,
    'Pt': 1.30,
    'Mo': 1.30,
    'Xe': 1.31,
    'Ti': 1.32,
    'I': 1.33,
    'Ta': 1.34,
    'Nb': 1.34,
    'Ag': 1.34,
    'Au': 1.34,
    'Te': 1.36,
    'Mg': 1.36,
    'Sn': 1.41,
    'Sb': 1.41,
    'U': 1.42,
    'In': 1.44,
    'Sc': 1.44,
    'Hf': 1.44,
    'Zr': 1.45,
    'At': 1.45,
    'Bi': 1.46,
    'Po': 1.46,
    'Pb': 1.47,
    'Cd': 1.48,
    'Tl': 1.48,
    'Hg': 1.49,
    'Na': 1.54,
    'Tm': 1.56,
    'Lu': 1.56,
    'Er': 1.57,
    'Ho': 1.58,
    'Dy': 1.59,
    'Tb': 1.59,
    'Gd': 1.61,
    'Y': 1.62,
    'Sm': 1.62,
    'Pm': 1.63,
    'Nd': 1.64,
    'Th': 1.65,
    'Ce': 1.65,
    'Pr': 1.65,
    'La': 1.69,
    'Yb': 1.74,
    'Ca': 1.74,
    'Eu': 1.85,
    'Pu': 1.87,
    'Sr': 1.91,
    'Ba': 1.98,
    'K': 2.03,
    'Rb': 2.16,
    'Cs': 2.35
}


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
        :param normal: normal vector of a mirror plane (np.array, [1,0,0] bu default)
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
        distance = abs(a * self.coord[0] + b * self.coord[1] + c * self.coord[3] + d) / np.sqrt(a**2 + b**2 + c**2)
        # provided normal vector can either face the same direction as mirrored point or the opposite
        test1 = self.coord + 2 * n * distance
        test2 = self.coord - 2 * n * distance
        distance_test1 = abs(a * test1[0] + b * test1[1] + c * test1[3] + d) / np.sqrt(a**2 + b**2 + c**2)
        distance_test2 = abs(a * test2[0] + b * test2[1] + c * test2[3] + d) / np.sqrt(a**2 + b**2 + c**2)
        if distance_test1 < distance_test2:
            # same direction
            self.coord = test1
        else:
            # opposite direction
            self.coord = test2

    def xyz_mirror(self, plane="xy", plane_point=np.zeros(3)):
        """
        Simplified Vector.mirror method to reflect vector in xy, xz and yz planes
        :param plane: string representing one of default planes: "xy", "xz", "yz" or "ab", "ac", "bc" for fractional coordinates.
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


class Atom(Vector):

    def assign_weight(self):
        if self.symbol in element_weight:
            self.weight = element_weight[self.symbol]
        else:
            kasuga_io.quit_with_error(f'Unrecognized {self.symbol} atom encountered!')

    def __init__(self, symbol=""):
        self.weight = float  # Atomic weight
        self.symbol = str  # Chemical symbol of an atom
        super().__init__()
        if symbol != "":
            self.assign_weight()


def atoms_distance(a: Atom, b: Atom, simplified=False):
    diff_vector = a.coord - b.coord
    if simplified:
        # Sometimes, when performing especially heavy calculations, it might be useful to take crude simplified approach
        return abs(diff_vector[1, 1]) + abs(diff_vector[2, 1]) + abs(diff_vector[3, 1])
    else:
        # But usually numpy is fast enough, so it's False by default
        return np.linalg.norm(diff_vector)


def atoms_connected(a: Atom, b: Atom, simplified=False, cutoff=0.025):
    distance = atoms_distance(a, b, simplified)
    covalent_sum = covalent_radius[a.symbol] + covalent_radius[b.symbol]
    if abs(covalent_sum - distance) <= cutoff:
        return True
    else:
        return False


class Molecule:
    atoms = [Atom]
    molecular_formula = ""
    mass_center = np.zeros((3, 1))
    connectivity_matrix = np.zeros((1, 1))
    inertia_tensor = np.zeros((3, 3))
    inertia_eigenvalues = np.zeros((3, 1))
    inertia_eigenvectors = np.zeros((3, 3))
    inertia_x, inertia_y, inertia_z = np.zeros((3, 1))

    def __add__(self, other):
        for i in other.atoms:
            self.atoms.append(i)

    @classmethod
    def get_mass_center(cls):
        if cls.mass_center == np.zeros((3, 1)):
            mass = 0.0
            mass_vector = np.zeros((3, 1))
            for atom in cls.atoms:
                mass += element_weight[atom.symbol]
                mass_vector += atom.coord * element_weight[atom.symbol]
            cls.mass_center = mass_vector / mass
        return cls.mass_center

    @classmethod
    def get_molecular_formula(cls):
        if cls.molecular_formula == "":
            atom_list = [str]
            count_list = [int]
            for i in range(len(cls.atoms)):
                if cls.atoms[i].symbol not in atom_list:
                    atom_list.append(cls.atoms[i].symbol)
                    count_list.append(0)
            atom_list.sort()
            for i1 in range(len(cls.atoms)):
                for i2 in range(len(atom_list)):
                    if cls.atoms[i1].symbol == atom_list[i2]:
                        count_list += 1
            result = ""
            for i in range(len(atom_list)):
                result += (atom_list[i] + str(count_list[i]))
            cls.molecular_formula = result
        return cls.molecular_formula

    @classmethod
    def get_connectivity_matrix(cls):
        if cls.connectivity_matrix == np.zeros((1, 1)):
            cls.connectivity_matrix = np.zeros((len(cls.atoms), len(cls.atoms)))
            for i1 in range(len(cls.atoms)):
                for i2 in range(i1, len(cls.atoms)):
                    if atoms_connected(cls.atoms[i1], cls.atoms[i2]):
                        cls.connectivity_matrix[i1, i2] = 1
                        cls.connectivity_matrix[i2, i1] = 1
        return cls.connectivity_matrix

    @classmethod
    def get_inertia_tensor(cls):
        if cls.inertia_z != np.zeros((3, 1)) and cls.inertia_y != np.zeros((3, 1)) and cls.inertia_x != np.zeros((3, 1)):
            # First, we translate origin to mass center
            mass_center = cls.get_mass_center()
            atoms_mc_system = [Atom]
            for a in cls.atoms:
                new_a = Atom
                new_a.coord = a.coord - mass_center
                new_a.symbol = a.symbol
                atoms_mc_system.append(new_a)
            # Calculate inertia tensor
            for a in atoms_mc_system:
                cls.inertia_tensor[0, 0] += element_weight[a.symbol] * (a.coord[1, 0]**2 + a.coord[2, 0]**2)  # xx
                cls.inertia_tensor[0, 1] += -1 * element_weight[a.symbol] * a.coord[0, 0] * a.coord[1, 0]  # xy
                cls.inertia_tensor[0, 2] += -1 * element_weight[a.symbol] * a.coord[0, 0] * a.coord[2, 0]  # xz
                cls.inertia_tensor[1, 1] += element_weight[a.symbol] * (a.coord[0, 0]**2 + a.coord[2, 0]**2)  # yy
                cls.inertia_tensor[1, 2] += -1 * element_weight[a.symbol] * a.coord[1, 0] * a.coord[2, 0]  # yz
                cls.inertia_tensor[2, 2] += element_weight[a.symbol] * (a.coord[0, 0]**2 + a.coord[1, 0]**2)  # zz
            cls.inertia_tensor[1, 0] = cls.inertia_tensor[0, 1]
            cls.inertia_tensor[2, 0] = cls.inertia_tensor[0, 2]
            cls.inertia_tensor[2, 1] = cls.inertia_tensor[1, 2]
            # Calculate eigenvalues and eigenvectors of the inertia tensor
            cls.inertia_eigenvalues, cls.inertia_eigenvectors = np.linalg.eig(cls.inertia_tensor)
            # Assign eigenvectors to Cartesian axis: highest for Z, lowest for X
            internal_e = cls.inertia_eigenvalues
            for i in range(3):
                index = np.where(internal_e == internal_e.max())
                if cls.inertia_z != np.zeros((3, 1)):
                    cls.inertia_z = (cls.inertia_eigenvalues[index[0]-1:index[0], :] /
                                     np.linalg.norm(cls.inertia_eigenvalues[index[0]-1:index[0], :]))
                elif cls.inertia_y != np.zeros((3, 1)):
                    cls.inertia_y = (cls.inertia_eigenvalues[index[0] - 1:index[0], :] /
                                     np.linalg.norm(cls.inertia_eigenvalues[index[0] - 1:index[0], :]))
                elif cls.inertia_x != np.zeros((3, 1)):
                    cls.inertia_x = (cls.inertia_eigenvalues[index[0] - 1:index[0], :] /
                                     np.linalg.norm(cls.inertia_eigenvalues[index[0] - 1:index[0], :]))
                internal_e[index[0], 0] = -1.0
        return cls.inertia_z, cls.inertia_y, cls.inertia_x


def molecules_equal(a: Molecule, b: Molecule, same_order=True, simplified=False):
    diff = 0.0
    # Simple tests first to potentially save the hassle
    if len(a.atoms) != len(b.atoms):
        return False
    if a.get_molecular_formula() != b.get_molecular_formula():
        return False
    if not same_order:
        # For symmetry cloned molecules it's safe to assume that the order of atoms is still the same
        # But generally it's not always the case, especially if molecules originate from different sources
        for i1 in range(len(a.atoms)):
            current_min = 1000
            for i2 in range(len(b.atoms)):
                # We look for the closest atom with the same symbol
                current_diff = atoms_distance(a.atoms[i1], b.atoms[i2], simplified)
                if current_diff < current_min:
                    current_min = current_diff
            diff += current_min
    else:
        for i1 in range(len(a.atoms)):
            diff = atoms_distance(a.atoms[i1], b.atoms[i1], simplified)
    if diff < 0.05:
        return True
    else:
        return False


def molecules_connected(a: Molecule, b: Molecule, simplified=False):
    for i1 in a.atoms:
        for i2 in b.atoms:
            if atoms_connected(i1, i2, simplified):
                return True
    return False


def molecules_match_rotation(rotated: Molecule, static: Molecule):
    # Extract principal axes for each molecule
    # Static stays in place, rotated is transformed
    rotated_x, rotated_y, rotated_z = rotated.get_inertia_tensor()
    static_x, static_y, static_z = static.get_inertia_tensor()
    # Since the vectors are stored in ((3,1)) shape we need transpose them first
    rotated_x = np.transpose(rotated_x)
    rotated_y = np.transpose(rotated_y)
    rotated_z = np.transpose(rotated_z)
    static_x = np.transpose(static_x)
    static_y = np.transpose(static_y)
    static_z = np.transpose(static_z)
    # Create matrices that rotate standard coordinate system 0 to molecule's principal axes
    rot_mat_rotated = np.array([[rotated_x],[rotated_y],[rotated_z]])
    rot_mat_static = np.array([[static_x], [static_y], [static_z]])
    # Combine two rotations: Rotated -> 0  (inverted 0 -> Rotated) and 0 -> Static
    final_rotation = np.matmul(np.linalg.inv(rot_mat_rotated), rot_mat_static)
    # We translate Rotated to 0 system, perform rotation, and translate it back
    mass_center = rotated.get_mass_center()
    for i in range(len(rotated.atoms)):
        rotated.atoms[i].coord -= mass_center
        rotated.atoms[i].coord = np.transpose(np.matmul(np.transpose(rotated.atoms[i].coord), final_rotation))
        rotated.atoms[i].coord += mass_center


class CifFile:
    # Parser and processor for .cif files according to CIF v1.1 standard
    tags = {}  # Single fields from CIF
    loops = []  # Looped fields from CIF
    # Cell parameters
    cell_length_a = float
    cell_length_b = float
    cell_length_c = float
    cell_angle_alpha = float
    cell_angle_beta = float
    cell_angle_gamma = float
    # Cartesian translation vectors
    translation_a = np.zeros((3, 1))
    translation_b = np.zeros((3, 1))
    translation_c = np.zeros((3, 1))
    # Transformation matrix from abc-system to Cartesian
    transform_matrix = np.zeros((3, 3))
    # Asymmetric unit of a primitive cell
    as_unit = Molecule

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

    @staticmethod
    def transform_abc_to_cartesian(vector: np.array, transform: np.array):
        if transform.shape != (3, 3):
            kasuga_io.quit_with_error(f'Wrong dimensions of a transformation matrix!')
        if vector.shape == (3, 1):
            return np.matmul(transform, vector)
        elif vector.shape == (1, 3):
            return np.transpose(transform, np.transpose(vector))
        else:
            kasuga_io.quit_with_error(f'Wrong dimensions of a transformed vector!')

    @staticmethod
    def transform_cartesian_to_abc(vector: np.array, transform: np.array):
        if transform.shape != (3, 3):
            kasuga_io.quit_with_error(f'Wrong dimensions of a transformation matrix!')
        transform_revert = np.linalg.inv(transform)
        if vector.shape == (3, 1):
            return np.matmul(transform_revert, vector)
        elif vector.shape == (1, 3):
            return np.transpose(transform_revert, np.transpose(vector))
        else:
            kasuga_io.quit_with_error(f'Wrong dimensions of a transformed vector!')

    @classmethod
    def read_raw(cls, file_path):
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
                                cls.loops.append(parsed_loop)
                                loop_parsed = True
                                break
                        else:
                            loop_pre_contents = cls.parse_line(file_contents[i])
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
                        cls.tags[split[0][1:]] = tag_content
                    else:
                        cls.tags[split[0][1:]] = cls.parse_line(tag_content)

    @classmethod
    def parse_raw(cls):
        # Primitive cell dimensions
        cls.cell_length_a = cls.tags['cell_length_a']
        cls.cell_length_b = cls.tags['cell_length_b']
        cls.cell_length_c = cls.tags['cell_length_c']

        # Primitive cell angles
        cls.cell_angle_alpha = cls.tags['cell_angle_alpha']  # Between c and b
        cls.cell_angle_beta = cls.tags['cell_angle_beta']  # Between c and a
        cls.cell_angle_gamma = cls.tags['cell_angle_gamma']  # Between a and b

        # Generate transformation matrix from abc to Cartesian
        cosa = np.cos(np.deg2rad(cls.cell_angle_alpha))
        cosb = np.cos(np.deg2rad(cls.cell_angle_beta))
        cosg = np.cos(np.deg2rad(cls.cell_angle_gamma))
        sing = np.sin(np.deg2rad(cls.cell_angle_gamma))
        volume = np.sqrt(1.0 - cosa ** 2.0 - cosb ** 2.0 - cosg ** 2.0 + 2.0 * cosa * cosb * cosg)
        cls.transform_matrix = np.array([[cls.cell_length_a, cls.cell_length_b * cosg, cls.cell_length_c * cosb],
                                         [0, cls.cell_length_b * sing, cls.cell_length_c * (cosa - cosb * cosg) / sing],
                                         [0, 0, cls.cell_length_c * volume / sing]])

        # Translation vectors in cartesian coordinates
        # Most of the symmetry operations are performed in the abc-system for the sake of simplicity
        # Yet, for some processing down the line we might need cartesian vectors as well
        cls.translation_a[0, 0] = cls.cell_length_a  # We assume that X-axis is aligned with a-axis
        cls.translation_b = CifFile.transform_abc_to_cartesian(np.array([0, 1, 0]), cls.transform_matrix)
        cls.translation_c = CifFile.transform_abc_to_cartesian(np.array([0, 0, 1]), cls.transform_matrix)

        # Extract fractional coordinates from CIF loops
        found_as = False
        for i1 in range(len(cls.loops)):
            if "atom_site_label" in cls.loops[i1][0]:
                if found_as:
                    kasuga_io.quit_with_error(f'Duplicated asymmetric units in CIF file!')
                else:
                    found_as = True
                for i2 in range(len(cls.loops[i1])):
                    a = Atom
                    a.symbol = cls.loops[i1][i2]['atom_site_type_symbol']
                    a.assign_weight()
                    a.coord_abc[0, 1] = cls.loops[i1][i2]['atom_site_fract_x']
                    a.coord_abc[1, 1] = cls.loops[i1][i2]['atom_site_fract_y']
                    a.coord_abc[2, 1] = cls.loops[i1][i2]['atom_site_fract_z']
                    cls.as_unit.atoms.append(a)


# Cluster is an array of molecules either natively generated from an associated CifFile or appended through other means
class Cluster:
    molecules = [Molecule]
    cif = CifFile


