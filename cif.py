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


class Atom:
    weight = 0.0  # Atomic weight
    symbol = ""  # Chemical symbol of an atom
    coord = np.zeros((3, 1))  # Cartesian coordinates
    coord_abc = np.zeros((3, 1))  # Coordinates in the abc-system

    @classmethod
    def assign_weights(cls):
        if cls.symbol in element_weight:
            cls.weight = element_weight[cls.symbol]
        else:
            kasuga_io.quit_with_error(f'Unrecognized {cls.symbol} atom encountered!')


class CifFile:
    # Parser and processor for .cif files according to CIF v1.1 standard
    tags = {}  # Single fields from CIF
    loops = []  # Looped fields from CIF
    # Cell parameters
    cell_length_a = 0.0
    cell_length_b = 0.0
    cell_length_c = 0.0
    cell_angle_alpha = 0.0
    cell_angle_beta = 0.0
    cell_angle_gamma = 0.0
    # Cartesian translation vectors
    translation_a = np.zeros((3, 1))
    translation_b = np.zeros((3, 1))
    translation_c = np.zeros((3, 1))
    # Transformation matrix from abc-system to Cartesian
    transform_matrix = np.zeros((3, 3))

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
        cls.translation_a[1, 1] = cls.cell_length_a  # We assume that X-axis is aligned with a-axis
        cls.translation_b = CifFile.transform_abc_to_cartesian(np.array([0, 1, 0]), cls.transform_matrix)
        cls.translation_c = CifFile.transform_abc_to_cartesian(np.array([0, 0, 1]), cls.transform_matrix)