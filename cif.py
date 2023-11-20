import os
import kasuga_io


def periodic_to_float(number, count=1):
    """
    Round periodic number to regular float with accuracy count (1 by default).
    """
    if "(" not in number and ")" not in number:
        return float(number)
    else:
        temp = float(number[0:number.find("(")])
        for i in range(count):
            l1 = len(number[number.find(".") + 1: number.find("(")])
            l2 = len(number[number.find("(") + 1: number.find(")")])
            temp += float(number[number.find("(") + 1:number.find(")")]) * 10 ** (-1 * (l1 + (i + 1) * l2))
        return temp


class CifFile:

    # Parser and processor for .cif files according to CIF v1.1 standard

    def __init__(self):
        self.tags = {}  # Dict containing all data fields found in .cif file
        self.loops = {}

    def read_raw(self, file_path=""):
        file_contents = []
        loop = []
        loop_tags = []
        loop_contents = []
        try:
            if "/" not in file_path:  # Try to open file in the same directory
                file = open(os.path.join(os.getcwd(), file_path), "r")
            else:
                file = open(file_path, "r")
            file_contents = file.readlines()
            file.close()
        except OSError:
            kasuga_io.quit_with_error(f'Error! Can`t open: {file_path}')

        for i in range(len(file_contents)):  # Initial survey for any Shelxl data to be expunged
            if "_shelx_res_file" in file_contents[i]:
                file_contents = file_contents[:i]  # Slice away everything below first Shelxl tag
                break

        in_loop = False

        for s in file_contents:
            if s == "loop_":
                in_loop = True
                loop = []
                loop_tags = []
                loop_contents = []
                continue
            if s == "" and in_loop:
                in_loop = False
                continue
            if in_loop:
                start_index = file_contents.index(s)
                reading_tags = True
                for i in range(start_index, len(file_contents)):
                    if file_contents[i].strip()[0] == "_":
                        split = file_contents[i].split()
                        loop_tags.append(split[0][1:])
                    else:
                        reading_tags = False
                        start_index = i
                if not reading_tags:
                    for i in range(start_index, len(file_contents)):
                        if file_contents[i] == "":
                            if len(loop_contents) % len(loop_tags) != 0:
                                kasuga_io.quit_with_error(f'Faulty loop block around "{s}"! Please verify "{file_path}" integrity.')
                            else:
                                for i1 in range(len(loop_contents) // len(loop_tags)):
                                    d = {}
                                    for i2 in range(len(loop_tags)):
                                        d[loop_tags[i2]] = loop_contents[i1 + i2]
                                    loop.append(d)
                                break
                        else:
                            loop_pre_contents = file_contents[i].split()
                            for l in loop_pre_contents:
                                loop_contents.append(l.strip())
            if s[0] == "_" and not in_loop:
                split = s.split()
                tag_content = ""
                if len(split) > 1:
                    for i in range(1, len(split)):
                        tag_content += split[i]
                elif file_contents[file_contents.index(s) + 1] == ";":
                    ind = file_contents.index(s) + 2
                    if file_contents[ind] == ";":
                        kasuga_io.quit_with_error(f'Faulty tag ;-; block encountered around "{s}"! Please verify "{file_path}" integrity.')
                    else:
                        for i in range(ind, len(file_contents)):
                            if file_contents[i] == ";":
                                break
                            else:
                                tag_content += file_contents[i]
                else:
                    tag_content = file_contents[file_contents.index(s) + 1].strip()
                if tag_content == "" or split[0] == "_":
                    kasuga_io.quit_with_error(f'Faulty CIF tag encountered around "{s}"! Please verify "{file_path}" integrity.')
                else:
                    self.tags[split[0][1:]] = tag_content
