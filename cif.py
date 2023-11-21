import os
import kasuga_io


class CifFile:
    # Parser and processor for .cif files according to CIF v1.1 standard
    tags = {}
    loops = []

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
