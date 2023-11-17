import os


class CifFile:

    # Parser and processor for .cif files according to CIF v1.1 standard

    def __init__(self):
        self.tags = {}  # Dict containing all data fields found in .cif file

    def read_raw(self, file_path=""):
        try:
            if "/" not in file_path:  # Try to open file in the same directory
                file = open(os.path.join(os.getcwd(), file_path), "r")
            else:
                file = open(file_path, "r")
            file_contents = file.readlines()
            file.close()
        except OSError:  # TODO: Native IO Error handing
            print("Error! Can't open: ", file_path)
            exit(-1)

        for i in range(len(file_contents) - 1):  # Initial survey for any Shelxl data to be expunged
            if "_shelx_res_file" in file_contents[i]:
                file_contents = file_contents[:i]  # Slice away everything below first Shelxl tag
                break

        for s in file_contents:
            if s[0] == "_":
                split = s.split()
                if len(split) > 1:
                    tag_content = ""
                    for i in range(1, len(split) - 1):
                        tag_content += split[i]
                    self.tags[split[0]] = tag_content
                else:
                    self.tags[split[0]] = file_contents[file_contents.index(s) + 1].strip()
                    
