import os

class CifFile:

    # Parser and processor for .cif files according to CIF v1.1 standard

    def __init__(self):
        self.tags = {}  # Dict containing all data fields found in .cif file

    def read(self, file_path=""):
        if "/" not in file_path:
            try:
                file = open(os.path.join(os.getcwd(), file_path), "r")
                file_contents = file.readlines()
                file.close()
            except OSError:
                print("Error! Can't open: ", file_path)
                exit(-1)