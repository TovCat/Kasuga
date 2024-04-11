import os


def quit_with_error(error_message):
    file = open("./error_output.out", "w")
    print(error_message)
    file.write(error_message)
    exit(-1)


def try_read(file_path):
    file_contents = []
    try:
        if "\\" not in file_path:  # Try to open file in the same directory
            file = open(os.path.join(os.getcwd(), file_path), "r")
        else:
            file = open(file_path, "r")
        file_contents = file.readlines()
        file.close()
    except OSError:
        quit_with_error(f'Can`t open: {file_path}')
    return file_contents
