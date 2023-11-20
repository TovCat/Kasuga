def quit_with_error(error_message):
    file = open("./error_output.out", "w")
    print(error_message)
    file.write(error_message)
    exit(-1)