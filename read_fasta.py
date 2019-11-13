def read_file(file_name):
    sequences = []
    with open(file_name) as f:
        lines = f.readlines()
    for line in lines:
        if line[0] != '>':
            sequences.append(line.rstrip())
    return sequences